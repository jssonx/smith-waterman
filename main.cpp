#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <chrono>
#include <immintrin.h> // AVX2
#include <climits>     // INT_MIN

#define MATCH_SCORE 3
#define MISMATCH_SCORE -3
#define GAP_PENALTY -2

std::string read_sequence_from_file(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Unable to open file " << filename << std::endl;
        return "";
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

int match_score(char a, char b)
{
    return (a == b) ? MATCH_SCORE : MISMATCH_SCORE;
}

int smith_waterman(const std::string &seq1, const std::string &seq2)
{
    int len1 = seq1.length();
    int len2 = seq2.length();
    std::vector<std::vector<int>> H(len1 + 1, std::vector<int>(len2 + 1, 0));

    int max_score = 0;
    int max_i = 0, max_j = 0;

    for (int i = 1; i <= len1; i++)
    {
        for (int j = 1; j <= len2; j++)
        {
            int score = std::max({H[i - 1][j - 1] + match_score(seq1[i - 1], seq2[j - 1]),
                                  H[i - 1][j] + GAP_PENALTY,
                                  H[i][j - 1] + GAP_PENALTY,
                                  0});

            H[i][j] = score;

            if (score > max_score)
            {
                max_score = score;
                max_i = i;
                max_j = j;
            }
        }
    }
    // print max score
    std::cout << "\nMaximum score: " << max_score << std::endl;
    return max_score;
}

int smith_waterman_avx(const std::string &seq1, const std::string &seq2)
{
    int len1 = seq1.length();
    int len2 = seq2.length();
    std::vector<int> prev(len2 + 1, 0);
    std::vector<int> current(len2 + 1, 0);
    std::vector<int> temp(len2 + 1, 0);

    int max_score = INT_MIN;

    for (int i = 1; i <= len1; i++)
    {
        for (int j = 1; j <= len2; j += 8)
        {
            if (j + 7 <= len2)
            {
                // Load previous scores
                __m256i diag_score = _mm256_loadu_si256((__m256i *)&prev[j - 1]);
                __m256i top_score = _mm256_loadu_si256((__m256i *)&prev[j]);
                __m256i left_score = _mm256_loadu_si256((__m256i *)&current[j - 1]);

                // Compute match/mismatch scores
                int scores[8];

                for (int k = 0; k < 8; ++k)
                {
                    scores[k] = (seq1[i - 1] == seq2[j + k - 1]) ? MATCH_SCORE : MISMATCH_SCORE;
                }

                __m256i match_scores = _mm256_loadu_si256((__m256i *)scores);
                __m256i gap_penalty = _mm256_set1_epi32(GAP_PENALTY);

                // Calculate scores
                __m256i score_diag = _mm256_add_epi32(diag_score, match_scores);
                __m256i score_top = _mm256_add_epi32(top_score, gap_penalty);
                __m256i score_left = _mm256_add_epi32(left_score, gap_penalty);
                __m256i zero = _mm256_setzero_si256();

                // Use max of scores
                __m256i score = _mm256_max_epi32(zero, _mm256_max_epi32(score_diag, _mm256_max_epi32(score_top, score_left)));

                // Store the result
                _mm256_storeu_si256((__m256i *)&current[j], score);

                // Update max score
                for (int k = 0; k < 8; ++k)
                {
                    max_score = std::max(max_score, current[j + k]);
                }
            }
            else
            {
                for (int k = 0; k < len2 - j + 1; ++k)
                {
                    int score = std::max({0,
                                          prev[j + k - 1] + ((seq1[i - 1] == seq2[j + k - 1]) ? MATCH_SCORE : MISMATCH_SCORE),
                                          prev[j + k] + GAP_PENALTY,
                                          current[j + k - 1] + GAP_PENALTY});
                    current[j + k] = score;
                    max_score = std::max(max_score, score);
                }
            }
        }
        std::swap(prev, current);
    }

    std::cout << "\nMaximum score: " << max_score << std::endl;
    return max_score;
}

int main()
{
    std::string seq1 = read_sequence_from_file("./data/GQ457487.txt");
    std::string seq2 = read_sequence_from_file("./data/GQ117044.txt");

    if (seq1.empty() || seq2.empty())
    {
        std::cerr << "Failed to load sequences from files." << std::endl;
        return 1;
    }

    // Output file
    std::ofstream output_file("cpp_output.txt");

    // Naive implementation
    auto start = std::chrono::high_resolution_clock::now();
    int max_score_naive = smith_waterman(seq1, seq2);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> run_time = end - start;
    output_file << "Naive SW score: " << max_score_naive << std::endl;
    std::cout << "Execution time: " << run_time.count() << " milliseconds" << std::endl;

    // Naive + AVX2
    start = std::chrono::high_resolution_clock::now();
    int max_score_avx = smith_waterman_avx(seq1, seq2);
    end = std::chrono::high_resolution_clock::now();
    run_time = end - start;
    output_file << "AVX2 SW score: " << max_score_avx << std::endl;
    std::cout << "Execution time: " << run_time.count() << " milliseconds" << std::endl;

    return 0;
}