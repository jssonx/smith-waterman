CC=g++
CFLAGS=-mavx2 -O3
TARGET=main

all: $(TARGET)

$(TARGET): main.cpp
	$(CC) $(CFLAGS) main.cpp -o $(TARGET)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)

.PHONY: all run clean
