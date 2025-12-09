## Executable name
TARGET = generic_app

## Folders
SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build
BIN_DIR = bin

## Compiler and flags
CC = gcc # Compiler
CFLAGS = -O2 -I$(INCLUDE_DIR) # Compilation flags (-Wall)?
LIBS = -lnetcdf -lgsl -lm # Librarys

## Arquivos
SRC = $(wildcard $(SRC_DIR)/*.c) # Source files
OBJ = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRC)) # Generated objects

## Main rules
all: $(BIN_DIR)/$(TARGET)

## Linking
$(BIN_DIR)/$(TARGET): $(OBJ)
	mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

## Compiling .c files to .o
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

## Cleaning the generated files
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

## Rule for recompiling from scratch
rebuild: clean all