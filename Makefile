# Makefile for BAM Statistics Tool

# Compiler and flags
CXX = g++
# Aggressive optimization: -O3 enables auto-vectorization, -march=native for CPU-specific optimizations
CXXFLAGS = -Wall -O3 -std=c++17 -march=native -mtune=native

ifdef DEBUG
  CXXFLAGS += -DDEBUG=true -g
endif

# Detect OS
OS := $(shell uname)

# OS-specific settings
ifeq ($(OS),Darwin)
  # macOS
  CPPFLAGS = -I/usr/local/include -I./ArgumentsParser
  LDFLAGS = -L/usr/local/lib
else
  # Linux/Unix
  HTSLIB = $(shell pwd)/htslib
  CPPFLAGS = -I$(HTSLIB)/include -I./ArgumentsParser
  LDFLAGS = -static -L$(HTSLIB)/lib
endif

# Libraries
LIBS = -lhts -ldeflate -llzma -lbz2 -lz -lm -lpthread

# Target executable
TARGET = bamstats
BIN_DIR = bin
OBJ_DIR = obj
INSTALL_PATH = $(BIN_DIR)/$(TARGET)

# Source files
SRC_DIR = src
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))

# Default target
all: $(INSTALL_PATH)

# Build target
$(INSTALL_PATH): $(OBJS) | $(BIN_DIR)
	@echo "Building BAM Statistics Tool for $(OS)"
ifeq ($(OS),Darwin)
	$(CXX) $(OBJS) $(LDFLAGS) -o $@ $(LIBS)
else
	$(CXX) -static $(OBJS) $(LDFLAGS) -o $@ $(LIBS)
endif
	strip $@
	@echo "Build complete: $@"

# Create bin directory
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Create obj directory
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Compile source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -f $(OBJS)
	rm -rf $(BIN_DIR) $(OBJ_DIR)

# Install (optional)
install: $(INSTALL_PATH)
	cp $(INSTALL_PATH) /usr/local/bin/

.PHONY: all clean install
