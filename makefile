# Makefile for Barrier Option Pricing
CXX = g++
CXXFLAGS = -O3 -Wall -std=c++11 -Iinclude
TARGET := barrier_option
SRC := main.cpp $(wildcard ./src/*.cpp)
OBJS := $(SRC:.cpp=.o)

#LDFLAGS

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)


#run: $(TARGET)
#    ./$(TARGET)

#.PHONY: all clean run
