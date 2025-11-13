#CPP = g++
CPP = clang++

CFLAGS = -std=c++20 -O3

SRC := q214203-gluing-tetrahedra-together.cpp

LFLAGS = -lgmp -lpthread
#LFLAGS = $(LFLAGS) -lboost_system

OBJ	= $(SRC:%.cpp=%.o)

BIN	= q214203-gluing-tetrahedra-together

all: $(BIN)

.cpp.o:
	$(CPP) $(CFLAGS) -c $< -o $@

$(BIN):\
	$(OBJ) 
	$(CPP) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ): config.h

clean:; rm -f $(OBJ) $(BIN) core
