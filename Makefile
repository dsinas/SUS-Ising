EXEC := susIsing

SDIR := src
SRC := $(SDIR)/main.cpp

OBJ := $(SRC:.cpp=.o)

CXX := g++
CXXFLAGS := -O2 -Wall

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@

$(SDIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(SDIR)/*.o 

clean_data:
	rm -rf data
