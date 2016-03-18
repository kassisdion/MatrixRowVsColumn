NAME= matrixRowVSColumn

CXX= g++

RM= rm -f

SRC= matrixRowVSColumn.cpp

OBJS= $(SRC:.cpp=.o)

CXXFLAGS= -std=c++11 -g -Wall

LDFLAGS= -lpthread

all: $(NAME)

$(NAME): $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS)

clean:
	$(RM) $(OBJS)

fclean: clean
	$(RM) $(NAME)

re: fclean all
