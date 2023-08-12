all: testing_submit

#For debugging
#OPT=-g -Wall
#For optimistaion
OPT=-O

#All objects (except main) come from cpp and hpp 
%.o:	%.cpp %.hpp
	g++ ${OPT} -c -o $@ $<
#use_vectors relies on objects which rely on headers
testing_submit:	testing_submit.cpp Vector.o Matrix.o Exception.o
		g++ ${OPT} -o testing_submit testing_submit.cpp Vector.o Matrix.o Exception.o

clean:
	rm -f *.o *~ testing_submit

	

	
	