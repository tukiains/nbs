FC = gfortran
FCFLAGS = -Wall -Wconversion -Wextra -pedantic -fbounds-check
LFLAGS = 
PROGRAMS = test_n

all: $(PROGRAMS)

test_n.o: math.o body.o velocity_verlet.o

test_n: math.o body.o velocity_verlet.o

%: %.o  
	$(FC) $(FCFLAGS) -o $@ $^ 

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod test_n *~

