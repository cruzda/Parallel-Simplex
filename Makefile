CC=gcc -g -fopenmp
LIBS=-lm
OBJECTS=amoeba.o simplex.o
SRCS=amoeba.c test.c simplex.c

all: test tags

test: $(OBJECTS) test.o
	$(CC) $(OBJECTS) test.o -o test $(LIBS)

MGH17: $(OBJECTS) MGH17.o
	$(CC) $(OBJECTS) MGH17.o -o MGH17 $(LIBS)
  
lanczos1: $(OBJECTS) lanczos1.o
	$(CC) $(OBJECTS) lanczos1.o -o lanczos1 $(LIBS)
lanczos2: $(OBJECTS) lanczos2.o
	$(CC) $(OBJECTS) lanczos2.o -o lanczos2 $(LIBS)

misra1a: $(OBJECTS) misra1a.o
	$(CC) $(OBJECTS) misra1a.o -o misra1a $(LIBS)

misra1b: $(OBJECTS) misra1b.o
	$(CC) $(OBJECTS) misra1b.o -o misra1b $(LIBS)

chwirut2: $(OBJECTS) chwirut2.o
	$(CC) $(OBJECTS) chwirut2.o -o chwirut2 $(LIBS)

chwirut1: $(OBJECTS) chwirut1.o
	$(CC) $(OBJECTS) chwirut1.o -o chwirut1 $(LIBS)

lanczos3: $(OBJECTS) lanczos3.o
	$(CC) $(OBJECTS) lanczos3.o -o lanczos3 $(LIBS)

gauss1: $(OBJECTS) gauss1.o
	$(CC) $(OBJECTS) gauss1.o -o gauss1 $(LIBS)

gauss2: $(OBJECTS) gauss2.o
	$(CC) $(OBJECTS) gauss2.o -o gauss2 $(LIBS)

danwood: $(OBJECTS) danwood.o
	$(CC) $(OBJECTS) danwood.o -o danwood $(LIBS)

kirby2: $(OBJECTS) kirby2.o
	$(CC) $(OBJECTS) kirby2.o -o kirby2 $(LIBS)

hahn1: $(OBJECTS) hahn1.o
	$(CC) $(OBJECTS) hahn1.o -o hahn1 $(LIBS)

amoeba.o          : amoeba.c  simplex.h
simplex.o         : simplex.c simplex.h
test.o            : test.c    simplex.h
test_momentum.o   : test_momentum.c simplex.h
cost_functions.o  : cost_functions.c simplex.h

tags     : $(SRCS)
	ctags $(SRCS)

clean    :
	rm *.o tags core*
