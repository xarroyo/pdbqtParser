
all: pdbqtParser pdbqtParser.pg

#-fomit-frame-pointer 
#-static
CC_FLAGS=-O3 -mtune=native -fstrict-aliasing -fno-schedule-insns -ffast-math -Wall -Wextra -lm
OMP=-fopenmp
DEBUG=-DDEBUG -lrt
PATH_TEST=../files
PATH_BIN=../bin

pdbqtParser: pdbqtParser.c
	gcc $(CC_FLAGS) -o pdbqtParser pdbqtParser.c

pdbqtParser.pg: pdbqtParser.c
	gcc $(CC_FLAGS) -DFLOAT -o pdbqtParser.pg pdbqtParser.c $(DEBUG);

pdbqtParser.omp: pdbqtParser.c
	gcc $(CC_FLAGS) -DFLOAT -o pdbqtParser.omp pdbqtParser.c $(DEBUG);

install:
	install -m 644 pdbqtParser $(PATH_BIN);
	install -m 644 pdbqtParser.pg $(PATH_BIN);
	install -m 644 pdbqtParser.omp $(PATH_BIN);

test:
	./pdbqtParser $(PATH_TEST)/target1.pdbqt $(PATH_TEST)/target2.pdbqt

clean:
	rm pdbqtParser pdbqtParser.pg pdbqtParser.omp
