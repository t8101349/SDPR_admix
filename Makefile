export CONDA_PREFIX=/opt/conda/envs/test
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

CC = g++
CFLAGS = -g3 -O3 -ftree-vectorize -DHAVE_INLINE -march=native -fopenmp -std=c++11 -Igsl/include -I$(CONDA_PREFIX)/include
LDFLAGS = -Lgsl/lib -L$(CONDA_PREFIX)/lib -lgsl -lgslcblas -lgomp -lz

all: SDPR_admix regress_prog score

SDPR_admix: parse_gen.o mcmc.o main.o
	${CC} ${CFLAGS} parse_gen.o mcmc.o main.o ${LDFLAGS} -o SDPR_admix

regress_prog: regress.o parse_gen.o mcmc.o
	${CC} ${CFLAGS} regress.o parse_gen.o mcmc.o ${LDFLAGS} -o regress_prog

parse_gen.o: parse_gen.cpp parse_gen.h
	${CC} ${CFLAGS} -c parse_gen.cpp

mcmc.o: mcmc.cpp mcmc.h parse_gen.h regress.h
	${CC} ${CFLAGS} -c mcmc.cpp

main.o: main.cpp mcmc.h parse_gen.h
	${CC} ${CFLAGS} -c main.cpp

regress.o: regress.cpp regress.h
	${CC} ${CFLAGS} -c regress.cpp

score: score.o
	${CC} ${CFLAGS} score.o /opt/conda/envs/test/lib/libz.so -o score

score.o: score.cpp score.h
	${CC} ${CFLAGS} -c score.cpp

clean:
	rm -f *.o SDPR_admix regress_prog score
