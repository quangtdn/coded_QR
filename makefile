EXECS=house 2d_house
MPICC?=mpicc

all: ${EXECS}

house: house.c
	${MPICC} -o house house.c -lm

2d_house: 2d_house.c
	${MPICC} -o 2d_house 2d_house.c -lm

clean:
	rm -f ${EXECS}
