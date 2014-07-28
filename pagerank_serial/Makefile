CC = gcc
FLG = -O4
NAME = pagerank_serial

pagerank_serial: pagerank_serial.c pagerank_serial.h

	$(CC) pagerank_serial.c -lm -o $(NAME)

clean:
	rm -f *.o *.out *.exe
	rm -f *.bin