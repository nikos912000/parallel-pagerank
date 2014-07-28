CC = gcc
FLG = -O4
NAME = pagerank_pthreads

pagerank_pthreads: pagerank_pthreads.c pagerank_pthreads.h

	$(CC) pagerank_pthreads.c -lpthread -lm -o $(NAME)

clean:
	rm -f *.o *.out *.exe
	rm -f *.bin