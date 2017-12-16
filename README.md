Parallel Pagerank using Pthreads
================================

A serial and a parallel implementation of Pagerank algorithm, in C, using Pthreads.

4th Course Assignment for Parallel and Distributed Computing Systems (2013).

How to use
----------
1. Download a graph to use (eg. web-Google.txt from [here][google-graph]). You might need to remove any headers from the file.
2. Run `make` command in a unix-based system
3. Run the executable given the appropriate arguments

##### Arguments for serial version
* graph filename
* nodes
* threshold
* d  

Example: 
```
./pagerank_serial web-Google.txt 916428 0.0001 0.85
```

##### Arguments for PTHREADS version
* graph filename
* nodes
* threshold
* d  
* threads 

Example: 
```
./pagerank_pthreads web-Google.txt 916428 0.0001 0.85 8
```

Output
------
Info messages in stdout.

[google-graph]: https://snap.stanford.edu/data/web-Google.html
