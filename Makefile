# Copyright (c) 2011 the authors listed at the following URL, and/or
# the authors of referenced articles or incorporated external code:
# http://en.literateprograms.org/Hash_table_(C)?action=history&offset=20100620072342
# 
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
# 
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# 
# Retrieved from: http://en.literateprograms.org/Hash_table_(C)?oldid=16749

#find and replace all 'condensed' with 'graph' to run full graph

all: main

main: hashtbl.o main.o cluster.o condensed.o
	cc -o main -Wall -pedantic -g -lm hashtbl.o main.o cluster.o condensed.o

hashtbl.o: hashtbl.c hashtbl.h
	cc -o hashtbl.o -Wall -pedantic -g -c hashtbl.c

main.o: main.c hashtbl.h cluster.h 
	cc -o main.o -Wall -pedantic -g -c main.c

cluster.o: cluster.c hashtbl.h cluster.h
	cc -o cluster.o -Wall -pedantic -g -c cluster.c

condensed.o: condensed.c hashtbl.h cluster.h
	cc -o condensed.o -Wall -pedantic -g -c condensed.c
