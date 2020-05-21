CC=gcc
CFLAGS= -O2 -lm
DEPS=
OBJS=main.o
OUT=main

$(OUT): $(OBJS)
	$(CC)  -o $@  $^ $(CFLAGS)
	
.c.o:
	$(CC) -c -o $@ $< $(CFLAGS) 
	
.PHONY: clean	
clean:
	rm -f *o
	rm -f *~
	rm -f *.dat
	rm -f *.eps
	rm -f $(OUT)
	
run:
	./main
