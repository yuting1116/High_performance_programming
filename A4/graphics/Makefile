CFLAGS=-Wall -O0 -g -pthread
INCLUDES=-I/opt/X11/include
LDFLAGS=-L/opt/X11/lib -lX11 -lm
GRAPHICS = .

galsim: $(GRAPHICS)/galsim.o $(GRAPHICS)/graphics.o
	gcc -o galsim $(GRAPHICS)/galsim.o $(GRAPHICS)/graphics.o $(LDFLAGS)

galsim_pthread_lock: $(GRAPHICS)/galsim_pthread_lock.o $(GRAPHICS)/graphics.o
	gcc -o galsim_pthread_lock $(GRAPHICS)/galsim_pthread_lock.o $(GRAPHICS)/graphics.o $(LDFLAGS)

galsim.o: $(GRAPHICS)/galsim.c $(GRAPHICS)/graphics.h
	gcc $(CFLAGS) $(INCLUDES) -c $(GRAPHICS)/galsim.c -o galsim.o

galsim_pthread_lock.o: $(GRAPHICS)/galsim_pthread_lock.c $(GRAPHICS)/graphics.h
	gcc $(CFLAGS) $(INCLUDES) -c $(GRAPHICS)/galsim_pthread_lock.c -o galsim_pthread_lock.o

graphics.o: $(GRAPHICS)/graphics.c $(GRAPHICS)/graphics.h
	gcc $(CFLAGS) $(INCLUDES) -c $(GRAPHICS)/graphics.c

clean:
	rm -f $(GRAPHICS)/*.o
	rm -f ./galsim
	rm -f ./result.gal
	rm -f ../result.gal
