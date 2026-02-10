CFLAGS=-Wall
INCLUDES=-I/opt/X11/include
LDFLAGS=-L/opt/X11/lib -lX11 -lm
GRAPHICS = ./graphics
galsim: $(GRAPHICS)/galsim.o $(GRAPHICS)/graphics.o
	gcc -o galsim $(GRAPHICS)/galsim.o $(GRAPHICS)/graphics.o $(LDFLAGS)

galsim.o: $(GRAPHICS)/galsim.c $(GRAPHICS)/graphics.h
	gcc $(CFLAGS) $(INCLUDES) -c $(GRAPHICS)/galsim.c

graphics.o: $(GRAPHICS)/graphics.c $(GRAPHICS)/graphics.h
	gcc $(CFLAGS) $(INCLUDES) -c $(GRAPHICS)/graphics.c

clean:
	rm -f $(GRAPHICS)/*.o
	rm -f ./galsim
