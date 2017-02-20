
target = iqat
objects = iqat.o
prefix = /usr/local
bindir = $(prefix)/bin

CC = g++
CFLAGS = -O3 -std=c++11 -D__STDC_CONSTANT_MACROS -fopenmp -Wall -Wno-sign-compare -Wno-write-strings -Wno-deprecated-declarations
LD = g++
LDFLAGS = -fopenmp -lavutil -lavcodec -lavformat -lswscale
STRIP = strip

all: $(objects)
	$(LD) $(LDFLAGS) $(objects) -o $(target)
	$(STRIP) $(target)

$(objects): %.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(target) $(objects)

install:
	cp $(target) $(bindir)

uninstall:
	rm -f $(bindir)/$(target)

