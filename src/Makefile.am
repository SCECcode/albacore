# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_FCFLAGS = ${FCFLAGS}
AM_LDFLAGS = ${LDFLAGS}

TARGETS = libalbacore.a libalbacore.so

all: $(TARGETS)

install:
	mkdir -p ${prefix}
	mkdir -p ${prefix}/lib
	mkdir -p ${prefix}/include
	cp libalbacore.so ${prefix}/lib
	cp libalbacore.a ${prefix}/lib
	cp albacore.h ${prefix}/include

libalbacore.a: albacore_static.o 
	$(AR) rcs $@ $^

libalbacore.so: albacore.o
	$(FC) -shared $(AM_FCFLAGS) -o libalbacore.so $^ $(AM_LDFLAGS)

albacore.o: albacore.c
	$(CC) -fPIC -DDYNAMIC_LIBRARY -o $@ -c $^ $(AM_CFLAGS)
	
albacore_static.o: albacore.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)
	
utm_geo.o: utm_geo.f90
	$(FC) -fPIC -o $@ -c $^ $(AM_FCFLAGS)
