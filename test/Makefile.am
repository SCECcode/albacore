# Autoconf/automake file

bin_PROGRAMS = test_albacore

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_LDFLAGS = ${LDFLAGS}

objects = test_api.o
TARGETS = $(bin_PROGRAMS)

all: $(bin_PROGRAMS)

install:
	mkdir -p ${prefix}/tests
	cp test_albacore ${prefix}/tests

test_albacore: $(objects)
	$(CC) -o $@ $^ $(AM_LDFLAGS)

$(objects): %.o: %.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS) 

run_unit : test_albacore
	./test_albacore

clean :
        rm -rf *~ *.o test_alabacore
