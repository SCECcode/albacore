#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(ALBACORE, 1.0, mei@usc.edu)
AC_CONFIG_SRCDIR([src])
AM_INIT_AUTOMAKE()

# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Set final CFLAGS and LDFLAGS
CFLAGS="$CHECK_CFLAGS $ETREE_INCL $PROJ4_INCL"
LDFLAGS="$CHECK_LDFLAGS $ETREE_LIB $PROJ4_LIB"
LDFLAGS="$LDFLAGS -lm"

AC_CONFIG_FILES([Makefile
		data/Makefile
		src/Makefile
		tests/Makefile])
AC_OUTPUT
