# Other include directory (for CFITSIO, libsla, which is in PRESTO)
OTHERINCLUDE = -I/usr/include/cfitsio  #-I/home.local/phil/svn/pdev/include
# Other link directory (for CFITSIO)
OTHERLINK = -L/usr/lib64 -lcfitsio #-L/home.local/phil/svn/pdev/libs

# Source directory
SRCDIR = $(shell pwd)

#BINDIR = /home/deneva/bin64
BINDIR = $(shell pwd)

# Which C compiler
CC = gcc
CFLAGS = $(OTHERINCLUDE) -DSRCDIR=\"$(SRCDIR)\"\
        -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64\
        -g -Wall -W -O1 \

# When modifying the CLIG files, the is the location of the clig binary
#CLIG = /usr/bin/clig
CLIG =/home/deneva/local/bin64/clig
# Rules for CLIG generated files
%_cmd.c : %_cmd.cli
	$(CLIG) -o $*_cmd -d $<
	cp $*_cmd.c ..

OBJS = fixbeampos.o read_psrfits.o write_psrfits.o alfaoff.o cal2mjd.o prec.o azza.o deg_trig.o radec.o

fixbeampos: $(OBJS)
	$(CC) $(CFLAGS) -o $(BINDIR)/$@ $(OBJS) -lm $(OTHERLINK)

clean:
	rm -f *.o *~ *#

