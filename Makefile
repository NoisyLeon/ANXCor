INST_DIR = $(HOME)/bin
EXECUTABLE  = ANXCor

cflags = -std=c++0x -O2 -fopenmp #-Wall -I${HOME}/usr/include

LDLIBS = -L${HOME}/usr/lib -fopenmp -lfftw3 

CFLAGS = $(DBG) $(cflags)

CC = g++

DBG = -g 

all : $(EXECUTABLE)

FOBJS = CCList.o FTNorm.o CCRec.o vincenty.o CCDatabase.o SeedRec.o SacRec.o SysTools.o $(EXECUTABLE).o

$(EXECUTABLE) : $(FOBJS)
	$(CC) -o $@ $^ $(LDLIBS) $(CFLAGS)

%.o : %.cpp
	$(CC) $(cflags) -c $<

install : $(EXECUTABLE)
	install -s $(EXECUTABLE) $(INST_DIR)

install_ariadne : $(EXECUTABLE)
	cp $(EXECUTABLE) $(INST_DIR)/$(EXECUTABLE)_ariadne

clean :
	rm -f $(EXECUTABLE) core $(FOBJS)
