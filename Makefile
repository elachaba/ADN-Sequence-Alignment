# Makefile qui genere l'executable distanceEdition et fait des tests de verification
#
#
CC=gcc
LATEXC=pdflatex
DOCC=doxygen
CFLAGS=-g -Wall 

REFDIR=.
SRCDIR=$(REFDIR)/src
BINDIR=$(REFDIR)/bin
DOCDIR=$(REFDIR)/doc
TESTDIR=$(REFDIR)/tests
REPORTDIR=$(REFDIR)/report
CURRTEST = Needleman-Wunsch-iter.o

LATEXSOURCE=$(wildcard $(REPORTDIR)/*.tex)
CSOURCE=$(wildcard $(SRCDIR)/*.c)
PDF=$(LATEXSOURCE:.tex=.pdf)

all: binary report doc binary_perf

binary: $(BINDIR)/distanceEdition

binary_perf: $(BINDIR)/distanceEdition-perf

$(BINDIR)/distanceEdition-perf: $(SRCDIR)/distanceEdition.c $(BINDIR)/$(CURRTEST)
	$(CC) $(OPT) -D__PERF_MESURE__ -I$(SRCDIR) -o $(BINDIR)/distanceEdition-perf $(BINDIR)/$(CURRTEST) $(SRCDIR)/distanceEdition.c

report: $(PDF) 

doc: $(DOCDIR)/index.html


$(BINDIR)/distanceEdition: $(SRCDIR)/distanceEdition.c $(BINDIR)/$(CURRTEST)
	$(CC) $(OPT) -I$(SRCDIR) $(CFLAGS) -o $(BINDIR)/distanceEdition $(BINDIR)/$(CURRTEST) $(SRCDIR)/distanceEdition.c

$(BINDIR)/Needleman-Wunsch-recmemo.o: $(SRCDIR)/Needleman-Wunsch-recmemo.h $(SRCDIR)/Needleman-Wunsch-recmemo.c $(SRCDIR)/characters_to_base.h
	$(CC) $(OPT) -I$(SRCDIR) -c  -o $(BINDIR)/Needleman-Wunsch-recmemo.o $(SRCDIR)/Needleman-Wunsch-recmemo.c

$(BINDIR)/Needleman-Wunsch-iter.o: $(SRCDIR)/Needleman-Wunsch-iter.h $(SRCDIR)/Needleman-Wunsch-iter.c $(SRCDIR)/characters_to_base.h
	$(CC) $(OPT) -I$(SRCDIR) -c  -o $(BINDIR)/Needleman-Wunsch-iter.o $(SRCDIR)/Needleman-Wunsch-iter.c

$(BINDIR)/Needleman-Wunsch-CA.o: $(SRCDIR)/Needleman-Wunsch-CA.h $(SRCDIR)/Needleman-Wunsch-CA.c $(SRCDIR)/characters_to_base.h
	$(CC) $(OPT) -I$(SRCDIR) -c  -o $(BINDIR)/Needleman-Wunsch-CA.o $(SRCDIR)/Needleman-Wunsch-CA.c

$(BINDIR)/Needleman-Wunsch-CO.o: $(SRCDIR)/Needleman-Wunsch-CO.h $(SRCDIR)/Needleman-Wunsch-CO.c $(SRCDIR)/characters_to_base.h
	$(CC) $(OPT) -I$(SRCDIR) -c -g -o $(BINDIR)/Needleman-Wunsch-CO.o $(SRCDIR)/Needleman-Wunsch-CO.c

$(BINDIR)/extract-fasta-sequences-size: $(SRCDIR)/extract-fasta-sequences-size.c
	$(CC) $(OPT) -I$(SRCDIR) -o $(BINDIR)/extract-fasta-sequences-size $(SRCDIR)/extract-fasta-sequences-size.c

clean:
	rm -rf $(DOCDIR) $(BINDIR)/* $(REPORTDIR)/*.aux $(REPORTDIR)/*.log  $(REPORTDIR)/rapport.pdf 

#$(BINDIR)/distanceEdition: $(CSOURCE)
#	$(CC) $(CFLAGS)  $^ -o $@ 

$(BINDIR)/distanceEditiondebug: $(CSOURCE)
	$(CC) $(CFLAGS)  $^ -o $@ -DDEBUG

%.pdf: $(LATEXSOURCE)
	$(LATEXC) -output-directory $(REPORTDIR) $^ 

$(DOCDIR)/index.html: $(SRCDIR)/Doxyfile $(CSOURCE)
	$(DOCC) $(SRCDIR)/Doxyfile


test: $(BINDIR)/distanceEdition $(TESTDIR)/Makefile-test
	cd $(TESTDIR) ; make -f Makefile-test all 
	
test-valgrind: $(BINDIR)/distanceEdition $(TESTDIR)/Makefile-test
	make -f $(TESTDIR)/Makefile-test all-valgrind
	
.PHONY: all doc bin report 

