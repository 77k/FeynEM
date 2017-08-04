CPP     	= 	g++ -O3  -std=c++1z -fexceptions -D__GLIBCXX_TYPE_INT_N_0=__int128 -D__GLIBCXX_BITSIZE_INT_N_0=128
STRIPPER 	= 	./removecomments.pl
EMCC    	=       ~/Extern/emscripten3/emscripten/emcc -Wc++11-extensions -std=c++11
OPT             = 	 -s PRECISE_I64_MATH=1
SET     	= 	-O2 -s ALLOW_MEMORY_GROWTH=1 -s ASM_JS=0
LIBPATH 	= 
SYSLIBS 	= 	-lstdc++ -lboost_regex -lboost_thread -lboost_system -lboost_date_time -lrt  -lpthread -lm
PRGS		=	feynemd


all: le_build feynemd feynem_test

feynem_test: le_build/feynem_test.o
	$(CPP) -fexceptions -D_BOOL  $^ -I$(INCPATH) -L$(LIBPATH) $(SYSLIBS) -o $@
le_build/feynem_test.o: feynem_test.cpp cw_complex.h spaces.h algebra.h math_helpers.h qed.h
	$(CPP) -fexceptions -D_BOOL -c $< -I$(INCPATH) -L$(LIBPATH)  -DHAVE_IOMANIP -DHAVE_IOSTREAM -DHAVE_LIMITS_H -o $@

feynemd: le_build/feynemd.o
	$(CPP) -fexceptions -D_BOOL  $^ -I$(INCPATH) -L$(LIBPATH) $(SYSLIBS) -o $@

clean:
	rm -rf le_build feynemd feynem_test;\

le_build/feynemd.o: feynemd.cpp feynemd.h simple_httpd.h cw_complex.h spaces.h algebra.h math_helpers.h qed.h
	$(CPP) -fexceptions -D_BOOL -c $< -I$(INCPATH) -L$(LIBPATH)  -DHAVE_IOMANIP -DHAVE_IOSTREAM -DHAVE_LIMITS_H -o $@

le_build:
	mkdir ./le_build

feynemcqednt-opt.js: feynemcqednt.js
	closure-compiler --language_in=ECMASCRIPT5  --compilation_level ADVANCED_OPTIMIZATIONS --js $< --warning_level QUIET > $@ && cp ./$@ ./data/
feynemcqednt.js: main.js webglfoo.js qed-stripped.js LA_helpers.js
	cat $^ > $@ && cp ./$@ ./data/

qed.js: qed.sym qed.cpp qed.h
	$(EMCC) $(OPT) $<  -o $@ -s EXPORTED_FUNCTIONS=$(shell cat $^) $(SET)
qed-stripped.js: qed.js
	$(STRIPPER) $^ > $@
qed.sym: qed.cpp
symbols=\"[`awk '$$1 ~/EMCEXPORT/{sub(/\(.*/,"");printf "\x27_"$$3"\x27,"}' $<`]\";echo $$symbols > $@; echo "Functions to export: " $$symbols;


