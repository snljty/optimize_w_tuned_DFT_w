# Makefile for optimize_DFT_w

SHELL = cmd
CC = gcc
CLINKER = $(CC)
CLINKERFLAGS = -static
CCFLAGS = -O3
ARCH = ar
ARCHFLAGS = -rsc

LIBNAME := brent_fmin
TARGETNAME = optimize_DFT_w

.PHONY: all
all: lib $(TARGETNAME)

.PHONY: lib
lib: lib$(LIBNAME).a

lib$(LIBNAME).a: $(LIBNAME).obj
	@echo Generating archive $@ from $^ ...
	$(ARCH) $(ARCHFLAGS) $@ $^

$(LIBNAME).obj: $(LIBNAME).c
	@echo Compiling $@ ...
	$(CC) -o $@ -c $< $(CCFLAGS)

.PHONY: $(TARGETNAME)
$(TARGETNAME): $(TARGETNAME).exe

$(TARGETNAME).exe: $(TARGETNAME).obj lib$(LIBNAME).a
	@echo Linking $@ against $^ ...
	$(CLINKER) -o $@ $< -L . -l $(LIBNAME) $(CLINKERFLAGS)

$(TARGETNAME).obj: $(TARGETNAME).c
	@echo Compiling $@ ...
	$(CC) -o $@ -c $< $(CCFLAGS)

.PHONY: clean
clean: clean_tmp
	-del /q $(TARGETNAME).exe 2> NUL

.PHONY: clean_tmp
clean_tmp:
	-del /q $(LIBNAME).obj 2> NUL
	-del /q $(TARGETNAME).obj 2> NUL
	-del /q lib$(LIBNAME).a 2> NUL

.PHONY: clean_$(TARGETNAME)
clean_$(TARGETNAME): clean_exe
	-del /q $(TARGETNAME).obj 2> NUL

