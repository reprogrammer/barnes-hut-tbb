CPP      := clang++
CPPFLAGS += $(CXXFLAGS) -std=c++11 #-stdlib=libc++
CHECKER  := $(CPP) --analyze -Xclang -analyzer-checker=alpha.SafeParallelismChecker -Xclang -analyzer-config -Xclang -asap-default-scheme=parametric
CHECKERFLAGS += $(CPPFLAGS)
INFERER  := $(CPP) --analyze -Xclang -analyzer-checker=alpha.SafeParallelismChecker -Xclang -analyzer-config -Xclang -asap-default-scheme=parametric-effect-inference -Xclang -verify
INFERERFLAGS += $(CPPFLAGS) -DASAP_INFER_EFFECTS
FULLINFERER := $(CPP) --analyze -Xclang -analyzer-checker=alpha.SafeParallelismChecker -Xclang -analyzer-config -Xclang -asap-default-scheme=inference -Xclang -verify
FULLINFERERFLAGS += $(INFERERFLAGS) -DASAP_INFER_ARGS
LD       := clang++
LIBS     += -ltbb
INCLUDES := 

# Remove these files when doing clean
OUTPUT +=

CPPFLAGS += 

PROG := driver/BarnesHut

ONEFILE := driver/BarnesHut.cc

SRCS += src/BarnesHut.cpp

# This order must support concatenating all the headers as one file with includes removed.
HEADERS = src/asap.h

OBJS := ${SRCS:.cpp=.o}


.PHONY: default
default: $(PROG)

.PHONY: clean
clean:
	$(RM) $(OBJS) $(PROG) $(OUTPUT) $(ONEFILE)

.PHONY: check
check: $(ONEFILE)
	$(CHECKER) $(CHECKERFLAGS) $(INCLUDES) $(ONEFILE)

.PHONY: infer
infer: $(ONEFILE)
	$(INFERER) $(INFERERFLAGS) $(INCLUDES) $(ONEFILE)

.PHONY: inferfull
inferfull: $(ONEFILE)
	$(FULLINFERER) $(FULLINFERERFLAGS) $(INCLUDES) $(ONEFILE)

$(PROG): $(OBJS)
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $(PROG)

driver/%.o: src/%.cpp src/*.h
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

$(ONEFILE): src/*.cpp src/*.h
	cat $(HEADERS) $(SRCS) | egrep -v "include \"[^.]" | egrep -v "\#line" > $(ONEFILE)

