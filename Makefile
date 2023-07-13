CXX=/opt/intel/oneapi/compiler/latest/linux/bin/icpx
CXXFLAGS = -xCORE-AVX512 -std=c++17 -Werror
DEBUGFLAGS = -O0 -g -DVTUNE -L/opt/intel/oneapi/vtune/latest/sdk/lib64/ -ldl -pthread
OPTFLAGS = -O3

main: nuovo_main.cpp temporal_graph.hpp cuckoo.hpp span.hpp Makefile
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) nuovo_main.cpp -o $@

debug: nuovo_main.cpp temporal_graph.hpp cuckoo.hpp span.hpp Makefile
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) nuovo_main.cpp -o main -littnotify

.PHONY: clean
clean:
	rm -rf main *.o
