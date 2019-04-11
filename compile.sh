g++ `root-config --cflags` -g -std=c++11 -c MatCalc.C
g++ `root-config --cflags` -g -std=c++11 main.C -o matcalc `root-config --libs` MatCalc.o

