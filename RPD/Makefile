CCC = g++
FLAGS = -g -O3 -Wall -Wextra

main: main.cpp path_equilibration.o path_equilibration_rho.o bbstep_PGM.o
	g++ -g -O3 -Wall -DNDEBUG -std=c++11 main.cpp path_equilibration.o path_equilibration_rho.o bbstep_PGM.o -o main -lquadmath -fopenmp
path_equilibration.o: src/path_equilibration.cpp
	g++ -O3 -Wall -DNDEBUG -lquadmath -fopenmp -std=c++11 -MMD -MP -c src/path_equilibration.cpp
path_equilibration_rho.o: src/path_equilibration_rho.cpp
	g++ -O3 -Wall -DNDEBUG -lquadmath -fopenmp -std=c++11 -MMD -MP -c src/path_equilibration_rho.cpp
bbstep_PGM.o: src/bbstep_PGM.cpp
	g++ -O3 -Wall -DNDEBUG -std=c++11 -MMD -MP -c src/bbstep_PGM.cpp