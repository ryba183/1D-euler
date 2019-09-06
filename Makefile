PROFILE_FLAGS= -O0 -pg -fprofile-arcs -ftest-coverage
OPT_FLAGS = -O3
PG_OPT_FLAGS = -Og
WARN_FLAGS = -Wall -Wextra -pedantic

1d_euler: 1d_euler.cpp
	g++ -o 1d_euler $(WARN_FLAGS) $(OPT_FLAGS) 1d_euler.cpp

1d_euler_pg: 1d_euler.cpp
	g++ -o 1d_euler_pg $(WARN_FLAGS) $(PROFILE_FLAGS) $(PG_OPT_FLAGS) 1d_euler.cpp
