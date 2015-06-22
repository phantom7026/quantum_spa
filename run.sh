
# USEGA:	mpirun -np A B file D E
#
# A : the number of threads
# B : execution name
# file : input file
# D : the number of tests
# E : GF extension, m of GF(2^m)

EXEC=spa
TEST=10
EXTENSION=1

SOURCE=("codes/nb_generator_50_deg1.txt")

for value in "${SOURCE[@]}"; do
	mpirun -np 4 $EXEC $value $TEST $EXTENSION
done
