#parameters
n="$(ls ./data -1 | wc -l)" #count the data files
Nfourier=6 #fourier terms
FEE=1 #first estimation epsilon (1=yes, 0=no)
N=10 #coarse graining for FEE
its=10 #how many iterations

#write parameters
g++ ./programs/write_params.c -o ./programs/write_params.out
./programs/write_params.out $n $Nfourier $FEE $N

#compile programs
g++ ./programs/getting_eps.c -o ./programs/getting_eps.out
g++ ./programs/getting_prc.c -o ./programs/getting_prc.out
g++ ./programs/improving_eps.c -o ./programs/improving_eps.out
g++ ./programs/improving_prc.c -o ./programs/improving_prc.out

let "n--" #n-- because we start counting from 0
for observed in `seq 0 $n`;
do

	echo "node $observed"
		
	echo -e "\tfirst approximation"
	./programs/getting_eps.out $observed
	./programs/getting_prc.out $observed
	#subsequent iterations
	for it in `seq 2 $its`;
	do
		echo -e "\titeration $it"
		./programs/improving_eps.out $observed
		./programs/improving_prc.out $observed
	done

	#save recovery
	scp EPS.txt reconstruction/EPS$observed.txt
	scp PRC.txt reconstruction/PRC$observed.txt
	
done

#compile connectivity matrix
g++ ./programs/compile_connectivity_matrix.c -o ./programs/compile_connectivity_matrix.out
./programs/compile_connectivity_matrix.out

#measure errors
g++ ./programs/error.c -o ./programs/error.out
./programs/error.out

#clean the main directory
rm EPS.txt
rm PRC.txt
rm parameters.h
