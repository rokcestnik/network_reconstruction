# Reconstructing networks of pulse-coupled oscillators

This is an implementation of our algorithm for reconstructing networks of pulse-coupled oscillators from passive observation of pulse trains of all nodes. It is assumed that units are described by their phase response curves and that their phases are instantaneously reset by incoming pulses. Using an iterative procedure, we recover the properties of all nodes, namely their phase response curves and natural frequencies, as well as strengths of all directed connections. For details refer to our [paper](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.96.012209). 

## Setup

Put the [Eigen](http://eigen.tuxfamily.org) header files into the 'Eigen' folder.

## Execution

Just present the data (see how bellow) and run the bash script [run.bash](run.bash) in a Unix terminal:
```
bash run.bash
```

## Data

Put the data in the folder 'data'. Each oscillator in its own '.txt' file, containing only the spike times in the format "t(i)\n" (strictly increasing). Name the files with numbers starting from 0. Do not put any other files in the folder.

### Data example:

|'0.txt' |'1.txt' |'2.txt' | . . .  |'20.txt'|
|--------|--------|--------|--------|--------|
|7.75    |3.05    |4.55    | . . .  |6.      |
|8.      |7.2     |6.7     | . . .  |6.65    |
|8.3     |10.25   |8.8     | . . .  |7.8     |
|8.85    |11.55   |10.05   | . . .  |8.85    |
|9.15    |12.05   |11.7    | . . .  |10.15   |
|9.8     |14.45   |13.9    | . . .  |13.55   |
|...     |...     |...     | . . .  |...     |


## Reconstruction

The reconstruction will be saved in the folder 'reconstruction'. 

* connectivity in the file 'connectivity.txt' with the element in the n-th row and m-th column representing the strength of the connections from unit m to unit n, 
* natural frequencies in 'frequencies.txt' with each line representing the natural frequency of a unit, 
* phase response curves in 'PRCs.txt' with each line representing the Fourier expansion of the PRC for one unit, 
* errors of reconstruction in 'err.txt' with the two columns representing the errors Δ<sub>ψ<sub>T</sub></sub> and Δ<sub>ψ</sub> respectively, defined in this [paper](https://www.nature.com/articles/s41598-018-32069-y), 
* reconstruction score in 'rec_score.txt' is just the logarithm of the quocient of the two erors and can be interpreted as "the bigger the better". 
