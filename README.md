# PCR

A .R file allowing to calculate errors in PCR. 
This is the code I am using to run this [online tool](https://biosoft-ipcms.fr/files/code.php).

## Background 
- ThermoFisher has a [PCR Calculator](https://www.thermofisher.com/uk/en/home/brands/thermo-scientific/molecular-biology/molecular-biology-learning-center/molecular-biology-resource-library/thermo-scientific-web-tools/pcr-fidelity-calculator.html), which is basic and incaccurate. It  just caculates ~ $1-l \times c \times u$ ( $l$, fragment length; $u$, error rate and $c$ the number of cycles), which is nothing else than the probality for a molecule of length $l \times c$ to have no errors. Thus, it assumes 100% efficiency and ignores that some molecules are not replicated and so contribute to the overall number of molecules without errors. 

- Here, I use a more rigourous approach (a good article can be found [here](https://doi.org/10.3929/ethz-a-006088024)).

- Note that, New England Biolabs also has an online  [PCR fidelity estimator](https://pcrfidelityestimator.neb.com) but they do not explain how the calculations are performed.


As for the error rates, you may use this :



| Enzyme| 	Error Rate| 	Fidelity rel. Taq| 
| ------------- | ------------- | ------------- |
| Taq	| 1.5E-4	| 1| 
| Q5  | 5.3E-7	| 280| 
| Phusion	| 3.9E-6	| 39| 
| Deep Vent	| 4.0E-6	| 44| 
| Pfu	| 5.1E-6| 	30| 
| PrimeSTAR| 	8.4E-6	| 18| 
| KOD| 	1.2E-5| 	12| 
| Kapa HIFI| 	1.6E-5	| 9.4| 
| DV Exo-	| 5.0E-4| 	0.3| 

Adapted from Potapov V, Ong JL [PLOS ONE 2017](https://doi.org/10.1371/journal.pone.0169774)

## How to use

Just change the numbers at the begining of the code.
```
# Rate
u<-1.5E-4
# Fragment length
l<-2000
# Number of Cycles
numcycles<-25
```

