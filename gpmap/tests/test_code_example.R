## Test code example 1 from 
## Gjuvsland AB, Wang Y, Plahte E, Omholt SW (2013) Monotonicity is a key feature of genotype-phenotype maps. Submitted to Frontiers in Genetics

library(gpmap)
data(GPmaps)
gp <- mouseweight

## Tabulate genotypic values
cbind(gp$genotype,gp$values)

## Plot the GP map
plot(gp)		

## Compute degree of monotonicity
gp <- degree_of_monotonicity(gp)
gp$degree.monotonicity.locus
print(gp)

## Quantify monotonicity by isotonic regression
gp <- decompose_monotone(gp)
print(gp)

## Plot decomposed GP map
plot(gp,decomposed=TRUE) 

