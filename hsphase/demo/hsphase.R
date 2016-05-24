# example of using hsphase functions
# example files consist of a genotypes file with 90 individuals across 2 families and 2 chromosomes

# if R's working directory has the 3 files - genotypes.txt, map.txt and pedigree.txt; the example below should run well 
# files can be downloaded from http://www-personal.une.edu.au/~cgondro2/hsphase/hsphaseExample.zip 

# to install the dependencies
# install.packages("gtools") 
# install.packages("snowfall")
# install.packages("Rcpp") 
# install.packages("RcppArmadillo")

# install.packages("hsphase") in case the package is not installed

#######################################################################
################ Running hsphase ######################################
#######################################################################

library(hsphase) # load library

# reads in a file of genotypes and a pedigree file, then splits the data into a list (of matrices - ID x SNP) of half-sibs, one for each family
# check the example files to see the format
# genotypes: ID x SNP, with rownames = ID and colnames = SNP
# genotypes: 0 - AA, 1 - AB, 2 - BB, 9 missing
# pedigree: 1st column ID, 2nd column Sire ID. Should not have a header

data(genotypes)

data(pedigree)

# split into list of half-sib groups 
halfsib <- hss(pedigree, genotypes)

# splits the family data into chromosomes
# also a list of matrices (ID x SNP): one for each chromosome X number of families
# per family, per chromosome matrices are ordered according to base pair position in the map file
# names in the list are the sireID_chromosome - can use this to find/parse subsets of data
data(map)
halfsib <- cs(halfsib,map)

# function para runs a parallel version of the main functions - it's probably what you want in practice
blocks <- para(halfsib, cpus=2, option="bmh", type = "SOCK") # bmh - builds blocks of relationship between paternal and maternal strand from sire (i.e. which chunks each offspring inherited from the sire)
sires <- para(halfsib, cpus=2, option="ssp", type = "SOCK") # ssp - phases and imputes the sire 
phased <- para(halfsib, cpus=2, option="aio", type = "SOCK") # aio - phases half-sib families (p - paternal strand from sire and m - maternal strand from dam)

# any of these function can also be run directly on a single matrix with e.g.
singleChromSirePhased <- aio(halfsib[[1]])


#######################################################################
################ ploting with hsphase #################################
#######################################################################

# visualize block structure for sire 1, chromosome 1
imageplot(blocks[[1]])

# plot number of recombinations between SNP
# needs distances between SNP - can get from map
distance <- map$Position[which(map$Chr==1)]
# plot
rplot(halfsib[[1]],distance) # plot recombinations for family one, chromosome 1


#######################################################################
################ Diagnostic tools #####################################
#######################################################################

# generates a matrix of 0/1 for recombinations
# dimension is 1 less than the number of SNP
# if a particular region cannot be resolved - all SNPs in the region are set to 1
# can help identify problems if e.g. too many recombinations are identified
recombinations <- pm(blocks[[1]])

# matches the phased haplotypes of the offspring with the phased haplotypes of its sire
# inputs are a matrix for each of the haplotypes in 0 (A), 1 (B) and 9 (missing) format
diagnostic <- hbp(phased[[1]],sires[[1]])
imageplot(diagnostic)  
# pedigree errors can be detected if there's an individual with too many recombinations
# problems with the phased data are evident if there are excessive numbers of recombinations
# this function can be used as a diagnostic tool for other phasing algorithms


#######################################################################
################ Pedigree inference and fix pedigree error ############
#######################################################################

oh <- ohg(genotypes)  # create a matrix of opposing homozygotes
hh(oh) # heatmap plot of halfsibs without pedigree colour coded sidebars

inferredPedigree <- rpoh(genotypes[, 1:ncol(halfsib[[1]])], oh, forwardVectorSize = 30, excludeFP = TRUE, nsap = 3, maxRec = 15, method = "recombinations")  # infer half-sib groups pedigree 

inferredPedigree <- pedigreeNaming(inferredPedigree, pedigree) # assigns inferred pedigree to original pedigree
hh(oh, inferredPedigree, pedigree)  # heatmap with colour coded bars for the inferred and original pedigrees











