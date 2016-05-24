# Function called within allelefreq.r, to calculate the relative
# proportion of each allele in all populations for a single locus
pop.allele.proportions <- function(anzahlk.pop.1pop){
allele.numbers=as.numeric(as.vector(anzahlk.pop.1pop$number))
hundred.percent=sum(allele.numbers)
percent=allele.numbers/hundred.percent
percent
}
