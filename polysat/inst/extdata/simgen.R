#### You can use this R code to generate a simulated dataset similar
#### to the simgen object provided in the package and described
#### in the tutorial.

# set up some allele frequencies
pop1loc1 <- c(rep(100, 9), rep(102, 13), rep(106, 24), rep(108, 5),
              rep(110, 15), rep(112, 18), rep(114, 10), rep(118, 6))
pop2loc1 <- c(rep(102, 20), rep(104, 10), rep(106, 7), rep(110, 15),
              rep(114, 9), rep(116, 11), rep(118, 28))
pop3loc1 <- c(rep(100, 17), rep(102, 4), rep(104, 30), rep(108, 12),
              rep(110, 9), rep(112, 19), rep(114, 11), rep(116, 28),
              rep(118, 7))
pop1loc2 <- c(rep(146, 14), rep(152, 2), rep(155, 18), rep(158, 25),
              rep(161, 19), rep(164, 7))
pop2loc2 <- c(rep(143, 20), rep(146, 3), rep(149, 15), rep(158, 11),
              rep(161, 9))
pop3loc2 <- c(rep(143, 8), rep(146, 18), rep(149, 12), rep(152, 23),
              rep(155, 14), rep(164, 17))
pop1loc3 <- c(rep(210, 8), rep(214, 12), rep(216, 7), rep(218, 30),
              rep(222, 15), rep(228, 9), rep(230, 21))
pop2loc3 <- c(rep(212, 11), rep(218, 18), rep(222, 10), rep(224, 7),
              rep(228, 2), rep(230, 5))
pop3loc3 <- c(rep(210, 25), rep(212, 6), rep(216, 17), rep(218, 12),
              rep(220, 19), rep(222, 8), rep(224, 10), rep(226, 22))

# set up object that will hold genotypes
simgen <- new("genambig", samples=c(paste("A", 1:100, sep=""),
                          paste("B", 1:100, sep=""), paste("C", 1:100, sep="")),
              loci= c("loc1", "loc2", "loc3"))
Description(simgen) <- "Simulated dataset"
Usatnts(simgen) <- c(2, 3, 2)
PopInfo(simgen) <- c(rep(1, 100), rep(2, 100), rep(3, 100))
simgen <- reformatPloidies(simgen, output="sample")
Ploidies(simgen) <- sample(c(2,4), 300, replace=TRUE)

# simulate unambiguous genotypes
unambiggen <- array(list(-9), dim=c(300, 3), dimnames=list(Samples(simgen),
                                             Loci(simgen)))
for(s in Samples(simgen, populations=1)){
    unambiggen[[s, "loc1"]] <- sample(pop1loc1, Ploidies(simgen)[s],
                                       replace=TRUE)
}
for(s in Samples(simgen, populations=2)){
    unambiggen[[s, "loc1"]] <- sample(pop2loc1, Ploidies(simgen)[s],
                                       replace=TRUE)
}
for(s in Samples(simgen, populations=3)){
    unambiggen[[s, "loc1"]] <- sample(pop3loc1, Ploidies(simgen)[s],
                                       replace=TRUE)
}
for(s in Samples(simgen, populations=1)){
    unambiggen[[s, "loc2"]] <- sample(pop1loc2, Ploidies(simgen)[s],
                                       replace=TRUE)
}
for(s in Samples(simgen, populations=2)){
    unambiggen[[s, "loc2"]] <- sample(pop2loc2, Ploidies(simgen)[s],
                                       replace=TRUE)
}
for(s in Samples(simgen, populations=3)){
    unambiggen[[s, "loc2"]] <- sample(pop3loc2, Ploidies(simgen)[s],
                                       replace=TRUE)
}
for(s in Samples(simgen, populations=1)){
    unambiggen[[s, "loc3"]] <- sample(pop1loc3, Ploidies(simgen)[s],
                                       replace=TRUE)
}
for(s in Samples(simgen, populations=2)){
    unambiggen[[s, "loc3"]] <- sample(pop2loc3, Ploidies(simgen)[s],
                                       replace=TRUE)
}
for(s in Samples(simgen, populations=3)){
    unambiggen[[s, "loc3"]] <- sample(pop3loc3, Ploidies(simgen)[s],
                                       replace=TRUE)
}

# write ambiguous genotypes to simgen
for(L in Loci(simgen)){
    for(s in Samples(simgen)){
        Genotype(simgen, s, L) <- unique(unambiggen[[s, L]])
    }
}

# add some missing genotypes
for(i in 1:5){
    Genotype(simgen, sample(Samples(simgen), 1), sample(Loci(simgen), 1)) <-
        Missing(simgen)
}

rm(pop1loc1, pop2loc1, pop3loc1, pop1loc2, pop2loc2, pop3loc2, pop1loc3,
   pop2loc3, pop3loc3, unambiggen)

# try: meandistance.matrix, cmdscale, plot (with colors for pops)
# deSilvaFreq
