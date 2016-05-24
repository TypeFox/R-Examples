checkFreqs = function(Freqs){
    sums = sapply(Freqs$freqs, sum)
    ranges = sapply(Freqs$freqs, function(x){any(x < 0 | x > 1)})

    if(any(sums!=1)){
        Loci = Freqs$loci[sums!=1]
        sums = sums[sums!=1]
        cat(paste("The locus", paste(Loci,sums), "does not sum to 1\n"))
    }

    if(any(ranges)){
        Loci = Freqs$loci[ranges]
        cat(paste("The loci", paste(Loci),
                  "have values less than 0 or greater than 1\n"))
    }
}
