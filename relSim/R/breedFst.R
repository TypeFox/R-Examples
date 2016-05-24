breedFst = function(Freqs, theta = 0.01, N = 10000, ns = 10,
                    DNAtools = FALSE){
    if(N<1000){
        stop("N must be >= 1000")
    }

    if(N %% ns != 0){
        stop("ns must divide N into a whole number\tThat is the subpopulation sizes must be integers")
    }

    if(theta <=0 || theta >= 0.5){
        stop("0 < theta < 0.5")
    }

    Ns = N / ns

    nGen = ceiling(log(1 - theta) / log(1 - 1/(2*Ns)))
    cat(paste("Breeding for", nGen, "generations\n"))

    nLoci = length(Freqs$freqs)

    ## generate the parental population

    parents = randomProfiles(Freqs$freqs, N)
    pb = txtProgressBar(min = 0, max = nGen, style = 3)

    for(t in 1:nGen){
        parents = .breed(parents, ns, Ns, nLoci)
        setTxtProgressBar(pb, t)
    }
    setTxtProgressBar(pb, nGen)

    if(DNAtools){
        parents = data.frame(cbind(1:N, matrix(parents, nrow = N, byrow = T)))
        colNames = paste(rep(Freqs$loci, rep(2, nLoci)), c('a','b'), sep = '_')
        names(parents) = c('id', colNames)
    }

    pop = list(profiles = parents, nProfiles = N,  nSubpops = ns,
               nLoci = nLoci, theta = theta, Freqs = Freqs)
    class(pop) = "population"

    return(pop)
}
