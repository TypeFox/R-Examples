gen2clear.2fis <- function(k, gen){
    ncand <- choose(k+length(gen),2)
    hilf <- combn(k+length(gen),2)
    sel <- 1:ncand
    ww <- words.all(k, gen, max.length=4)
    if (length(ww$WLP)==0) return(list(nclear.2fis=ncand, clear.2fis=hilf))
    else{
                woerter <- ww[[2]]
                nclear.2fis <- 0
                spalte <- numeric(0)
                for (i in 1:ncol(hilf))
                    if (all(sapply(woerter, function(obj) sum(hilf[,i] %in% obj)<2))){
                      nclear.2fis <- nclear.2fis+1
                      spalte <- c(spalte,sel[i])
                    }
                return(list(nclear.2fis=nclear.2fis, clear.2fis=hilf[,spalte]))
    }
}

gen2CIG <- function(nruns, gen){
    if (!is.numeric(nruns)) stop("nruns must be a number")
    if (!length(nruns)==1) stop("nruns must be a single number")
    k <- round(log2(nruns))
    if (!nruns==2^k) stop("nruns must be a power of 2")
    if (!is.numeric(gen)) stop("gen must be a numeric vector")
    if (!all(gen==round(gen))) stop("gen must contain integer numbers only")
    hilf <- gen2clear.2fis(k, gen)
    if (hilf$nclear.2fis==0) return(graph(numeric(0), n=k+length(gen)))
    else{
    return(graph(hilf$clear.2fis, n=k+length(gen),directed=FALSE))
    }
}