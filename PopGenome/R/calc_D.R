calc_D <- function(bial,populations,outgroup){

# Patterson's D
# f

# pop 1: populations[[1]] :
# pop 2: populations[[2]] : 
# pop 3: populations[[3]] : archaic population

if(length(populations)!=3){
stop("This statistic requires 3 populations; the third population is the archaic population")
}

if(outgroup[1]==FALSE || length(outgroup[1])==0){
stop("This statistic needs an outgroup ! (set.outgroup)")
}

# calculate frequencies of the derived alleles.
freqs  <- jointfreqdist(bial,populations,outgroup)
freqs  <- freqs$jfd

# calc D
ABBA   <- (1-freqs[1,])*freqs[2,]*freqs[3,]
BABA   <- freqs[1,]*(1-freqs[2,])*freqs[3,]

sum_ABBA <- sum(ABBA,na.rm=TRUE)
sum_BABA <- sum(BABA,na.rm=TRUE) 

D <- (sum_ABBA - sum_BABA)/(sum_ABBA + sum_BABA)

# calc f
freqs23    <- freqs[2:3,]
maxfreqs23 <- apply(freqs23,2,max)

maxABBA <- (1-freqs[1,])*maxfreqs23*maxfreqs23
maxBABA <- freqs[1,]*(maxfreqs23)*maxfreqs23

sum_maxABBA <- sum(maxABBA, na.rm=TRUE)
sum_maxBABA <- sum(maxBABA, na.rm=TRUE)

f <- (sum_ABBA - sum_BABA)/(sum_maxABBA + sum_maxBABA)


return(list(D=D,f=f))

}
