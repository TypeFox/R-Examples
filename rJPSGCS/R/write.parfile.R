write.parfile=function(snp.data,map,file="out.par"){

    MAFs <-  chopsticks::summary(snp.data)$MAF
    nSNPs=ncol(snp.data)
    
    #n loci, risk locus, sexlinked, program code
    cat(c(nSNPs,0,0,5), file = file,"\n")
    
    #mutation locus, mutation rates male, female, linkage disequil
    cat(c(0, 0.0, 0.0, 1), file = file,"\n",append=T)
    
    #Marker Order
    cat(1:nSNPs, file = file,"\n",append=T)
    
    #snps freq
    freqmat=cbind(1-MAFs,MAFs)
    
    x <- paste("3 2 ",colnames(snp.data),"\n",
    round(freqmat[,1],6)," ",round(freqmat[,2],6),"\n",
    sep="",collapse="")
    
    cat(x,file=file,append=T)

    #sexdifference, interference
    cat(c(0, 0), file = file,"\n",append=T)
    
    #recombination values
    dcm=rep(0,(nSNPs-1))
    for(i in 1:(nSNPs-1)){
        dcm[i]=abs(map[i]-map[i+1])
    }
    rf=0.5*(1-exp(-2*dcm/100))#Haldane mapping function
    
    cat(rf, file = file,"\n",append=T)
    
    #last line
    cat(c(1, 5, 0.2, 0.1 ), file = file, "\n",append=T)
}

