LS2Wsim.cddews <-function (spectrum,innov=rnorm,...) 
{
    datadim <- spectrum$datadim
    spec2d <- cdtoimwd(spectrum)  

    nlev <- spec2d$nlev
    len <- (2^nlev)*(2^nlev)
    newspec2d <- spec2d
 
    for (i in (nlev - 1):0) {

    tmpV<-lt.to.name(i,"DC")
    tmpH<-lt.to.name(i,"CD")
    tmpD<-lt.to.name(i,"DD")  

     v1 <- spec2d[[tmpV]]
     v2 <- spec2d[[tmpH]]
     v3 <- spec2d[[tmpD]]
 
     v1[v1<0] <- 0
     v2[v2<0] <- 0
     v3[v3<0] <- 0

     simvDC <- sqrt(v1) * 4^(nlev - i) * innov(len,...)
     simvCD <- sqrt(v2) * 4^(nlev - i) * innov(len,...)
     simvDD <- sqrt(v3) * 4^(nlev - i) * innov(len,...)
            
     newspec2d[[tmpV]] <- simvDC
     newspec2d[[tmpH]] <- simvCD
     newspec2d[[tmpD]] <- simvDD   
}
    
sim<-AvBasis.wst2D(convertimwd(newspec2d))

return(sim)
    
}

