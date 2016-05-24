findBest <-
function(x,y,bd,di,method){
    n <- length(y); bestC <- 0; cacc <- 0;i <- 0;addCut <- NULL
    iBest <- NULL; newDi <- di
    nb <- numeric()
    for(i in 1:length(di)) {
        iw <- which(di[i]>bd)
        nb[iw] <- bd[iw]
        nl <- length(iw)
        nb[nl+1] <- di[i]
        nb[(nl+2):(length(bd)+1)] <- bd[(nl+1):length(bd)]
        bd1 <- nb
        
        dff <- findInterval(x,bd1,rightmost.closed=TRUE) ##faster than cut
        
        tb <- table(dff,y)
        if(method==1) cacc <- caim(tb)
        if(method==2) cacc <- cacc(tb)
        if(method==3) cacc <- ameva(tb)
        if(cacc>bestC){
            bestC <- cacc
            iBest <- i
            addCut <- di[i]
        }
    }  
    if(!is.na(iBest)) newDi <- di[-iBest] 
    return(list(addCut=addCut,cacc=bestC,newDi=newDi, bd=bd))
}
