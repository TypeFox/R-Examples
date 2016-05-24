
LF.gof<- function(X, rmin=NULL, rmax=NULL, na.rm=TRUE){
    
   # controls that is using an envelope object and that it has been computed with the argument savefuns=TRUE
    if(sum(!is.na(match(class(X),"envelope")))==0) stop("LF.gof only accept spatstat envelope objects")
    if(is.null(attributes(X)$simfuns)) stop ("did you compute the envelopes with the argument 'savefuns=TRUE'?") 
    
    # get data
    r <- as.data.frame(attributes(X)$simfuns)[,1]
    H1 <- X$obs
    Hi <- as.data.frame(attributes(X)$simfuns)[,-1]
    
    #Join observed and simulated functions
    H <- cbind(H1,Hi) 
    
    
    # select the range of integrable results
    if (!is.null(rmin) & !is.null(rmax)){ H <- H[r >= rmin & r <= rmax,]}
    if (!is.null(rmin) & is.null(rmax)) { H <- H[r >= rmin,]}
    if (is.null(rmin) & !is.null(rmax)) { H <- H[r <= rmax,]}
    
    # count how may NA values for each distance r
    na.count<-apply(H,1, function(x) sum(is.na(x)))
    na.count.by.r <-na.count[na.count>0] 

    s <- dim(H)[2]

    # compute summary statistics for the observed (i.e., first) pattern 
    Hmeans <- rowSums (H[,-1],  na.rm = na.rm)/(s-1)
    u <- sum((H[,1]-Hmeans)^2, na.rm=T)
     
     # compute summary statistics for the simulated patterns
     for ( i in 2:s) {
         Himeans <- rowSums (H[,-i],na.rm=na.rm)/(s-1) 
         u <- c(u, sum((H[,i]-Himeans)^2, na.rm=na.rm)) 
       }
     p <-  1-((rank(u)[1]-1)/(s))
     
return(list(u=u[1], p=p,na.count.by.r=na.count.by.r))
}
