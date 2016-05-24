performpcls<-function(x, y, nKnots = 6){

    f.ug <- gam(y ~ s(x, k = nKnots, bs = "cr"),family=binomial())
    dat <- data.frame(x = x, y = y)
    sm <- smoothCon(s(x, k = nKnots, bs = "cr"), dat, knots = NULL)[[1]]
    if (length(sm$xp) < nKnots) 
        warning("Few than 6 nKnots were specified!\n Possible inaccurate fitting!\n")
    F <- mono.con(sm$xp,TRUE)
     
    G<-list(X=sm$X,C=matrix(0,0,0),sp=f.ug$sp,p=sm$xp,y=y,w=y*0+1)
    G$Ain<-F$A;G$bin<-F$b;G$S<-sm$S;G$off<-0
    p <- pcls(G)
    fv <- Predict.matrix(sm, data.frame(x = x)) %*% p
    fv <- as.vector(fv)
    #If pcls still fails for some reason and returns negative values (or values that are so low that they'll have an effective relative weight of Inf
    #then we need to add a little bit to every weight to make things non-negative, the choice of minimum weight is somewhat arbitrary and serves only to
    #ensure that we don't have positive weights that will be infinitely prefered as a ratio.
    lower_bound=10^-3
    if(min(fv)<lower_bound)
    	fv=fv-min(fv)+lower_bound
    fv2<-as.matrix(fv)
    return(fv2)

}######function