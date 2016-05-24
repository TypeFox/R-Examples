vimpint <- function(obj,nseq.partial=51,nseq.inter=20,trace=FALSE,percentage=TRUE) #obj is object of boost function
{
dat1.glm <- obj$res.fit
x <- obj$x
degree <- obj$degree
learner <- obj$learner
coef.bj <- obj$coef.bj
interactions <- NULL
    d <- ncol(x)
    if(learner=="linear.regression" && degree==2) d <- -0.5+0.5*sqrt(1+8*d)  #p+p*(p-1)/2 = d for main + interaction terms
    vim <- rep(0,d)
    #construct a one-to-one corespondence table between input column and the ensemble index
    ind <- matrix(NA,ncol=2,nrow=d+d*(d-1)/2)
    kk <- 1
    for(i in 1:d)
     for(j in i:d){
      ind[kk,1] <- i; ind[kk,2] <- j
      kk <- kk + 1
}
    if(learner=="tree"){
     vim <- summary(dat1.glm,plotit=FALSE,order=FALSE)[,2]
     interactions <- vim.interactions(dat1.glm,pred.data=x,x.pair=subset(ind,ind[,1]!=ind[,2]),learner="tree",verbose=trace)
}
    if(learner=="mars"){
    stop("There is a problem here\n")
    vim <- evimp(update(dat1.glm),sqrt.=TRUE,trim=FALSE)[,c(1,4)]
    vim <- vim[order(vim[,1]),]
    vim <- vim[,2]
    interactions <- vim.interactions(dat1.glm,pred.data=x,x.pair=subset(ind,ind[,1]!=ind[,2]),learner="mars",verbose=trace)
}
    if(learner=="smooth.spline"){
     fpartial <- gamplot.out(dat1.glm)
     for(j in 1:d)
       vim[j] <- rif(x=fpartial$xp[,j],y=fpartial$yp[,j])
    }
    if(learner=="linear.regression" && degree!=2)
    vim <- (abs(coef.bj)*apply(x,2,sd))[-1] #standardized coefficients without intercept
    if((learner=="linear.regression" && degree==2) || learner=="tensor.product"){
    if(learner=="linear.regression"){
    xselect <- unique(dat1.glm$ensemble[,"xselect"])
    xselect <- subset(xselect,xselect!=1) - 1 #the first column is for intercept, so "2" means x[,1] instead
    } else xselect <- unique(dat1.glm$ensemble)
    x.part <- unique(as.vector(ind[xselect,])) #for partial dependence
    x.pair <- ind[xselect,] #for paired interactions
    x.pair <- subset(x.pair,x.pair[,1]!=x.pair[,2]) #for paired interactions
    fpartial <- vector("list",d)
    cat("\nlearner is",learner,"\n")
    for (i in x.part){
     if(learner=="linear.regression")
      fpartial[[i]] <- partialPlot(dat1.glm,pred.data=x[,-1],x.var=i,plot=FALSE,learner=learner) #remove the first column 1s, for intercept
     else fpartial[[i]] <- partialPlot(dat1.glm,pred.data=x,x.var=i,plot=FALSE,learner=learner,nseq=nseq.partial)
       vim[i] <- rif(x=fpartial[[i]]$x,y=fpartial[[i]]$y)
      if(trace) cat(i,"",vim[i],"\n")
    }
    interactions <- vim.interactions(dat1.glm,pred.data=x,x.pair=x.pair,learner=learner,verbose=trace,nseq=nseq.inter);
}
   if(percentage) vim <- round(100*vim/sum(vim),2)
    list(vim=vim,interactions=interactions)
}
