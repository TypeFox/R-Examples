rne<-function (y, env, geodist,fun, neighbours, subsets = c(1,0.75, 0.5, 0.25, 0.1),...) {
  dots<-list(...)
  if (inherits(geodist, "dist")) 
      geodist = as.matrix(geodist)
  rne <- list()
  N = nrow(y)
  rne$random <- t(sapply(subsets, function(ss) {
    print(paste("random subset = ", ss))
    r2 <- replicate(10, {  
      est<-sapply(1:N, function(n) {
        retain <- sample(1:(N - 1), size = round((N - 1) * (ss)))
        y2<-y[-n,][retain,]
        keepcols<-colSums(y2)!=0
        mod<-do.call("fun", c(list(y=y2[,keepcols],x=env[-n][retain]),dots))
        predict(mod,y[n,keepcols,drop=FALSE])$fit
      })      
      if(is.null(dim(est))){cor(est, env)^2}
      else{ apply(est,1,cor,env)^2}
    })
    if(is.null(dim(r2))){r2<-mean(r2)} #average across replicates
    else{r2<-rowMeans(r2)}             #average across replicates
    print(r2)
    c(prop = ss, r2 = r2)
  }))
  print(rne$random)
  rne$neighbour<-lapply(neighbours, function(neighbour) {
    print(paste("neighbourhood = ", neighbour,"km"))
    en <- sapply(1:nrow(y), function(n) {
        sum(geodist[n, ] >= neighbour)
    })
    effn = (nrow(y) - mean(en))/(nrow(y) - 1)
    hb <- sapply(1:nrow(y), function(n) {
        y1<-y[-n,]
        env1<-env[-n]
        exneigh <- geodist[n, -n] >= neighbour
        y2<-y1[exneigh,]
        keepcols<-colSums(y2)>0
        mod<-do.call("fun", c(list(y=y2[,keepcols],x=env1[exneigh]),dots))
        predict(mod,y[n,keepcols,drop=F])$fit        
    })
    if(is.null(dim(hb))){hbr<-cor(hb, env)^2}
    else{ hbr<-apply(hb,1,cor,env)^2}

    eb <- sapply(1:nrow(y), function(n) {
        y1<-y[-n,]
        env1<-env[-n]
        neigh <- which(rank(-abs(env1 - env[n]), ties.method = "random") <= (en[n]))
        y2<-y1[neigh,]
        keepcols<-colSums(y2)!=0
        mod<-do.call("fun",c(list(y=y2[,keepcols],x=env1[neigh]),dots))
        predict(mod,y[n,keepcols,drop=FALSE])$fit        
    })
    if(is.null(dim(eb))){ebr<-cor(eb, env)^2}
    else{ ebr<-apply(eb,1,cor,env)^2}
    
    list(neighbour = neighbour, effn = effn, hb.r2 =hbr, eb.r2 = ebr)
  })

  class(rne) <- "RNE"
  return(rne)
}
