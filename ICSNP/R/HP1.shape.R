`HP1.shape` <-
function(X,location="Estimate",na.action=na.fail,...)
    {
    X<-na.action(X)
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric") 
    X<-as.matrix(X)
    
    n<-dim(X)[1]
    p<-dim(X)[2]

    if (p<2) stop("'X' must be at least bivariate")  

    if (!is.numeric(location))
     {
     loc<-match.arg(location,c("Estimate","Origin"))
     spat.signs<-switch(location,
                       "Estimate"=spatial.sign(X),
                       "Origin"=spatial.sign(X,center=F))
     }
    else 
     {
     if (length(location)!=p) stop("'location' is of wrong dimension")
     spat.signs<-spatial.sign(X,center=location)
     }
    V<-attr(spat.signs,"shape")
    center<-attr(spat.signs,"center")
    V.sqrt<-mat.sqrt(V)
    radius<-sqrt(mahalanobis(X,center,V))
    scores.radius<-qchisq(rank(radius)/(n+1),p)
    W<-V.sqrt %*% ((1/n)*(t(scores.radius*spat.signs) %*% spat.signs)) %*% V.sqrt
    colnames(W) <- colnames(X)
    rownames(W) <- colnames(X)
    W/det(W)^(1/p)
    }
