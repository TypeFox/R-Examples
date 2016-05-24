coxseiInt <- function(dat,parest,hessian=NULL,vcovmat=solve(hessian),m=2,
                       gfun=function(x,pa){
                         ifelse(x <= 0, 0,
                                pa[1]*pa[2]*exp(-pa[2]*x))
                       },
                       gfungrd=function(x,pa){
                         if(length(x)==0)return(matrix(0,2,0));
                         rbind(pa[2]*exp(-pa[2]*x),
                               pa[1]*exp(-pa[2]*x)*(1-pa[2]*x)
                               )
                       }
                       ){
  var.warn <- is.null(hessian)&&is.null(vcovmat) ||
  !is.null(hessian)&&(any(dim(hessian)!=length(parest))) ||
  !is.null(vcovmat)&&(any(dim(vcovmat)!=length(parest)))
  if(var.warn)
    warning("Both hessian & vcovmat missing or dimension is wrong!\n");

  ids <- unique(dat$id)
  ng <- length(ids)
  ##   gs <- sapply(ids, function(i){sum(dat$id==ids[i])})
  gs <- as.numeric(table(dat$id))
  gofst <- cumsum(gs)-gs
  C <- dat$Y[!dat$delta]
  Z <- as.matrix(dat[,setdiff(colnames(dat),c("Y","delta","id"))],
                 nrow=nrow(dat)) 
  Zs <- array(0,dim=c(dim(Z),ng)) ## we want
  for(l in 1:ng){
    for(i in 1:ng){
      for(j in 1:NCOL(Z)){
        Zs[dat$id==ids[i],j,l] <- if(diff(range(Z[dat$id==ids[l],j]))==0){
          Z[dat$id==ids[l],j][1]
        }else{
          approxfun(dat$Y[dat$id==ids[l]],
                      Z[dat$id==ids[l],j],method="constant",
                      f=1,rule=2)(dat$Y[dat$id==ids[i]])
        }
      }
    }
  }

  np <- length(parest);
  par <- parest[1:(np-2)];
  parg <- parest[np- 1:0];

  ##   this gives the jumpsizes of \hat U(t) at Y_{i1}, \dots,
  ##   Y_{i,N_i} and the jumptizes of the jump component in its
  ##   variance estimator, i.e. R1/R0^2
  js <- function(i){
    if(gs[i]==1)return(numeric(0))
    res <- numeric(gs[i]-1) # for R0
    res1 <- matrix(0,nrow=np,ncol=gs[i]-1);# for R1/R0^2
    for(j in 1:(gs[i]-1) ){
      posij <- gofst[i]+j;
      for(l in (1:ng)[dat$Y[posij]<=C]){
        yl <- dat$Y[dat$id==l & (dat$delta==1)]
        R0part <- exp(Zs[posij,,l]%*%par +
                      sum(gfun(dat$Y[posij]-
                               tail(yl[yl<dat$Y[posij]],m),
                               parg)))
        res[j] <- res[j]+R0part;
        res1[,j] <- res1[,j]+
          c(Zs[posij,,l], ##separated paresteters, so not summation
            rowSums(gfungrd(dat$Y[posij]-tail(yl[yl<dat$Y[posij]],m),
                            parg))
            ) * R0part;
      }
    }
    res1 <- apply(res1,1,"/",y=res^2);
    dim(res1) <- c(gs[i]-1,np); # now it's high and thin
    res1 <- t(res1); # convert res1 to a matrix of gs[i]-1 columns
    rbind(1/res,res1) # size = (np+1) x (gs[i]-1)
  }
  x <- dat$Y[dat$delta==1];
  ox <- order(x);
  jmps <- unlist(lapply((1:ng)[gs>1],js))
  lj <- length(jmps)
  jo <- jmps[seq(1,lj,by=1+np)][ox]
  varestjo <- jmps[-seq(1,lj,by=1+np)]
  dim(varestjo) <- c(np,length(x))
  varestjo <- varestjo[,ox]
  varest <- apply(varestjo,1,cumsum)
  dim(varest) <- c(length(x),np)
  varest <- if(!is.null(hessian)){
    cumsum(jo^2)+
      apply(varest,1,function(x)x%*%solve(hessian,x))
  }else{
    cumsum(jo^2)+
      apply(varest,1,function(x)x%*% (vcovmat %*% x))
  }
  list(x=x[ox],y=cumsum(jo), varest=varest)
}

## coxseiInt2 <- function(dat,parest,hessian=NULL,vcovmat=solve(hessian),m=2,
##                        gfun=function(x,pa){
##                          ifelse(x <= 0, 0,
##                                 pa[1]*pa[2]*exp(-pa[2]*x))
##                        },
##                        gfungrd=function(x,pa){
##                          if(length(x)==0)return(matrix(0,2,0));
##                          rbind(pa[2]*exp(-pa[2]*x),
##                                pa[1]*exp(-pa[2]*x)*(1-pa[2]*x)
##                                )
##                        }
##                        ){
##   var.warn <- is.null(hessian)&&is.null(vcovmat) ||
##   !is.null(hessian)&&(any(dim(hessian)!=length(parest))) ||
##   !is.null(vcovmat)&&(any(dim(vcovmat)!=length(parest)))
##   if(var.warn)
##     warning("Both hessian & vcovmat missing or dimension is wrong!\n");

##   ##   ret <- numeric(sum(dat$delta))
##   ids <- unique(dat$id)
##   ng <- length(ids)
##   ##   gs <- sapply(ids, function(i){sum(dat$id==ids[i])})
##   gs <- as.numeric(table(dat$id))
##   gofst <- cumsum(gs)-gs
##   C <- dat$Y[!dat$delta]
##   Z <- as.matrix(dat[,setdiff(colnames(dat),c("Y","delta","id"))],
##                  nrow=nrow(dat)) 
##   np <- length(parest);
##   par <- parest[1:(np-2)];
##   parg <- parest[np - 1:0];

##   ## this function gives the [R0(t),R1(t)']':
##   R01 <- function(tt){
##     R0 <- 0;
##     R1 <- numeric(np); # R1 initialized to (0,...,0)'.
##     for(l in (1:ng)[tt <= C]){
##       Zlt <- numeric(np-2);
##       for(j in 1:(np-2)){
##         Zs[dat$id==ids[i],j,l] <- if(diff(range(Z[dat$id==ids[l],j]))==0){
##           Z[dat$id==ids[l],j][1]
##         }else{
##           approxfun(dat$Y[dat$id==ids[l]],
##                       Z[dat$id==ids[l],j],method="constant",
##                       f=1,rule=2)(tt)
##         }
## ##         Zlt[j] <- approxfun(dat$Y[dat$id==ids[l]],Z[dat$id==ids[l],j],
## ##                               method="constant",f=1,rule=2)(tt);
##       }
##       yl <- dat$Y[dat$id==l & dat$delta==1];
##       R0part <- exp(Zlt%*%par + sum(gfun(tt-tail(yl[yl<tt],m), parg)));
##       R1 <- R1 + c(Zlt,rowSums(gfungrd(tt- tail(yl[yl<tt],m),parg)))*R0part;
##       R0 <- R0 + R0part;
##     }
##     c(R0,R1)
##   }
##   x <- dat$Y[dat$delta==1];
##   ox <- order(x);
##   R01mat <- sapply(x,R01)
##   dim(R01mat) <- c(1+np,length(x));
##   R01mat <- R01mat[,ox]
##   Intjmp <- 1/R01mat[1,]
##   IntVjmp <- apply(R01mat[-1,,drop=FALSE],1,"/",y=R01mat[1,]^2)
##   dim(IntVjmp) <- c(length(x),np);
##   ##   IntVjmp <- t(IntVjmp)
##   IntVjmpCum <- apply(IntVjmp,2,cumsum);
##   dim(IntVjmp) <- c(length(x),np);
  
##   varest <- if(!is.null(hessian)){
##     cumsum(Intjmp^2)+
##       apply(IntVjmpCum,1,function(x)x%*%solve(hessian,x))
##   }else{
##     cumsum(Intjmp^2)+
##       apply(IntVjmpCum,1,function(x)x%*% (vcovmat %*% x))
##   }
##   list(x=x[ox],y=cumsum(Intjmp),varest=varest)
## }

