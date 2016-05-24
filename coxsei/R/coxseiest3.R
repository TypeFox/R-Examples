coxseiest3 <-
function(dat,par.init,m=2,mit=1000,tr=TRUE,method="L-BFGS-B",
         lower=c(rep(-Inf,ncol(dat)-3),-Inf,0),
         upper=rep(Inf,ncol(dat)-3 + 2)){
  ids <- unique(dat$id) # group ids
  ng <- length(ids) # no. of groups
  ## group sizes (no. of events in each group + 1 (censoring)):
  gs <- as.numeric(table(dat$id))#sapply(ids, function(i){sum(dat$id==ids[i])}) 
  ## offset to locate the starting position of each groups of data
  gofst <- cumsum(gs)-gs
  ## e.g.: dat[gofst[2]+1:gs[2],] contains data for group/individual 2,
  ## which is same as dat[dat$id==ids[2],] contains data for individual 2
  ## same as dat[dat$id]
  C <- dat$Y[!dat$delta] # censoring times
  Z <- as.matrix(dat[,setdiff(colnames(dat),c("Y","delta","id"))],
                 nrow=nrow(dat)) 
  Zs <- array(0,dim=c(dim(Z),ng)) ## we want
  ##Z[gofst+1:gs[i]&dat$id==ids[i], l] gives Z_l(Y_{i1}, ..., Y_{iN_i})
##   for(l in 1:ng){
##     for(i in 1:ng){
##       for(j in 1:NCOL(Z))
##         Zs[gofst[i]+1:gs[i],j,l] <- approxfun(dat$Y[gofst[l]+1:gs[l]],
##              Z[gofst[l]+1:gs[l],j],method="constant",
##              f=1,rule=2)(dat$Y[gofst[i]+1:gs[i]])
##     }
##   }
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

  loglik <- function(param){
    .Call("ll2",dat$Y, as.double(as.vector(Z)), as.double(as.vector(Zs)),
          as.double(C),  as.integer(gs), as.integer(gofst), as.double(param),
          as.integer(m), PACKAGE="coxsei")
  }
  ret <- optim(par.init,loglik,control=list(trace=tr,maxit=mit),
               hessian=TRUE,method=method,
               lower=lower,upper=upper)
  ret
}

