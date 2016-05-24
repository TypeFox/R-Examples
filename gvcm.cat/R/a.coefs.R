a.coefs <-
function(indices, control, beta, x) {
  
# functions

  a <- function (p = NULL, m = 1, ...) # order(beta) => for high coding, 1:length(beta) => else
      { k <- length(p)
        m <- min(k-1,m)
        if (k>1)
        {
          if (k > 2)
          {
            mat <- matrix(nrow=k, ncol=0)
            for (j in 1:m){
              mat2 <- diag(1, ncol=k, nrow=k)
              for (i in 1:(k-j)){
                mat2[p[i],p[i+j]] <- -1
                }
              mat2 <- mat2[,-p[1:j]]
              mat <- cbind (mat, mat2)
            }
          } else { mat <- matrix(c(-1,1), ncol=1) }  # k == 2
        } else {
        mat <- matrix(ncol=0,nrow=1)                
        }
        colnames(mat) <- c()
        return (mat)
      }

# definitions
  indexs <- as.matrix(indices[c("index1", "index2", "index3", "index4", "index5", "index6", "index7", "index8", "index9"), ])
  index.appro <- colnames(indices)
  ncov <- ncol(indices)
  pairwise <- control$pairwise
  high <- control$high
  
# blank output
  A <- matrix(0,ncol=0,nrow=0)
  phis.v       <- c() # 1 for phi, -1 for 1-phi, 0 for no weighting
  phis.gf      <- c() # grouped.fused
  phis.vs      <- c() # vspline
  phis.sp      <- c() # sp line
  w.categ      <- c()
  w.pairwise   <- c()
  which.a      <- c()
  which.covariate <- c()
  
# pairwise & NULL
  if (pairwise==TRUE) {  
      p. <- function(x){1:length(x)}
      m.n <- function(x){x}
      } else {
      p. <- function(x){order(x)} 
      m.n <- function(x){1}
      }
  if (!is.null(high)) {  
      p. <- function(x){order(x)}
      m.n <- function(x){high}
      }
      
# build A, weights, etc.
for (i in 1:ncov){

  vonbis <- (sum(indexs["index1",1:i])-indexs["index1",i]+1):sum(indexs["index1",1:i])
  ki <- indexs["index1",i]

  if(sum(indexs[-c(1),i])==0){ # covariate i not regularized
    A <- bdiag(A, matrix(0,ncol=0,nrow=indices["index1",i]))   

  } else { # covariate i regularized

    if (indexs["index2",i]!=0) {   # index2 - indicator for special v

          if (indices["index2b",i]==1) { # indicator b_j = 1, abs von coefs penalisiert!

              if (indices["index2",i]<0 ) Ai <- a(p.(beta[vonbis]), m.n(ki))
              if (indices["index2",i]>0 ) Ai <- a(p.(beta[vonbis]), 1) 
              phis.v <- c(phis.v, rep(1,ncol(Ai)))
    
              Ai     <- cbind(Ai,  diag(ki))
              phis.v <- c(phis.v, rep(-1, ki))
              
              phis.gf <- c(phis.gf, rep(0, ncol(Ai)))
              phis.vs <- c(phis.vs, rep(0, ncol(Ai)))
              phis.sp <- c(phis.sp, rep(0, ncol(Ai)))
                          
              A <- bdiag(A, Ai)
              
              w.categ <- if (indices["index2",i]<0) c(w.categ, rep(4/(1+ki), ncol(Ai))) else
                                  c(w.categ, rep(ki/(ki-.5), ncol(Ai)))
                  
          } else { # indicator b_j = 0, abs von coefs nicht penalisiert!
          
              if (indices["index2",i]<0 ) Ai <- a(p.(beta[vonbis]), m.n(ki))
              if (indices["index2",i]>0 ) Ai <- a(p.(beta[vonbis]), 1) 

              phis.v  <- c(phis.v, rep(1,ncol(Ai)))
              phis.gf <- c(phis.gf, rep(0, ncol(Ai)))
              phis.vs <- c(phis.vs, rep(0, ncol(Ai)))
              phis.sp <- c(phis.sp, rep(0, ncol(Ai)))
                          
              A <- bdiag(A, Ai)
              
              w.categ <- if (indices["index2",i]<0) c(w.categ, rep(4/(ki-1), ncol(Ai))) else c(w.categ, rep(2*ki/(ki-1), ncol(Ai)))
          
          }    

          w.pairwise <- if (!pairwise) {
                rang <- rank(beta[vonbis], ties.method="first")  
                w.p.rank <-  rang*(ki-rang)  # without "+ 1" (ki - 1 difference)
                w.p.order <- c()
                for (j in 1:(ki-1)) {
                     w.p.order <- c(w.p.order, which(Ai[,j]==-1))
                }
                c(w.pairwise, w.p.rank[w.p.order], rep(1, ncol(Ai)-length(w.p.order))) # ordered according to differences in Ai
          } else {c(w.pairwise, rep(1, ncol(Ai)))}
          
          which.a <- c(which.a, rep(index.appro[i], ncol(Ai)))
          which.covariate <- c(which.covariate, rep(i, ncol(Ai)))      
          
        }
        
    if (indexs["index3",i]!=0) { # indicator for special p
          bd <- c(beta[vonbis], 0)
    
          if (indices["index3",i]<0) Ai <- a(p.(bd), m.n(ki+1))
          if (indices["index3",i]>0) Ai <- a(p.(bd), 1)               
                                 
          phis.v  <- c(phis.v,  rep(0, ncol(Ai)))
          phis.gf <- c(phis.gf, rep(0, ncol(Ai)))
          phis.vs <- c(phis.vs, rep(0, ncol(Ai)))
          phis.sp <- c(phis.sp, rep(0, ncol(Ai)))
          
          w.categ <- if (indices["index3",i]<0) c(w.categ, rep(2/(ki+1), ncol(Ai))) else
                             c(w.categ, rep(1, ncol(Ai)))
                                                                                                                                            
          w.pairwise <- if (!pairwise) {
                rang <- rank(bd, ties.method="first")  
                w.p.rank <-  rang*(ki-rang+1)  # ki differences
                w.p.order <- c()
                for (j in 1:ki) {
                     w.p.order <- c(w.p.order, which(Ai[,j]==-1))
                }
                c(w.pairwise, w.p.rank[w.p.order]) # ordered according to differences in Ai
          } else {c(w.pairwise, rep(1, ncol(Ai)))}
          
#          Ai <- -1 * Ai
#          A  <- bdiag(A, Ai[1:ki, ])
          A  <- bdiag(A, Ai[-p.(bd)[1], ]) # -1 # 2:(ki+1)  
          
          which.a <- c(which.a, rep(index.appro[i], ncol(Ai)))
          which.covariate <- c(which.covariate, rep(i, ncol(Ai)))      
          
        }
        
    if (indexs["index4",i]!=0) { # indicator for grouped 

          #Ai <- diag(ki) # hier auch differenzen matrix moeglich!!
          #A <- bdiag(A, Ai)
          #phis.v  <- c(phis.v,  rep(0, ki))
          #phis.gf <- c(phis.gf, rep(0, ki))
          #phis.vs <- c(phis.vs, rep(0, ki))
          #phis.sp <- c(phis.sp, rep(0, ki))
          #w.categ <- c(w.categ, rep(1, ki))    ##
          #w.pairwise <- c(w.pairwise, rep(1, ki)) 
          #which.a <- c(which.a, rep(paste("grouped", i, sep=""), ki))
          #which.covariate <- c(which.covariate, rep(i, ki))

          if (!control$grouped.cat.diffs) {
               Ai <- diag(ki) 
          } else {
              bd <- c(beta[vonbis], 0)
              if (indices["index4",i]<0) Ai <- a(p.(bd), m.n(ki+1))
              if (indices["index4",i]>0) Ai <- a(p.(bd), 1)               
#              Ai <- -1 * Ai
#              Ai <- Ai[1:ki, ]
              Ai <- Ai[-p.(bd)[1], ]      
          }      
          A <- bdiag(A, Ai)

          phis.v  <- c(phis.v,  rep(0, ncol(Ai)))
          phis.gf <- c(phis.gf, rep(0, ncol(Ai)))
          phis.vs <- c(phis.vs, rep(0, ncol(Ai)))
          phis.sp <- c(phis.sp, rep(0, ncol(Ai)))

          w.categ <- c(w.categ, rep(1, ncol(Ai))) 
#                     if (indices["index3",i]<0 && !control$grouped.cat.diffs) {
#                            c(w.categ, rep(2/(ki+1), ncol(Ai)))
#                          } else {
#                            c(w.categ, rep(1, ncol(Ai)))
#                          }

          w.pairwise <- c(w.pairwise, rep(1, ncol(Ai))) 
          which.a <- c(which.a, rep(paste("grouped", i, sep=""), ncol(Ai)))
          which.covariate <- c(which.covariate, rep(i, ncol(Ai)))
    
        } 
        
    if (indexs["index5",i]!=0) { # indicator for grouped.fused, fuer stephi
    
          if (indices["index5",i]<0 ) Ai <- a(p.(beta[vonbis]), m.n(ki))
          if (indices["index5",i]>0 ) Ai <- a(p.(beta[vonbis]), 1) 

          A <- bdiag(A, cbind(Ai, diag(ki)))       

          phis.v  <- c(phis.v,  rep(0, ncol(Ai)+ki) )
          phis.gf <- c(phis.gf, rep(1,ncol(Ai)), rep(-1, ki))
          phis.vs <- c(phis.vs, rep(0, ncol(Ai)+ki) )
          phis.sp <- c(phis.sp, rep(0, ncol(Ai)+ki) )

          w.categ <- if (indices["index5",i]<0) c(w.categ, rep(4/(ki+1), ncol(Ai)+ki)) else
                              c(w.categ, rep(ki/(ki-.5), ncol(Ai)+ki))

          w.pairwise <- if (!pairwise) {
                rang <- rank(beta[vonbis], ties.method="first")  
                w.p.rank <-  rang*(ki-rang)  # without "+ 1" (ki - 1 difference)
                w.p.order <- c()
                for (j in 1:(ki-1)) {
                     w.p.order <- c(w.p.order, which(Ai[,j]==-1))
                }
                c(w.pairwise, w.p.rank[w.p.order], rep(1, ki)) # ordered according to differences in Ai
          } else {c(w.pairwise, rep(1, ncol(Ai)+ki))}

          which.a <- c(which.a, rep("L1", ncol(Ai)), rep(paste("grouped", i, sep=""), ki))
          which.covariate <- c(which.covariate, rep(i, ncol(Ai)+ki))      

        } 
        
    if (indexs["index6",i]==1) { # indicator for sp
          if (colnames(indices)[i]=="L2"){ # transformed B-Splines
              Ai <- diag(ki-1)
              A <- bdiag(A, matrix(0,ncol=0,nrow=1), Ai)
              phis.v  <- c(phis.v,  rep(0, ki-1))
              phis.gf <- c(phis.gf, rep(0, ki-1))
              phis.vs <- c(phis.vs, rep(0, ki-1))
              phis.sp <- c(phis.sp, rep(0, ki-1))
              w.categ <- c(w.categ, rep(1, ki-1))  
              w.pairwise <- c(w.pairwise, rep(1, ki-1))
              which.a <- c(which.a, rep(paste("L2", sep=""), ki-1))
              which.covariate <- c(which.covariate, rep(i, ki-1))
          } else { # transformed B-Splines, L1 on slope, grouped on deviations from linear trend
              Ai <- diag(ki) 
              A <- bdiag(A, Ai) 
              phis.v  <- c(phis.v,  rep(0, ki))
              phis.gf <- c(phis.gf, rep(0, ki))
              phis.vs <- c(phis.vs, rep(0, ki))
              phis.sp <- c(phis.sp, 1, rep(-1, ki-1)) # 1*sp on slope, (1-sp) on deviations from slope
              w.categ <- c(w.categ, rep(1, ki))  
              w.pairwise <- c(w.pairwise, rep(1, ki))
              which.a <- c(which.a, "L1", rep(paste("grouped", i, sep=""), ki-1))
              which.covariate <- c(which.covariate, rep(i, ki))
          }          
        }

    if (indexs["index7",i]!=0) { # indicator for special SCAD
    
          bd <- c(beta[vonbis], 0)
    
          if (indices["index7",i]<0) Ai <- a(p.(bd), m.n(ki+1))
          if (indices["index7",i]>0) Ai <- a(p.(bd), 1)               
                                 
          phis.v  <- c(phis.v,  rep(0, ncol(Ai)))
          phis.gf <- c(phis.gf, rep(0, ncol(Ai)))
          phis.vs <- c(phis.vs, rep(0, ncol(Ai)))
          phis.sp <- c(phis.sp, rep(0, ncol(Ai)))
          
          w.categ <- if (indices["index3",i]<0) c(w.categ, rep(2/(ki+1), ncol(Ai))) else
                             c(w.categ, rep(1, ncol(Ai)))
          
          w.pairwise <- if (!pairwise) {
                rang <- rank(bd, ties.method="first")  
                w.p.rank <-  rang*(ki-rang+1)  # ki differences
                w.p.order <- c()
                for (j in 1:ki) {
                     w.p.order <- c(w.p.order, which(Ai[,j]==-1))
                }
                c(w.pairwise, w.p.rank[w.p.order]) # ordered according to differences in Ai
          } else {c(w.pairwise, rep(1, ncol(Ai)))}
          
#          Ai <- -1 * Ai
#          A  <- bdiag(A, Ai[1:ki, ])
          A  <- bdiag(A, Ai[-p.(bd)[1], ])
          
          which.a <- c(which.a, rep(index.appro[i], ncol(Ai)))
          which.covariate <- c(which.covariate, rep(i, ncol(Ai)))      
          
        }

    if (indexs["index8",i]!=0) { # indicator for special elastic
          bd <- c(beta[vonbis], 0)
    
          if (indices["index8",i]<0) Ai <- a(p.(bd), m.n(ki+1))
          if (indices["index8",i]>0) Ai <- a(p.(bd), 1)               
                                 
          phis.v  <- c(phis.v,  rep(0, ncol(Ai)))
          phis.gf <- c(phis.gf, rep(0, ncol(Ai)))
          phis.vs <- c(phis.vs, rep(0, ncol(Ai)))
          phis.sp <- c(phis.sp, rep(0, ncol(Ai)))
          
          w.categ <- if (indices["index3",i]<0) c(w.categ, rep(2/(ki+1), ncol(Ai))) else
                             c(w.categ, rep(1, ncol(Ai)))
          
          w.pairwise <- if (!pairwise) {
                rang <- rank(bd, ties.method="first")  
                w.p.rank <-  rang*(ki-rang+1)  # ki differences
                w.p.order <- c()
                for (j in 1:ki) {
                     w.p.order <- c(w.p.order, which(Ai[,j]==-1))
                }
                c(w.pairwise, w.p.rank[w.p.order]) # ordered according to differences in Ai
          } else {c(w.pairwise, rep(1, ncol(Ai)))}
          
#          Ai <- -1 * Ai
#          A  <- bdiag(A, Ai[1:ki, ])
          A  <- bdiag(A, Ai[-p.(bd)[1], ])
          
          which.a <- c(which.a, rep(index.appro[i], ncol(Ai)))
          which.covariate <- c(which.covariate, rep(i, ncol(Ai)))      
          
        }

    if (indexs["index9",i]!=0) {   # indicator for special vspline
          
          levels.u <- indices["index2b",i] # index2b = number levels of u; 
          k.spline <- ki/levels.u # k.spline = number of coefficients per spline and level
                                   # ki = number of overall coefficients
          k.nl <- k.spline - 1 # number of coefficients per spline and level  for non-linear deviation
          
          # A
          A.diffs <- if (indices["index9",i]<0 ) { a(1:levels.u, levels.u-1) } else { a(1:levels.u, 1) }
              #if (indices["index9",i]<0 ) a(1:levels.u, levels.u-1)
              #if (indices["index9",i]>0 ) a(1:levels.u, 1)
          n.diffs <- ncol(A.diffs)
          colnames(A.diffs) <- as.numeric(1:n.diffs)
          dia <- diag(k.spline)
          colnames(dia) <- rep(paste("grouped:", i, sep=""), k.spline)
          A.diffs <- kronecker(A.diffs, dia, make.dimnames = TRUE)
          
          L1abs <- (0:(levels.u-1))*k.spline + 1
          L1diffs <- (0:(n.diffs-1))*k.spline + 1

          A.smooth <- diag(ki)[,-L1abs]
          A.slope.abs <- diag(ki)[,L1abs]
          A.slope.diffs <- A.diffs[, L1diffs]
          A.grouped <- A.diffs[, -L1diffs]

          A <- bdiag(A, cbind(A.smooth, A.slope.abs, A.slope.diffs, A.grouped))
          
          # factors
          phis.v  <- c(phis.v,  rep(0, ki-levels.u), rep(-1, levels.u), rep(1, ncol(A.slope.diffs)), rep(0, ncol(A.grouped)))
          phis.gf <- c(phis.gf, rep(0, ki+ncol(A.diffs)))
          phis.vs <- c(phis.vs, rep(-1, ki-levels.u), rep(0, levels.u), rep(0, ncol(A.slope.diffs)), rep(1, ncol(A.grouped)))
          phis.sp <- c(phis.sp, rep(0, ki+ncol(A.diffs)))
#          phis.v  <- c(phis.v,  rep(0, ki-levels.u), rep(-1, levels.u), rep(1, n.diffs), rep(0, ki-levels.u))
#          phis.gf <- c(phis.gf, rep(0, 2*ki))
#          phis.vs <- c(phis.vs, rep(-1, ki-levels.u), rep(0, levels.u + n.diffs), rep(1, ki-levels.u))
#          phis.sp <- c(phis.sp, rep(0, 2*ki))
          w.categ <- if (indices["index9",i]<0 ) {
                     c(w.categ, rep(2*k.nl*(1/levels.u + .5*(levels.u-1)), ki-levels.u), rep(4/(levels.u+1), levels.u),         rep(4/(levels.u+1), n.diffs),         rep(2*k.nl*(1/levels.u + .5*(levels.u-1)), ki-levels.u))
                     } else {
                     c(w.categ, rep(2*k.nl, ki-levels.u),                                rep(levels.u/(levels.u-.5), levels.u), rep(levels.u/(levels.u-.5), n.diffs), rep(2*k.nl, ki-levels.u))                     
                     }
          w.pairwise <- c(w.pairwise, rep(1, ki+ncol(A.diffs) )) # rep(1, 2*ki)
          which.a <- c(which.a, rep(paste("grouped", i , 1:levels.u, sep="."), each=k.spline-1), 
                   rep("L1", levels.u + n.diffs), colnames(A.grouped))
          which.covariate <- c(which.covariate, rep(i, ki+ncol(A.diffs))) # rep(i, 2*ki)

        }
 
  } # else 

} # for

A <- as.matrix(A)

# vector fused wie in Fixed effects
inona <- indices[1,]; names(inona) <- c() 
if (control$subjspec.gr==TRUE && all.equal(inona, rep(indices[1,1], ncol(indices)))){
    nr <- length(which.a)/ncol(indices)
    which.a <- rep(paste("grouped.", 1:nr, sep=""), ncol(indices)) 
}

# adaptive weights
weighted.which <- c("index2", "index3", "index5", "index7", "index8")
if (control$subjspec.gr==TRUE && all.equal(inona, rep(indices[1,1], ncol(indices)))) { # v raus falls mit grouped...
    weighted.which <- c("index3", "index5", "index7", "index8")
}
weighted.covs <- which(colSums(as.matrix(indices[weighted.which,]))!=0) # v, p, grouped.fused, SCAD, elastic
weighted <- which(which.covariate %in% weighted.covs)
w.ada <- w.categories <- w.cases <- rep(1, ncol(A)) 

if (control$adapted.weights) {
            adaptive.weights <- abs(t(A)%*%beta)
            if (any(adaptive.weights==0)) adaptive.weights[which(adaptive.weights==0)] <- 10^(-control$digits)
            adaptive.weights <- as.vector(adaptive.weights^(-1)) # as.vector(abs(t(A)%*%beta)^(-1))
            lengths <- rle(which.covariate)[[1]]
            ends <- cumsum(lengths)
            if (control$adapted.weights.adj){
                adaptive.overall <- c()
                for (i in 1:length(ends)){
                     adaptive.overall <- c(adaptive.overall, sum(adaptive.weights[c(ends[i]-lengths[i]+1, ends[i])]))
                     } 
                } else {
                adaptive.overall <- rep(1, length(ends))
                }    
            adaptive.adj <- rep(1/adaptive.overall, times=lengths)
            w.ada[weighted]  <- (adaptive.adj * adaptive.weights)[weighted]
}

if (control$level.control)   {w.categories <- w.categ}
if (control$case.control)    {w.cases[weighted] <- (sqrt(abs(t(A))%*%as.matrix(colSums(x!=0))/nrow(x)))[weighted]}

which.continuous <- as.numeric(colSums(x!=0)==nrow(x))
which.continuous[1] <- 0
which.continuous <- as.vector(t(A) %*% which.continuous)

# check?!
if (nrow(A) != ncol(x))
        stop("Inconsistency in arguments!")

# output
output <- list(
  A = A,
  w.adaptive = w.ada,
  w.cases = w.cases,
  w.categories = w.categories,
  w.pairwise = w.pairwise,
  continuous = which.continuous,
  which.a = which.a,
  which.covariate = which.covariate,
  phis.v = phis.v,
  phis.gf = phis.gf,
  phis.vs = phis.vs,
  phis.sp = phis.sp
  )

  return(output)
}
