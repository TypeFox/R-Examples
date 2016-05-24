"rreslife.lmoms" <- function(f,para, nmom=5) {
    if(! are.par.valid(para)) return()
    
    if(length(f) > 1) {
       warning("f is a vector, only first value will be used")
       f <- f[1]
    }

    if(nmom < 1) {
       warning("nmom is < 1, which is meaningless")
       return()
    }

    L <- R <- vector(mode="numeric", length=nmom)
    
    "afunc" <- function(p, u=NULL, r=NULL, k=NULL) {
          C <- u
          A <- (p/u)^(r-k-1)
          B <- (1 - p/u)^k
          Qp <- par2qua(p,para,paracheck=FALSE)
          v <- Qp*A*B/C
          #cat(c("U:",   u,  "\n"))
          #cat(c("B:",   B,  "\n"))
          #cat(c("A:",   A,  "\n"))
          #cat(c("Qp:", Qp,  "\n"))
          #cat(c("v:", v,  "\n"))
          return(v)
    }

    for(r in 1:nmom) {
       if(f == 0) { 
          if(r == 1) {
             L[1] <- par2qua(0,para,paracheck=FALSE)
          } else {
             L[r] <- NA
          }
          next
       }
       vals <- sapply(0:(r-1),
       function(k) {
          tmp <- NULL
          try(tmp <- integrate(afunc, lower=0, upper=f, u=f, r=r, k=k))
          INTv <- ifelse(is.null(tmp), NA, tmp$value)
          A <- (-1)^k
          B <- choose(r - 1, k)^2
          return(A*B*INTv)
       })
       L[r] <- sum(vals)
    }

    R[1] <- NA
    if(nmom > 1) {
       R[2] <- L[2]/L[1]
    }
    if(nmom >= 3) { 
       R[3:nmom] <- sapply(3:nmom, function(r) { return(L[r]/L[2]) })
    }
    z <- list(lambdas=L, ratios=R,
              life.notexceeds=par2qua(f,para,paracheck=FALSE),
              life.percentile=100*(f),
              trim=NULL, leftrim=NULL, rightrim=NULL,
              source="lmoms.rreslife")
    return(z)
}


