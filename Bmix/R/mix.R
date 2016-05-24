
"mix" <- function(Z, alpha, g0params, times=NULL,  rho=NULL,  cat=0,
                  state=NULL, read=FALSE, print=FALSE, N=100, niter=0)
{
  if(is.null(state)){ state <- sample(seq(0,999), 3) }

  Z <- as.matrix(Z)
  dim <- ncol(Z)
  Total <- nrow(Z)
  if(cat > dim) cat <- 0
  dim <- dim-cat

    if(cat>0){ levels <- c()
               for(i in 1:cat){
                 levels <- c(levels,length(unique(Z[,dim+i])))
                 if(min(Z[,dim+i] != 0)){ stop("cat min starts at zero")  }
               }}
    else levels = 0

  if(N < 1 || niter < 0)
    { stop("nums parameter is invalid") }

  if(is.null(rho)){ rho <- 0 }
  if(alpha <= 0 || rho > 1 || rho < 0)
    { stop("alpha/rho is invalid") }

  if( length(times) != Total ){ times <- rep(0,Total) }

  if(niter > 10) k <- integer(niter*Total)
  else k <- integer(N*Total)

  out <- .C("mixsample",
            rstate=as.integer(state),
            iopar=as.integer(c(read,print)),
            nums=as.integer(c(N,niter)),
            dims=as.integer(c(dim, cat, levels)),
            Total=as.integer(Total),
            times=as.integer(times),
            Z=as.double(t(Z)),
            params=as.double(c(alpha, rho, g0params)),
            margllhd=double(Total),
            m=double(Total),
            k=k,
            PACKAGE="Bmix")

  dput(bmixrun <- list(dim=dim, cat=cat, levels=levels,
                       g0params=g0params, alpha = alpha, rho = rho),
       file=".bmixrun.robj")

  k <- array(out$k, dim=c(Total,N,niter+1),
                       dimnames=list(paste("o", 1:Total, sep=""),
                         paste("p", 1:N, sep=""),
                         paste("i", 1:(niter+1), sep="")))

  invisible( list(state=state, read=read, print=print, Z=Z,
                  N=N, niter=niter, dim=dim, cat=cat, times=times,
                  levels=levels, Total=Total, g0params=g0params,
                  alpha = alpha, rho = rho, logprob = out$margllhd, m=out$m, k=k) )

}



"particle" <- function(i, mixobj=NULL, t=1, rho=0){
  if(is.null(mixobj)) mixobj <- dget(".bmixrun.robj")
  pmat <- read.table(paste(".particle", i,".", t,".", rho,".txt", sep=""))
  counts <- c()
  if(mixobj$cat > 0){ for(b in 1:mixobj$cat){
    counts <- c(counts, paste("counts",b,1:mixobj$levels[b], sep=".")) }}
  names(pmat)  <- c("n", paste("mean", 1:mixobj$dim,sep="."), paste("S", 1:mixobj$dim^2,sep="."),
                    counts, "p", paste("a",1:mixobj$dim,sep="."), paste("B",1:mixobj$dim^2,sep="."),"c")
  row.names(pmat) <- c(0:(nrow(pmat)-1))
  return(pmat)
}

"getmu" <- function(prt, j){
  if((j+1)>nrow(prt)){  return(print("error: j>m")) }
  else{
    return(as.numeric(prt[(j+1),grep("a.",names(prt),fixed=TRUE)]))
  }
}

"getsigma" <- function(prt, j){
  if((j+1)>nrow(prt)){  return(print("error: j>m")) }
  else{
    sig <- as.numeric(prt[(j+1),grep("B.",names(prt),fixed=TRUE)])
    return(matrix(sig, ncol=sqrt(length(sig))))  }
}

"rwsh" <- function(nu, Omega, state=NULL){
  if(is.null(state)){ state <- sample(seq(0,999), 3) }
  if(ncol(Omega) != nrow(Omega)){ stop("invalid omega") }

  out <- .C("rwish",
            rstate=as.integer(state),
            dim=as.integer(ncol(Omega)),
            nu=as.integer(nu),
            psi=as.double(Omega) )

  return(matrix(out$psi, ncol=ncol(Omega)))
}





