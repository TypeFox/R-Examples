MCMCgrm <- function(model,prior,data,GRM,eps=0,n.thin=10,n.burnin=3000,n.iter=13000,...)
{
  for(p in c("Matrix", "MCMCglmm")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!require(p, quietly = TRUE, character.only=TRUE))
        warning(paste("MCMCgrm needs package `", p, "' to be fully functional; please install", sep=""))
     }
  }
  N <- dim(data)[1]
  GRM <- GRM + diag(eps,N)
  i <- rep(1:N,rep(N,N))
  j <- rep(1:N,N)
  s <- Matrix::spMatrix(N,N,i,j,as.vector(GRM))
  Ginv<-Matrix::solve(s)
  class(Ginv) <- "dgCMatrix"
  rownames(Ginv) <- Ginv@Dimnames[[1]] <- with(data,id)
  m <- MCMCglmm::MCMCglmm(as.formula(model), random=~id, ginverse=list(id=Ginv), data=data,
                prior=prior, thin=n.thin, burnin=n.burnin, nitt=n.iter, ...)
  summary(m)
  return(m)
}

# 8-1-2016 MRC-Epid JHZ
