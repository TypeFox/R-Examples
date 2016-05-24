pmatrix <- function(x, FUN, ..., mc.cores = setCores(,reprt=FALSE)){
  #require(parallel)
  cl <- makeCluster(mc <- mc.cores)
  #
  #   parCapply does not pass dimension attributes to FUN !!!!
  #
  z <- parCapply(cl, x, FUN, ...)
  stopCluster(cl)
  z
}

plmatrix <- function(x, FUN, ..., mc.cores = setCores(,reprt=FALSE)){
  dx <- dim(x)[2]
  if(mc.cores>dx) mc.cores <- dx
  n <- trunc((dx-1)/mc.cores)+1
  lx <- list(NULL)
  for(i in 1:(mc.cores-1)) lx[[i]] <- x[,(i-1)*n+1:n]
  lx[[mc.cores]] <- x[,((mc.cores-1)*n+1):dx]
  cl <- makeCluster(mc <- mc.cores)
  lz <- parLapply(cl, lx, FUN , ...)
  stopCluster(cl)
  z <- matrix(0,length(lz[[1]])/n, dx)
  for(i in 1:(mc.cores-1)) z[,(i-1)*n+1:n] <- lz[[i]]
  z[,((mc.cores-1)*n+1):dx] <- lz[[mc.cores]]
  z
}

pnlrdtirg <- function(si,btb,sdcoef,s0ind,ngrad){
  ns0 <- length(s0ind)
  nvox <- length(si)/ngrad
  dim(si) <- c(ngrad,length(si)/ngrad)
  s0 <- if(ns0>1) rep(1/ns0,ns0)%*%si[s0ind,] else si[s0ind,]
  #   ngrad <- dim(btb)[2]
  z <- .Fortran("nlrdtirp",
                as.double(si),
                as.integer(ngrad),
                as.integer(nvox),
                as.double(btb),
                as.double(sdcoef),
                as.double(s0),
                as.integer(200),
                as.double(1e-6),
                res=double((8+ngrad)*nvox),
                as.integer(8+ngrad),
                double(ngrad),
                PACKAGE="dti")$res
  z
}
ptensnl <- function(x,ngrad,btb,sdcoef,maxit=1000,reltol=1e-7){
  nvox <- dim(x)[2]
  matrix(.C("dtens",
            as.integer(nvox),
            param=as.double(x[1:7,]),
            as.double(x[-(1:7),]),
            as.integer(ngrad),
            as.double(btb),
            as.double(sdcoef),
            as.double(rep(0,ngrad)),#si
            as.double(rep(1,ngrad)),#var                         
            as.integer(maxit),#maxit
            as.double(reltol),#reltol
            PACKAGE="dti")$param,7,nvox)
}

pmixtens <- function(x,ngrad0,maxcomp,maxit,pen,grad,reltol,th,penIC,vert){
  nvox <- length(x)/(ngrad0+3+maxcomp)
  dim(x) <- c(ngrad0+3+maxcomp,nvox)
  z <- .C("mixture", 
          as.integer(nvox),
          as.integer(x[(ngrad0+2):(ngrad0+3+maxcomp),]),  
          as.integer(ngrad0),
          as.integer(maxcomp),
          as.integer(maxit),
          as.double(pen),
          as.double(t(grad)),
          as.double(reltol),
          as.double(th),
          as.double(penIC),
          as.double(x[ngrad0+1,]),
          as.double(vert),
          as.double(x[1:ngrad0,]),
          sigma2  = double(nvox),# error variance 
          orient  = double(2*maxcomp*nvox), # phi/theta for all mixture tensors
          order   = integer(nvox),   # selected order of mixture
          lev     = double(2*nvox),         # logarithmic eigenvalues
          mix     = double(maxcomp*nvox),   # mixture weights
          PACKAGE="dti")[c("sigma2","orient","order","lev","mix")]
  rbind(z$order,z$sigma2,matrix(z$lev,2,nvox),matrix(z$mix,maxcomp,nvox),matrix(z$orient,2*maxcomp,nvox))
}
pdkiQP <- function(x,TA,Dmat,Amat){
  ##
  ##  dkiTensor CLLS-QP parallel version
  ##
  #require(quadprog)
  nvox <- dim(x)[2]
  param <- matrix(0,21,nvox)
  for(i in 1:nvox){
    dvec <- -as.vector(t(TA) %*% x[,i])
    resQPsolution <- solve.QP(Dmat, dvec, Amat)$solution
    param[1:6, i] <- resQPsolution[1:6]
    param[7:21, i] <- resQPsolution[7:21] / mean(resQPsolution[1:3])^2
  }
  param
}
pmixtns0 <- function(x,ngrad0,maxcomp,maxit,grad,bv,lambda,alpha,factr,penIC,vert){
  nvox <- length(x)/(ngrad0+1+maxcomp)
  dim(x) <- c(ngrad0+1+maxcomp,nvox)
  z <- .C("mixtrl0", 
          as.integer(nvox),#n1
          as.integer(x[(ngrad0+2):(ngrad0+1+maxcomp),]),#siind 
          as.integer(ngrad0),#ngrad0
          as.integer(maxcomp),#maxcomp
          as.integer(maxit),#maxit
          as.double(t(grad)),#grad_in
          as.double(bv),#bv_in
          as.double(lambda),#lambda_in
          as.double(alpha),#alpha_in
          as.double(factr),#factr
          as.double(penIC),#penIC
          as.double(x[ngrad0+1,]),#sigma2
          as.double(vert),#vert
          as.double(x[1:ngrad0,]),#siq_in
          sigma2  = double(nvox),#sigma2_ret error variance 
          orient  = double(2*maxcomp*nvox),#orient_ret phi/theta for all mixture tensors
          order   = integer(nvox),#order_ret selected order of mixture
          mix     = double(maxcomp*nvox),#mixture weights
          PACKAGE="dti")[c("sigma2","orient","order","mix")]
  lev <- matrix(0,2,nvox)
  lev[1,] <- (alpha+1)*lambda
  lev[2,] <- lambda
  rbind(z$order,z$sigma2,matrix(lev,2,nvox),matrix(z$mix,maxcomp,nvox),matrix(z$orient,2*maxcomp,nvox))
}

pmixtns1 <- function(x,ngrad0,maxcomp,maxit,grad,bv,lambda,alpha,factr,penIC,vert){
  nvox <- length(x)/(ngrad0+1+maxcomp)
  dim(x) <- c(ngrad0+1+maxcomp,nvox)
  z <- .C("mixtrl1", 
          as.integer(nvox),#n1
          as.integer(x[(ngrad0+2):(ngrad0+1+maxcomp),]),#siind 
          as.integer(ngrad0),#ngrad0
          as.integer(maxcomp),#maxcomp
          as.integer(maxit),#maxit
          as.double(t(grad)),#grad_in
          as.double(bv),#bv_in
          as.double(lambda),#lambda_in
          as.double(alpha),#alpha_in
          as.double(factr),#factr
          as.double(penIC),#penIC
          as.double(x[ngrad0+1,]),#sigma2
          as.double(vert),#vert
          as.double(x[1:ngrad0,]),#siq_in
          sigma2  = double(nvox),#sigma2_ret error variance 
          orient  = double(2*maxcomp*nvox),#orient_ret phi/theta for all mixture tensors
          order   = integer(nvox),#order_ret selected order of mixture
          lambda  = double(nvox),#lambda_ret lambda_2 
          mix     = double(maxcomp*nvox),#mixture weights
          PACKAGE="dti")[c("sigma2","orient","order","lambda","mix")]
  lev <- matrix(0,2,nvox)
  lev[1,] <- (alpha+1)*z$lambda
  lev[2,] <- z$lambda
  rbind(z$order,z$sigma2,matrix(lev,2,nvox),matrix(z$mix,maxcomp,nvox),matrix(z$orient,2*maxcomp,nvox))
}

pmixtns2 <- function(x,ngrad0,maxcomp,maxit,grad,bv,lambda,alpha,factr,penIC,vert){
  nvox <- length(x)/(ngrad0+1+maxcomp)
  dim(x) <- c(ngrad0+1+maxcomp,nvox)
  z <- .C("mixtrl2", 
          as.integer(nvox),#n1
          as.integer(x[(ngrad0+2):(ngrad0+1+maxcomp),]),#siind 
          as.integer(ngrad0),#ngrad0
          as.integer(maxcomp),#maxcomp
          as.integer(maxit),#maxit
          as.double(t(grad)),#grad_in
          as.double(bv),#bv_in
          as.double(lambda),#lambda_in
          as.double(alpha),#alpha_in
          as.double(factr),#factr
          as.double(penIC),#penIC
          as.double(x[ngrad0+1,]),#sigma2
          as.double(vert),#vert
          as.double(x[1:ngrad0,]),#siq_in
          sigma2  = double(nvox),#sigma2_ret error variance 
          orient  = double(2*maxcomp*nvox),#orient_ret phi/theta for all mixture tensors
          order   = integer(nvox),#order_ret selected order of mixture
          alpha   = double(nvox),#alpha_ret alpha=(lambda_1-lambda_2)/lambda_2 
          lambda  = double(nvox),#lambda_ret lambda_2 
          mix     = double(maxcomp*nvox),#mixture weights
          PACKAGE="dti")[c("sigma2","orient","order","alpha","lambda","mix")] 
  lev <- matrix(0,2,nvox)
  lev[1,] <- (z$alpha+1)*z$lambda
  lev[2,] <- z$lambda
  rbind(z$order,z$sigma2,lev,matrix(z$mix,maxcomp,nvox),matrix(z$orient,2*maxcomp,nvox))
}
pmixtn0b <- function(x,ngrad,maxcomp,maxit,grad,bv,lambda,alpha,factr,penIC,vert){
  nvox <- length(x)/(ngrad+2+2*maxcomp)
  dim(x) <- c(ngrad+2+2*maxcomp,nvox)
  z <- .C("mixtrl0b", 
          as.integer(nvox),#n1
          as.integer(x[(ngrad+2):(ngrad+1+maxcomp),]),#siind 
          as.integer(x[-(1:(ngrad+1+maxcomp)),]),#wi 
          as.integer(ngrad),#ngrad
          as.integer(maxcomp),#maxcomp
          as.integer(maxit),#maxit
          as.double(t(grad)),#grad_in
          as.double(bv),#bv_in
          as.double(lambda),#lambda_in
          as.double(alpha),#alpha_in
          as.double(factr),#factr
          as.double(penIC),#penIC
          as.double(x[ngrad+1,]),#sigma2
          as.double(vert),#vert
          as.double(x[1:ngrad,]),#siq_in
          sigma2  = double(nvox),#sigma2_ret error variance 
          orient  = double(2*maxcomp*nvox),#orient_ret phi/theta for all mixture tensors
          order   = integer(nvox),#order_ret selected order of mixture
          mix     = double(maxcomp*nvox),#mixture weights
          PACKAGE="dti")[c("sigma2","orient","order","mix")]
  lev <- matrix(0,2,nvox)
  lev[1,] <- (alpha+1)*lambda
  lev[2,] <- lambda
  rbind(z$order,z$sigma2,matrix(lev,2,nvox),matrix(z$mix,maxcomp,nvox),matrix(z$orient,2*maxcomp,nvox))
}

pmixtn1b <- function(x,ngrad,maxcomp,maxit,grad,bv,lambda,alpha,factr,penIC,vert){
  nvox <- length(x)/(ngrad+2+2*maxcomp)
  dim(x) <- c(ngrad+2+2*maxcomp,nvox)
  z <- .C("mixtrl1b", 
          as.integer(nvox),#n1
          as.integer(x[(ngrad+2):(ngrad+1+maxcomp),]),#siind 
          as.integer(x[-(1:(ngrad+1+maxcomp)),]),#wi 
          as.integer(ngrad),#ngrad
          as.integer(maxcomp),#maxcomp
          as.integer(maxit),#maxit
          as.double(t(grad)),#grad_in
          as.double(bv),#bv_in
          as.double(lambda),#lambda_in
          as.double(alpha),#alpha_in
          as.double(factr),#factr
          as.double(penIC),#penIC
          as.double(x[ngrad+1,]),#sigma2
          as.double(vert),#vert
          as.double(x[1:ngrad,]),#siq_in
          sigma2  = double(nvox),#sigma2_ret error variance 
          orient  = double(2*maxcomp*nvox),#orient_ret phi/theta for all mixture tensors
          order   = integer(nvox),#order_ret selected order of mixture
          lambda  = double(nvox),#lambda_ret lambda_2 
          mix     = double(maxcomp*nvox),#mixture weights
          PACKAGE="dti")[c("sigma2","orient","order","lambda","mix")]
  lev <- matrix(0,2,nvox)
  lev[1,] <- (alpha+1)*z$lambda
  lev[2,] <- z$lambda
  rbind(z$order,z$sigma2,matrix(lev,2,nvox),matrix(z$mix,maxcomp,nvox),matrix(z$orient,2*maxcomp,nvox))
}

pmixtn2b <- function(x,ngrad,maxcomp,maxit,grad,bv,lambda,alpha,factr,penIC,vert){
  nvox <- length(x)/(ngrad+2+2*maxcomp)
  dim(x) <- c(ngrad+2+2*maxcomp,nvox)
  z <- .C("mixtrl2b", 
          as.integer(nvox),#n1
          as.integer(x[(ngrad+2):(ngrad+1+maxcomp),]),#siind 
          as.integer(x[-(1:(ngrad+1+maxcomp)),]),#wi 
          as.integer(ngrad),#ngrad
          as.integer(maxcomp),#maxcomp
          as.integer(maxit),#maxit
          as.double(t(grad)),#grad_in
          as.double(bv),#bv_in
          as.double(lambda),#lambda_in
          as.double(alpha),#alpha_in
          as.double(factr),#factr
          as.double(penIC),#penIC
          as.double(x[ngrad+1,]),#sigma2
          as.double(vert),#vert
          as.double(x[1:ngrad,]),#siq_in
          sigma2  = double(nvox),#sigma2_ret error variance 
          orient  = double(2*maxcomp*nvox),#orient_ret phi/theta for all mixture tensors
          order   = integer(nvox),#order_ret selected order of mixture
          alpha   = double(nvox),#alpha_ret alpha=(lambda_1-lambda_2)/lambda_2 
          lambda  = double(nvox),#lambda_ret lambda_2 
          mix     = double(maxcomp*nvox),#mixture weights
          PACKAGE="dti")[c("sigma2","orient","order","alpha","lambda","mix")] 
  lev <- matrix(0,2,nvox)
  lev[1,] <- (z$alpha+1)*z$lambda
  lev[2,] <- z$lambda
  rbind(z$order,z$sigma2,lev,matrix(z$mix,maxcomp,nvox),matrix(z$orient,2*maxcomp,nvox))
}
pgetsii30 <- function(x,maxcomp,dgrad,th,isample0,nsi,nth,nvico,nguess){
  nvox <- length(x)/(nsi+2)
  dim(x) <- c(nsi+2,nvox)
  z <- .Fortran("pgtsii30",
                as.double(x[1:nsi,]),
                as.double(x[nsi+1,]),#sigma2
                as.integer(nsi),
                as.integer(nvox),
                as.integer(maxcomp),
                as.double(dgrad),
                as.integer(nvico),
                as.double(th),
                as.integer(nth),
                as.integer(x[nsi+2,]),#indth
                double(nsi*nvico),
                as.integer(isample0),
                as.integer(nguess),
                double(nsi),
                double(nsi*(maxcomp+2)),
                siind=integer((maxcomp+2)*nvox),
                krit=double(nvox),
                as.integer(maxcomp+2),
                PACKAGE="dti")[c("siind","krit")]
  dim(z$siind) <- c(maxcomp+2,nvox)
  rbind(z$krit,z$siind)
}
pgetsii31 <- function(x,maxcomp,dgrad,th,isample1,nsi,nth,nvico,nguess,dgradi,maxc){
  nvox <- length(x)/(nsi+3)
  dim(x) <- c(nsi+3,nvox)
  z <- .Fortran("pgtsii31",
                as.double(x[1:nsi,]),
                as.double(x[nsi+1,]),#sigma2
                as.integer(nsi),
                as.integer(nvox),
                as.integer(maxcomp),
                as.double(dgrad),
                as.integer(nvico),
                as.integer(x[nsi+3,]),#iandir
                as.double(th),
                as.integer(nth),
                as.integer(x[nsi+2,]),#indth
                double(nsi*nvico),
                as.integer(isample1),
                as.integer(nguess),
                double(nsi),
                double(nsi*(maxcomp+2)),
                siind=integer((maxcomp+2)*nvox),
                krit=double(nvox),
                as.integer(maxcomp+2),
                as.double(dgradi),
                as.double(maxc),
                PACKAGE="dti")[c("siind","krit")]
  dim(z$siind) <- c(maxcomp+2,nvox)
  rbind(z$krit,z$siind)
}
pgetsiind2 <- function(x,grad,bv,nvico,dgrad,dgradi,isample,alpha,lambda,
                       maxcomp,maxc,nguess){
  #x contains  
  # si:      x[1:nsi,]
  # sigma2:  x[nsi+1,]
  nvox <- dim(x)[2]
  nsi <- dim(x)[1]-1
  siind <- .Fortran("getsii",
                    as.double(x[1:nsi,]),
                    as.double(x[nsi+1,]),
                    as.integer(nsi),
                    as.integer(nvox),
                    as.integer(maxcomp),
                    as.double(dgrad),
                    as.double(bv),
                    as.integer(nvico),
                    as.double(alpha),
                    as.double(lambda),
                    double(nsi*nvico),
                    as.integer(isample),
                    as.integer(nguess),
                    double(nsi),
                    double(nsi),#z0
                    double(nsi*(maxcomp+1)),
                    siind=integer((maxcomp+1)*nvox),
                    krit=double(nvox),
                    as.integer(maxcomp+1),
                    PACKAGE="dti")[c("siind","krit")]
  z <- matrix(0,maxcomp+2,nvox)
  z[-1,] <- siind$siind
  z[1,] <- siind$krit
  z
}
pgetsiindbv <- function(x,grad,bv,nvico,dgrad,dgradi,isample,alpha,lambda,
                       maxcomp,maxc,nguess){
  #x contains  
  # si:      x[1:nsi,]
  # sigma2:  x[nsi+1,]
  nvox <- dim(x)[2]
  nsi <- dim(x)[1]-1
  siind <- .Fortran("getsiibv",
                    as.double(x[1:nsi,]),
                    as.integer(nsi),
                    as.integer(nvox),
                    as.integer(maxcomp),
                    as.double(dgrad),
                    as.double(bv),
                    as.integer(nvico),
                    as.double(alpha),
                    as.double(lambda),
                    double(nsi*nvico),
                    as.integer(isample),
                    as.integer(nguess),
                    double(nsi),
                    double(nsi),#z0
                    double(nsi*(maxcomp+1)),
                    siind=integer((maxcomp+1)*nvox),
                    wi=double((maxcomp+1)*nvox),
                    krit=double(nvox),
                    as.integer(maxcomp+1),
                    PACKAGE="dti")[c("siind","krit","wi")]
  z <- matrix(0,2*maxcomp+3,nvox)
  z[2:(maxcomp+2),] <- siind$siind  ## vertex indices 
  z[-(1:(maxcomp+2)),] <- siind$wi  ## weights
  z[1,] <- siind$krit  ## kriterion
  z
}
