"ipd" <- function(mat){ # "ipd" == "Is Positive Definite"
  if(nrow(mat) != ncol(mat)){
    print("not square")
    return(FALSE)
  } else if(max(abs(mat - t(Conj(mat))))>1e-10){
    print("not symmetric")
    return(FALSE)
  } else if(any(Re(eigen(mat,TRUE,TRUE)$values) < 0)){
    print("at least one negative eigenvalue")
    return(FALSE)
  } else {
    return(TRUE)
  }
}

"ss" <- function(A,B,Ainv,Binv){
  if(missing(Ainv) | missing(Binv)){
    return((det(crossprod((A+B)/2,(solve(A)+solve(B))/2)))^(-0.25))
  } else {
    return((det(crossprod((A+B)/2,(Ainv    + Binv   )/2)))^(-0.25))
  }
}

"ss_matrix_simple" <- function(hp,useM=TRUE){
  k <- length(levels(hp))
  B <- B(hp)
  M <- M(hp)
  out <- matrix(0,k,k)
  for(i in seq_len(k)){
    for(j in seq_len(k)){
      if(useM){
        out[i,j] <- ss(B[,,i],B[,,j])*M[i,j]
      } else {
        out[i,j] <- ss(B[,,i],B[,,j])
      }
    }
  }
  return(out)
}

"ss_matrix" <- function(hp,useM=TRUE){
  stopifnot(is.mhp(hp))
  M <- M(hp)
  B <- B(hp)
  
  if(useM){ # multiply by M
    f <- function(i,j){
      ss(B[,,i],B[,,j])* M[i,j]
    }
  } else { # use M==matrix(1,k,k): the highest possible correlation
    f <- function(i,j){
      ss(B[,,i],B[,,j]) 
    }
  }

  jj <- levels(hp)
  
  x <- seq_along(jj)
  y <- seq_along(jj)

  out <-
  do.call("cbind",
          sapply(x, function(x) {
            do.call("rbind",
                    sapply(y,
                           function(y) {f(x, y)},
                           simplify=FALSE
                           )
                    )
          },
                 simplify=FALSE )
          )
  rownames(out) <- jj
  colnames(out) <- jj
  return(out)
}

"default_LoF" <- function(x){
  f <- function(y){
      out <- c(1,y)
      names(out) <- c("const" , names(x))
      return(out)
    }
  
  LoF <- rep(list(f),length(levels(x)))
  names(LoF) <- levels(x)
  return(LoF)
}

"as.separate" <- function(expt){
  x <- get_mdm(expt)
  d <- get_obs(expt)
  stopifnot(is.mdm(x))
  jj <- as.list(x)
  if(missing(d)){return(jj)}  
  dd <- split(d,types(x))
    out <- list()
    for(i in seq_along(levels(types(x)))){
      out[[i]] <- list(
                       val = as.matrix(jj[[i]]),
                       obs = dd[[i]]
                       )
    }
    names(out) <- levels(x)
    return(out) 
}

"regressor"  <- function(x, LoF=NULL){

  if(is.null(LoF)){  # as opposed to a list of functions
    LoF <- default_LoF(x)
  }

  xin <- x
  x <- xold(xin)
  k <- length(levels(xin))
  stopifnot(identical(names(LoF),levels(xin)))
  
  jj <- x[1,,drop=TRUE]
  names(jj) <- NULL
  
  named_thing <- NULL
  n <- rep(0,k)
  for(i in seq_len(k)){
    func_out <- LoF[[i]](jj)
    n[i] <- length(func_out)
    named_thing <- c(named_thing,func_out)
  }
  cn <- cumsum(n)
  
  f <- function(i){
    f_out <- rep(FALSE,sum(n))
    if(i==1){
      f_out[seq(from=1,to=cn[1])] <-  TRUE


          } else {
      f_out[seq(from=1+cn[i-1],to=cn[i])] <- TRUE
    }
    return(which(f_out))  # return from f(), not regressor()
  }
     
  out <- matrix(0,nrow(x),sum(n))
  colnames(out) <- names(named_thing)
  rownames(out) <- rownames(xin)
  
  for(i in seq_along(levels(xin))){
    sub <- x[levels(xin)[i] == types(xin),,drop=FALSE]
    ind <- which(types(xin)==levels(xin)[i])
    out[ind,f(i)] <- t(apply(sub,1,LoF[[i]]))
  }
  return(out)
}

"var.matrix" <- function(x1,x2=x1,hp, ...){

  stopifnot(is.mdm(x1))
  stopifnot(is.mdm(x2))
  stopifnot(compatible(x1,x2))
  stopifnot(compatible(x1,hp))
    
  t1 <- types(x1)
  t2 <- types(x2)

  Sigma <- matrix(NA,nrow(xold(x1)),nrow(xold(x2)))
  rownames(Sigma) <- rownames(x1)
  colnames(Sigma) <- rownames(x2)

  B <- B(hp)

  for(i in levels(t1)){
    for(j in levels(t2)){
        ni <- which(i == levels(t1))
        nj <- which(j == levels(t2))
        ii <- which(i==t1)
        jj <- which(j==t2)


        if(FALSE){
          print(i)
          print(j)
          print("-------")
        }
        Sigma[ii,jj] <-
          ss(B[,,ni],B[,,nj]) * M(hp)[ni,nj]  *   
            corr.matrix(xold(x2)[jj,,drop=FALSE],xold(x1)[ii,,drop=FALSE],
                  pos.def.matrix = solve(solve(B[,,i])/2+solve(B[,,j])/2), ...)
      }
  }
  return(Sigma)
}

"beta_hat" <- function(expt,hp,LoF,...){
  if(missing(LoF)){LoF <- default_LoF(expt)}
  mm <- get_mdm(expt)
  d  <- get_obs(expt)
  betahat_mult_Sigma(H=regressor(mm,LoF), Sigma=var.matrix(x1=mm,hp=hp,...), d=d)
}

"betahat_mult"  <- function(H, Sigmainv, d){
  out <- as.vector(solve(quad.form(Sigmainv, H), crossprod(crossprod(Sigmainv, H), d)))
  names(out) <- colnames(H)
  return(out)
}

"betahat_mult_Sigma" <- function(H, Sigma, d){
  out <- as.vector(solve(quad.form.inv(Sigma, H), crossprod(H, solve(Sigma, d))))
  names(out) <- colnames(H)
  return(out)
}

"eq2.36" <- function(H, Sigmainv, d, log=TRUE){
  f <- function(m){ c(determinant(m,logarithm=TRUE)$modulus) }
    
  betahat <- betahat_mult(H, Sigmainv, d)
  out <- f(Sigmainv)-f(quad.form(Sigmainv, H))-quad.form(Sigmainv,d - H%*%betahat)
  out <- out/2
  if(log){
    return(out)
  } else {
    return(exp(out))
  }
}

"eq2.36_Sigma" <- function(H, Sigma, d){
betahat <- betahat_mult_Sigma(H,Sigma,d)
sqrt((1/det(Sigma)) / det(quad.form.inv(Sigma,H))) * 
  exp( -0.5*quad.form.inv(Sigma,   d-H %*% betahat))
}

"compatible"  <- function(x1, x2){
  if(
     identical( names(x1), names(x2)) & 
     identical(levels(x1),levels(x2))
     ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

"cstar" <- function(x1, x2=x1, expt, hp, LoF=NULL, Sigmainv=NULL, ...){ # x -> x1;  xdash -> x2

  x <- get_mdm(expt)
  d <- get_obs(expt)

  stopifnot(is.mdm(x1))
  stopifnot(is.mdm(x2))
  stopifnot(is.mdm(x ))
  
  stopifnot(compatible(x1,x2))
  stopifnot(compatible(x1,x))
  stopifnot(compatible(x1,hp))

  if(is.null(Sigmainv)){
    Sigmainv <- solve(var.matrix(x,x,hp,...))
  }
  
  TEEx1 <- var.matrix(x1=x,x2=x1,hp,...)
  TEEx2 <- var.matrix(x1=x,x2=x2,hp,...)

  H   <- regressor(x ,LoF=LoF)
  Hx1 <- regressor(x1,LoF=LoF)
  Hx2 <- regressor(x2,LoF=LoF)
  
  bit1 <- var.matrix(x1=x1,x2=x2,hp,...)
  bit2 <- quad.3form(Sigmainv, TEEx1, TEEx2)

  jjleft  <- Hx1 - quad.3form(Sigmainv, TEEx1, H)
  jjright <- Hx2 - quad.3form(Sigmainv, TEEx2, H)

  bit3 <- jjleft %*% tcrossprod(solve(quad.form(Sigmainv, H)), jjright)
  
  return(bit1-bit2+bit3)
}


"multem" <- function(x, expt, hp=NULL, LoF=NULL, give=FALSE, Sigmainv=NULL, ...){  # 'multem' == 'MULTivariate EMulator'

  x_known <- get_mdm(expt)
  d       <- get_obs(expt)
  
  stopifnot(is.mdm(x_known))
  stopifnot(is.mdm(x))
  stopifnot(compatible(x_known,hp))

  
  if(is.null(Sigmainv)){
    Sigmainv <- solve(var.matrix(x_known,x_known,hp, ...))
  }

  H <- regressor(x_known, LoF)
  bhat <- betahat_mult(H, Sigmainv, d)

  Hx <- regressor(x, LoF)
  TEE <- var.matrix(x1=x_known,x2=x,hp=hp, ...)

  res <- d-regressor(x_known, LoF)%*%bhat  # "res" == residual

  out <- drop(Hx %*%  bhat +  crossprod(TEE, Sigmainv) %*% res)
  names(out) <- rownames(x)
  
  if(give){
    jj <-
      cstar(x1=x, x2=x, expt=expt, hp=hp, LoF=LoF, Sigmainv=NULL)
    rownames(jj) <- rownames(x)
    colnames(jj) <- rownames(x)
    return(list(mstar=out, cstar=jj))
  } else {
    return(out)
  }
} 

"obs_maker" <- function(x, hp, LoF, beta, Sigma=NULL, ...){
  if(is.null(Sigma)){Sigma <- var.matrix(x1=x, hp=hp, ...)}
  d <- rmvnorm(n=1, mean=regressor(x, LoF)%*% beta, sigma=Sigma)
  d <- as.vector(d)
  names(d) <- rownames(x)
  return(d)
}

"toy_mm_maker" <- function(na,nb,nc,include_first=TRUE){
  
  out <- latin.hypercube(na+nb+nc,4)
  if(include_first){
    out <- rbind(matrix(0.5, 3,4),out)
  }
  
  colnames(out) <-letters[1:4]
  
  tt <- c("temp","rain","humidity")
  tt1 <- c("t","r","h")

  jjind <- c(rep(1,na),rep(2,nb),rep(3,nc))
  
  if(include_first){
    jjtypes <- factor(tt[c(1:3,jjind)], levels=tt)
    jj1 <-           tt1[c(1:3,jjind)]
  } else {
    jjtypes <- factor(tt[jjind], levels=tt)
    jj1 <-           tt1[jjind]
  }

  jjnum <- rep(0,length(jjtypes))

  w <- which(jjtypes == tt[1])
  jjnum[w] <- seq_along(w)
  
  w <- which(jjtypes == tt[2])
  jjnum[w] <- seq_along(w)
  
  w <- which(jjtypes == tt[3])
  jjnum[w] <- seq_along(w)
  rownames(out) <- paste(jj1,jjnum,sep='')
  mdm(out,jjtypes)
}

"optimal_B" <- function(expt, LoF, start_hp, option='a', verbose=FALSE, ...){ # returns a B

  mm <- get_mdm(expt)
  
  if(missing(start_hp)){start_hp <- as.mhp(mm)}
  stopifnot(compatible(mm,start_hp))
  seps <- as.separate(expt)
  
  B <- B(start_hp)
  
  if(option == 'a'){    ## each B[,,i] a multiple of the Identity matrix (Id)
    for(i in seq_along(names(seps))){ #care! names(seps) =
                                      #c("temp","rain","humidity")
      if(verbose){print(paste("calculating ",i," / ",length(names(seps)),sep=""))}
      jj <- optimal.scale(seps[[i]]$val, seps[[i]]$obs, func=LoF[[i]], ...)
      B[,,i] <- diag(jj,nrow=dim(B)[1])
      if(verbose){print(paste("calculated ",i," / ",length(names(seps)),sep=""))}
    }
  } else if(option == 'b'){   ## each B[,,i] a diagonal matrix
    for(i in seq_along(names(seps))){  #care! names(seps) = c("temp","rain","humidity")
      jj <- optimal.scales(val          = seps[[i]]$val,
                           scales.start = diag(B[,,i]),
                           d            = seps[[i]]$obs,
                           func         = LoF[[i]], ...
                           )
      B[,,i] <- diag(jj,nrow=dim(B)[1])
    }
  } else if(option =='c'){  ## each B[,,i] diagonal; all the B[,,i] identical
    B <- optimal_identical_B(expt, LoF, start_hp, verbose=FALSE, ...)
  } else {
    stop("option must be 'a' or 'b' or 'c'")
  }
  return(B)
}

"optimal_identical_B" <- function(expt, LoF, start_hp, verbose=FALSE, ...){
  seps <- as.separate(expt)
  p <- length(seps)
  objective <- function(diag_of_B){# objective function for optimizer
    out <- 0
    for(i in seq_len(p)){
      jj <- seps[[i]]
      # following line will not work with emulator_1.1-8
      out <- out + scales.likelihood(scales=exp(diag_of_B), xold=jj$val,d=jj$obs, func=LoF[[i]],give_log=TRUE)
    }
    return(-out)  #return negative because optim() minimizes
  }

  start_B <- log(diag(apply(B(start_hp),1:2,mean)))  #ie a vector
  out <- optim(par=start_B,fn=objective,...)
  if(verbose){
    return(out)
  } else {
    return(kronecker(diag(exp(out$par)),array(1,c(1,1,p))))
  }
}

"optimal_diag_M" <- function(expt, LoF, start_hp){
                                        # This means the diagonal of
                                        # the *M* matrix: the marginal
                                        # variances (NB: analysis
                                        # conditional upon B)

  jj <- as.separate(expt)
  B_cond <- B(start_hp)
  
  shs <- rep(NA,length(levels(expt))) # 'SHS' == 'SigmaHatSquared'
  for(i in seq_along(shs)){
    val <-  jj[[i]]$val
    d1   <- jj[[i]]$obs
    shs[i] <-
    sigmahatsquared(
                    H    = t(apply(val, 1, LoF[[i]])),
                    Ainv = solve(corr.matrix(xold=val, scales=diag(B_cond[,,i]))),
                    d    = d1
                    )

  }
  return(shs) # ie a vector
}
 
"optimal_M" <- function(expt, LoF, start_hp, ...){
                                        #analysis conditional on
                                        #B(start_hp) [determined from
                                        #optimal_B()]; and the
                                        #diagonal elements of M
                                        #[determined from
                                        #optimal_diag_M()]


  mm <- get_mdm(expt)
  d  <- get_obs(expt)

  make_M <- function(vec){  #returns an 'M' matrix from a vector.
    jj_M <- M(start_hp)
    jj_M[upper.tri(jj_M,diag=TRUE )] <- vec
    jj_M[lower.tri(jj_M,diag=FALSE)] <- 0
    jj_M <- jj_M + t(jj_M)
    diag(jj_M) <- diag(jj_M)/2
    return(jj_M)
  }

  jj_hp <- start_hp
  f <- function(vec){  #objective function to be minimized, ie -log(likelihood)
    jj <- make_M(vec)
    if(any(Re(eigen(jj,TRUE,TRUE)$values) < 0)){
      return(Inf)
    } else {
      M(jj_hp) <- make_M(vec)
      out <- 
      eq2.36(
             H        = regressor(mm, LoF),
             Sigmainv = solve(var.matrix(x1=mm, hp=jj_hp, ...)),
             d        = d,
             log      = TRUE)

#      jjM <- M(jj_hp)
#      jjM <- jjM / sqrt(outer(diag(jjM),diag(jjM)))
      
      return(-out)
    }
  }

  jj <- M(start_hp)
  start_vec <- jj[upper.tri(jj,diag=TRUE)]
#  opt_vec <- optim(par=start_vec,fn=f, control=list(maxit=10000,trace=9))

  opt_vec <- optim(par=start_vec,fn=f, ...)
  
  return(make_M(opt_vec$par))
}
    
"optimal_params" <- function(expt, LoF, start_hp, option='a', ...){

  mm <- get_mdm(expt)
  d  <- get_obs(expt)

  if(missing(start_hp)){
    out <- as.mhp(mm)
  } else {
    out <- start_hp
  }

  if(missing(LoF)){
    LoF <- default_LoF(mm)
  }
  
  B(out)  <- optimal_B(expt, LoF, start_hp=out, option=option, ...)  # step 1

  {
    jj <- optimal_diag_M(expt, LoF, start_hp=out)  
    jjM <- diag(M(out))
    covs <- M(out) / sqrt(kronecker(jjM,t(jjM)))
    M(out) <- sqrt(kronecker(jj,t(jj)))*covs                          # step 2
  }
  
  M(out)  <- optimal_M(expt, LoF, start_hp=out, ...)                 # step 3
  return(out)
}

"apart" <- function(X, dependent,use_rownames=TRUE){
  if(is.logical(dependent)){dependent <- which(dependent)}
  xold <- X[,-dependent,drop=FALSE]
  jj <- colnames(X)[dependent]
  x <- do.call("rbind", rep(list(xold),length(jj)))
  types <- as.factor(rep(jj,each=nrow(xold)))

  if(use_rownames){
    rownames(x) <- paste(rep(rownames(xold),length(jj)),
                         rep(jj,each=nrow(xold)),
                         sep="_")
  } else {
    rownames(x) <- NULL
  }
  return(experiment(
              mm  = mdm(x,types),
              obs = as.vector(X[,dependent])
              ))
}

"showmap" <- function(z, pc,  ...){
  long <- seq(from=2.81,to=357,length.out=64)
  lat  <- c(-79.811531,seq(from=-74.81,to=86,len=30),86.6)
  z <- t(matrix(z,32,64))
  
  
  if(pc){  #ie Principal Component
    cp <- terrain.colors
    text <- " "
    z <-  z/sd(c(z))
  } else {  #a regular map of temperature
    cp <- heat.colors
    text <- "temp (C)"
  }

  filled.contour(x = long,
                 y = lat,
                   z = z,
                 color.palette = cp,
                 axes = TRUE,
                 key.title=title(main=text),
                 xlab = 'Longitude', ylab = 'Latitude',
                 plot.axes={axis(1);axis(2);
                              contour(long,lat, landmask, level = 0.5, drawlabels = FALSE, method = 'simple',add = TRUE, lwd = 1.2, col = 'black')},
                   ...)
  return(invisible(z))
}
