# Rcsdp.R
# Interface to the CSDP semidefinite programming library by Brian Borchers
# https://projects.coin-or.org/Csdp/
#
# Created: 22 February 2008
# Author: Hector Corrada Bravo

csdp <- function(C,A,b,K,control=csdp.control()) {
  prob.info <- get.prob.info(K,length(b))
  validate.data(C,A,b,prob.info)
  prob.data <- prepare.data(C,A,b,prob.info)
  write.control.file(control)
  
  ret <- .Call("csdp",
               as.integer(sum(prob.info$block.sizes)),
               as.integer(prob.info$nconstraints),
               as.integer(prob.info$nblocks),
               as.integer(c(0,prob.info$block.types)),
               as.integer(c(0,prob.info$block.sizes)),
               prob.data$C,
               prob.data$A,
               prob.data$b,
               PACKAGE="Rcsdp")

  unlink("param.csdp")
  ret[1:3] <- get.solution(ret[[1]],ret[[2]],ret[[3]],prob.info)
  structure(ret,names=c("X","Z","y","pobj","dobj","status"))
}

get.solution <- function(X,Z,y,prob.info)
  {
    list(X=blkmatrix_csdp2R(X,prob.info),
         Z=blkmatrix_csdp2R(Z,prob.info),
         y=vector_csdp2R(y))
      }

prepare.data <- function(C,A,b,prob.info)
  {
    list(C=blkmatrix_R2csdp(C,prob.info),
         A=constraints_R2csdp(A,prob.info),
         b=as.double(vector_R2csdp(b)))
  }

get.prob.info <- function(K,m) {
 if (!all.equal(names(K),c("type","size")))
   stop("Invalid cone specification 'K': elements must be 'type' and 'size'")

 if (!all(K$type %in% c("s","l")))
   stop("Invalid cone specification 'K': types must be 's' or 'l'")

 if (length(K$type) != length(K$size))
   stop("Invalid conse specification 'K': type and size elements must be of the same length")

 
 block.types <- ifelse(K$type == "s",1,2)
 nblocks <- length(K$type)
 
 block.sizes <- K$size
 nconstraints <- m;

 ret <- list(nblocks=nblocks,
             nconstraints=nconstraints,
             block.types=block.types,
             block.sizes=block.sizes)
 return(ret)
}

validate.data <- function(C,A,b,prob.info) {
  nblocks <- prob.info$nblocks
  nconstraints <- prob.info$nconstraints
  block.types <- prob.info$block.types
  block.sizes <- prob.info$block.sizes
  
  # Validate number of blocks in C
  if (length(C) != nblocks)
    stop("Number of blocks in C disagrees with K")

  # Validate number of constraint matrices in A
  if (length(A) != nconstraints)
    stop("Number of constraint matrices in A disagrees with b")

  # Validate each block of C
  for (j in 1:nblocks) {
    if ((block.types[j] == 1) && (block.sizes[j] != nrow(C[[j]]) || block.sizes[j] != ncol(C[[j]])))
      stop("Size of block ",j," in C does not agree with K")
    if ((block.types[j] == 2) && (block.sizes[j] != length(C[[j]])))
      stop("Size of block ",j," in C does not agree with K")
  }

  # Validate constraint matrices
  for (i in 1:nconstraints) {
    Ai <- A[[i]];
    if (length(Ai) != nblocks)
        stop("Number of blocks in constraint matrix ",i," in A disagrees with K")
             
    for (j in 1:nblocks) {
      if ((block.types[j] == 1) && (block.sizes[j] != nrow(Ai[[j]]) || block.sizes[j] != ncol(Ai[[j]])))
        stop("Size of block ",j," in constraint matrix ",i," in A does not agree with K")
      if ((block.types[j] == 2) && (block.sizes[j] != length(Ai[[j]])))
        stop("Size of block ",j," in constraint ",i," of A does not agree with K")
    }
  }
}

csdp.control <- function(axtol=1e-8,
                          atytol=1e-8,
                          objtol=1e-8,
                          pinftol=1e8,
                          dinftol=1e8,
                          maxiter=100,
                          minstepfrac=0.90,
                          maxstepfrac=0.97,
                          minstepp=1e-8,
                          minstepd=1e-8,
                          usexzgap=1,
                          tweakgap=0,
                          affine=0,
                          printlevel=1,
                          perturbobj=1,
                          fastmode=0)
  {
    list(axtol=axtol,
         atytol=atytol,
         objtol=objtol,
         pinftol=pinftol,
         dinftol=dinftol,
         maxiter=maxiter,
         minstepfrac=minstepfrac,
         maxstepfrac=maxstepfrac,
         minstepp=minstepp,
         minstepd=minstepd,
         usexzgap=usexzgap,
         tweakgap=tweakgap,
         affine=affine,
         printlevel=printlevel,
         perturbobj=perturbobj,
         fastmode=fastmode)
  }

write.control.file <- function(control)
  {
    fileptr <- file("param.csdp","w")
    for (i in 1:length(control))
      cat(names(control)[i],"=",control[[i]],"\n",sep="",file=fileptr)
    close(fileptr)
  }
