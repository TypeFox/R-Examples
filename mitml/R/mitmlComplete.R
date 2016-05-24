mitmlComplete <- function(x, print=0, force.list=FALSE){

  # .completeOne <- mitml:::.completeOne
  # .stripDataAttributes <- mitml:::.stripDataAttributes

  if(sum(print<=0)>1) stop("Only one negative or zero value is allowed in 'print'.")

  dat <- x$data
  srt <- attr(x$data,"sort")
  labs <- attr(x$data,"labels")
  method <- class(x)[2]

  m <- x$iter$m
  ind <- x$index.mat
  rpm <- x$replacement.mat

  if(class(print)%in%c("integer","numeric")){

    if(length(print)==1){
      if(print>0){
        com <- .completeOne(dat,print,ind,rpm,method)
        out <- com[srt,]
      }else{
        out <- .stripDataAttributes(dat[srt,])
      }
      if(force.list) out <- list(out)
    }else{
      out <- list()
      for(ii in print){
        if(ii>0){
          com <- .completeOne(dat,ii,ind,rpm,method)
          out <- c(out,list(com[srt,]))
        }else{
          out <- c(out,list(.stripDataAttributes(dat[srt,])))
        }
      }
    }

  }else{

    if(!print%in%c("list","all")) stop("Invalid 'print' argument.")
    out <- list()
    for(ii in 1:m){
      com <- .completeOne(dat,ii,ind,rpm,method)
      out <- c(out,list(com[srt,]))
    }

  }
  
  if(class(out)=="list") class(out) <- c("mitml.list","list")
  out

}

.completeOne <- function(x,i,ind,rpm,method){

  if(method=="jomo"){

    fac <- which(colnames(x) %in% names(attr(x,"labels")))
    nofac <- !(ind[,2] %in% fac)
    if(any(nofac)) x[ ind[nofac,] ] <- rpm[nofac,i]

    for(ff in fac){
      fi <- ind[,2]==ff
      lev <- attr(x,"labels")[[colnames(x)[ff]]]
      x[ ind[fi,] ] <- lev[rpm[fi,i]]
    }

  }else{

    x[ind] <- rpm[,i]
    
  }
  .stripDataAttributes(x)

}

.stripDataAttributes <- function(x){

  attr(x,"sort") <- NULL
  attr(x,"group") <- NULL
  attr(x,"levels") <- NULL
  attr(x,"labels") <- NULL

  x

}
