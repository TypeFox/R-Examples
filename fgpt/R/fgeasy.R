fgeasy <-
function(xy, group=1, marks, iter=999, ratio=1, scale.seq=seq(from=0, to=max(dist(xy)), length.out=21)[2:21], bootstrap=FALSE, pairwise=FALSE, correlate=FALSE){  
  
  # check data input
  
  if(dim(xy)[2]!=2) {stop("xy does not have the right dimensions")}
  if(!is.numeric(xy)) {stop("xy must be numeric")}  
  if(!is.numeric(marks)) {stop("marks must be numeric")}
  if(is.vector(marks)){if (group!=1) {if (!(dim(xy)[1]==length(group) && dim(xy)[1]==length(marks) && length(group)==length(marks))){stop("one or more of the following arrays does not have the appropriate length or dimensions: xy, group, marks")}}}
  if(is.vector(marks)){if (group==1) {if (dim(xy)[1]!=length(marks)){stop("one or more of the following arrays does not have the appropriate length or dimensions: xy, marks")}}}
  if(!is.vector(marks)){if (group!=1) {if (!(dim(xy)[1]==length(group) && dim(xy)[1]==dim(marks)[1] && length(group)==dim(marks)[1])){stop("one or more of the following arrays does not have the appropriate length or dimensions: xy, group, marks")}}}
  if(!is.vector(marks)){if (group==1) {if (dim(xy)[1]!=dim(marks)[1]){stop("one or more of the following arrays does not have the appropriate length or dimensions: xy, marks")}}}
  if(!is.vector(marks)){if (dim(marks)[1]!=dim(marks)[2]){if (dim(marks)[2]!=2){stop("marks needs to be a vector or a matrix with either 2 columns or equal numbers of rows and columns")}}}
  if(!is.numeric(iter)) {stop("iter needs to be numeric")}  
  if(!is.numeric(ratio)) {stop("ratio needs to be numeric")}  
  if(!is.logical(bootstrap)) {stop("bootstrap needs to be logical")} 
  if(!is.logical(pairwise)) {stop("pairwise needs to be logical")}  
  if(length(iter)!=1) {stop("iter is not allowed to have more than one value")}  
  if(length(ratio)!=1) {stop("ratio is not allowed to have more than one value")}  
  if(length(bootstrap)!=1) {stop("bootstrap is not allowed to have more than one value")}   
  if(length(pairwise)!=1) {stop("pairwise is not allowed to have more than one value")}  
  if(pairwise==FALSE & correlate!=FALSE) {stop("correlate needs to be FALSE, when pairwise is FALSE")}
  
  time1 <- Sys.time()
  
  xy <- as.matrix(xy)
  group <- as.integer(as.factor(group))
  
  
  if(pairwise==FALSE){
    perms <- c.list <- list()
    dist.mat <- as.matrix(dist(xy,upper = TRUE))
    diag(dist.mat) <- NA
    w.mat <- ifelse(dist.mat<=max(apply(dist.mat,1,min, na.rm=TRUE)),1,0)
    diag(w.mat) <- 0
    
    mori <- function(x,w.mat){
      x <- x-mean(x)
      length(x)*sum((x %*% t(x))*w.mat,na.rm=TRUE)/(sum(w.mat,na.rm=TRUE)*sum(x^2))
      }
    }else if(correlate==FALSE){
      perms <- m.list <- v.list <- list()
    }else{
      if(correlate==TRUE){correlate="pearson"}
      perms <- c.list <- list()
      cor.est <- function(x,y, ...){cor.test(x,y, ...,method=correlate)$estimate}
    }
  
  cat("\n")
  cat("=====================================","\n")
  cat(" Floating Grid Permutation Technique","\n")
  cat("=====================================","\n")
  cat("\n")
  cat("Progress bar for spatially restricted permutations","\n")
  cat("0% |------------------------------------------------| 100%","\n")
  cat("   ")
 
  for(i in 1:length(scale.seq)){
    
    scale <- scale.seq[i]
    
    if(bootstrap==FALSE){
      perms[[i]] <- fgperm(xy=xy,z=1:dim(xy)[1], scale=scale, group=group, iter=iter, FUN=fyshuffle, add.obs=TRUE)
    }else{
      perms[[i]] <- fgperm(xy=xy,z=1:dim(xy)[1], scale=scale, group=group, iter=iter, FUN=function(x){x[sample.int(length(x),replace=TRUE)]}, add.obs=TRUE)      
    }

    if(pairwise==FALSE){
      c.list[[i]] <- fgstat(perms[[i]],marks, FUN=mori, w.mat=w.mat)
    }else if(correlate==FALSE){
      m.list[[i]] <- fgstat(perms[[i]],marks,mean, na.rm=TRUE)
      v.list[[i]] <- fgstat(perms[[i]],marks,var, na.rm=TRUE)
    }else{
      c.list[[i]] <- fgstat(perms[[i]],marks,cor.est, na.rm=TRUE)
      }

    if(length(scale.seq)<50){
     cat(rep("=",diff(round(seq(from=0,to=50,length.out=length(scale.seq)+1)))[i]), sep="")
    }else{
      if(i %in% ceiling(1:50*(length(scale.seq)/50))){
        cat("=")
        }  
      }
    } # end of iterations
 
  cat("\n")
  cat("\n")
  cat("Non-spatial permutation test is now running.......","\n")
  
  perms[[length(scale.seq)+1]] <- c(observed=list(1:dim(xy)[1]),replicate(iter,list(fyshuffle(1:dim(xy)[1]))))
  
  if(pairwise==FALSE){
    c.list[[length(scale.seq)+1]] <- fgstat(perms[[i]],marks,FUN=mori, w.mat=w.mat)
  }else if(correlate==FALSE){
    m.list[[length(scale.seq)+1]] <- fgstat(perms[[i]],marks,mean, na.rm=TRUE)
    v.list[[length(scale.seq)+1]] <- fgstat(perms[[i]],marks,var, na.rm=TRUE)
  }else{
    c.list[[length(scale.seq)+1]] <- fgstat(perms[[i]],marks,cor.est, na.rm=TRUE)
    }
 
  cat("\n")
  time2 <- round(Sys.time()-time1,2)
  cat("\n")
  cat(paste("The analysis took", time2,attr(time2,"unit")),"\n")
  cat("\n")
  
  if(pairwise==TRUE & correlate==FALSE){
    output <- list(m.list=m.list, v.list=v.list, iter=iter, scales=c(scale.seq,Inf),correlate=correlate, pairwise=pairwise)  
  }else{
    output <- list(c.list=c.list, iter=iter, scales=c(scale.seq,Inf),correlate=correlate, pairwise=pairwise)    
    }
  
  class(output)<-"fg"
  
  return(output)
  
  }
