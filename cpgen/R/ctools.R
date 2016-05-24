#
# ctools.R
# Claas Heuer, June 2014
#
# Copyright (C)  2014 Claas Heuer
#
# This file is part of cpgen.
#
# cpgen is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# cpgen is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# in file.path(R.home("share"), "licenses").  If not, see
# <http://www.gnu.org/licenses/>.
#

# check_openmp

check_openmp <- function() {

if(.Call("check_openmp",PACKAGE="cpgen")) {
  
  ver = .Call("check_openmp_version",PACKAGE="cpgen")
  cat("OpenMP is available - Version: ",ver,"\n") 

} else {
      cat("OpenMP is NOT available \n")
  }

} 


# check_max_threads

get_max_threads <- function() .Call("get_max_threads",PACKAGE="cpgen")

# get_num_threads

get_num_threads <- function() options()$cpgen.threads

# check_limit_threads

#get_limit_threads <- function() .Call("get_limit_threads",PACKAGE="cpgen")

# set.threads

set_num_threads <- function(x, silent=FALSE, global=FALSE) {

if(!.Call("check_openmp",PACKAGE="cpgen")) { if(!silent) cat("OpenMP is not available, threads set to 1 \n") } else {
  if(length(x)!=1) stop("must be a scalar")
  t = as.integer(round(x,digits=0))
  if(t<=0) t=1

  options(cpgen.threads=t)
# this allows using the function to control other openmp functions in R also (MKL,...)
  if(global) .Call("set_num_threads",t,PACKAGE="cpgen")
  if(!silent) cat("Number of threads set to ",t,"\n")
  }

}

# ccov 

ccov <- function(X,lambda=0, w=NULL, compute_cor=FALSE){

  allowed=c("matrix","numeric")
  if(!class(X)%in%allowed) { stop("objects must match one of the following types: 'matrix' , 'numeric'") }
  if(is.vector(X)) X<-matrix(X)
  if(anyNA(X)) stop("No NAs allowed")
  if(missing(w)) { w<-rep(1/nrow(X),nrow(X))} 
  else { w <- w/sum(w) }
  if(lambda>1) lambda=1
  if(lambda<0) lambda=0

  .Call( "ccov", X,lambda, as.numeric(w), as.integer(compute_cor), options()$cpgen.threads, PACKAGE = "cpgen" )

}


# csolve


csolve <- function(X,y=NULL){

 allowed=c("matrix","numeric", "array", "dgCMatrix")
 a = class(X)

# this is just the dense inverse
 if(a == "matrix" & is.null(y) & get_num_threads()==1) {

   .Call( "cinverse_dense", X, PACKAGE = "cpgen" ) 

 } else {


     if(is.null(y)) y = diag(nrow(X)) 
     b = class(y)
 
     if(!a%in%allowed) stop("X must match one of the following types: 'matrix' , 'numeric', 'array', 'dgCMatrix'") 
     if(!b%in%allowed[c(1,2,3)]) stop("y must match one of the following types: 'matrix' , 'numeric', 'array'") 
     if(anyNA(X) | anyNA(y))  stop("no NAs allowed") 

     if(is.vector(X) | is.array(X)) { X = as.matrix(X); a = "matrix" } 
     if(is.vector(y) | is.array(y)) { y = as.matrix(y); b = "matrix" } 

     if(dim(X)[2]!=dim(y)[1]) {stop("ncol(X) doesn't match nrow(y)")}

     if(a == "matrix") {

       .Call( "csolve", X,y,options()$cpgen.threads ,PACKAGE = "cpgen" )

     } else {

         .Call( "csolve_sparse", X,y,options()$cpgen.threads ,PACKAGE = "cpgen" )

       }

  }

}


# cscanx

cscanx <- function(path){
	.Call( "cscanx", path ,PACKAGE = "cpgen" )
}




# cgrm.A 

cgrm.A <- function(X, lambda=0, yang=FALSE){

         if(class(X)!= "matrix") stop("X must be an object of class 'matrix'")
         if(anyNA(X)) stop("No NAs allowed in X")
         if(lambda>1) lambda=1
         if(lambda<0) lambda=0        
	 .Call( "camat", X, lambda, yang, options()$cpgen.threads ,PACKAGE = "cpgen" )
}

# cgrm.D

cgrm.D <- function(X, lambda=0){

         if(class(X)!= "matrix") stop("X must be an object of class 'matrix'")
         if(anyNA(X)) stop("No NAs allowed in X")
         if(lambda>1) lambda=1
         if(lambda<0) lambda=0
	.Call( "cdmat", X, lambda, options()$cpgen.threads ,PACKAGE = "cpgen" )
}

# cgrm

cgrm <- function(X, w = NULL, lambda=0){

         if(class(X)!= "matrix") stop("X must be an object of class 'matrix'")
         if(anyNA(X)) stop("No NAs allowed in X")
         isw = TRUE
	 if(missing(w)) {w = rep(1,ncol(X)); isw = FALSE} else {
	   if(!is.vector(w) | !is.numeric(w)) stop("weights must be passed as a numeric vector")
           if(length(w)!=ncol(X)) stop("weight vector must have as many items as columns in X") 
           vars = ccolmv(X,compute_var=T)
           var_zero = sum(vars==0)
           if(var_zero>0) {
             cat(paste(var_zero," Columns with zero variance (non-polymorphic) omitted\n",sep=""))
	     X = X[,vars>0]
	     w = w[vars>0]
           }
	 }
         if(lambda>1) lambda=1
         if(lambda<0) lambda=0

	 .Call( "cgrm", X, w,isw, lambda, options()$cpgen.threads ,PACKAGE = "cpgen" )
}


# parallel crossproduct operator

 `%c%` <- function(X,Y) {

 allowed=c("matrix","dgCMatrix","numeric","array")
 a = class(X)
 b = class(Y)

 if(sum(c(a,b)%in%allowed)!=2) { stop("objects must match one of the following types: 'matrix', 'dgCMatrix', 'numeric', 'array'") }
 if(anyNA(X) | anyNA(Y)) { stop("no NAs allowed") }

 if(is.vector(X) | is.array(X)) { X = as.matrix(X); a = "matrix" } 
 if(is.vector(Y) | is.array(Y)) { Y = as.matrix(Y); b = "matrix" } 

 if(dim(X)[2]!=dim(Y)[1]) {stop("ncol(X) doesn't match nrow(Y)")}

 if(a == "matrix"){
   if(b == "matrix") {.Call( "ccp_dense_dense", X,Y,options()$cpgen.threads ,PACKAGE = "cpgen" ) } 
     else { .Call( "ccp_dense_sparse", X,Y,PACKAGE = "cpgen" ) } 
 } else {
   if(b == "matrix") { .Call( "ccp_sparse_dense", X,Y,PACKAGE = "cpgen" ) }
     else { .Call( "ccp_sparse_sparse", X,Y,PACKAGE = "cpgen" ) } 
 } 
 

}


# cCV


cCV <- function(y,folds=5,reps=1,matrix=FALSE,seed=NULL) {

if(!is.vector(y)) stop("y must be a vector")
if(folds < 2) stop("'folds' must be larger than one")
if(reps < 1) reps = 1

reps = as.integer(reps)
folds = as.integer(folds)

isy <- (1:length(y))[!is.na(y)]

n_vec <- rep(NA,folds)
start <- rep(NA,folds)
end <- rep(NA,folds)

n_vec[1:(folds-1)] <- round(length(isy)/folds,digits=0)
n_vec[folds] <- length(isy) - sum(n_vec[1:(folds-1)])


for(j in 1:folds) {
if(j==1) { start[j] = 1 } else { start[j] = end[j-1]+1 }
end[j] = sum(n_vec[1:j]) 
}

count=0
cv_pheno <- list(folds*reps)

if(missing(seed)) { seed = as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) }
set.seed(seed)

  for(i in 1:reps) {

    ids <- sample(isy,length(isy),replace=F)
 
    for(j in 1:folds) {

      count=count+1
      y_temp = y    
      y_temp[ids[start[j]:end[j]]] <- NA 
      cv_pheno[[count]] <- y_temp  

    }

  }


if(matrix) { return(matrix(unlist(cv_pheno),nrow=length(y),ncol=length(cv_pheno),byrow=FALSE)) } else {
return(cv_pheno) }

}




# cmaf

cmaf <- function(X)	{
   maf <- .Call( "cmaf", X ,PACKAGE = "cpgen" )[1,]
   maf[maf>0.5] = 1-maf[maf>0.5]
   return(maf)
}



# ccross

ccross <- function(X,D=NULL){

         if(class(X)!= "matrix") stop("X must be an object of class 'matrix'")
         if(anyNA(X)) stop("No NAs allowed in X")
	 if(missing(D)) { D = rep(1,ncol(X)) } else {
	   if(!is.vector(D) | !is.numeric(D)) stop("D must be passed as a numeric vector") }
         if(length(D)!=ncol(X)) stop("vector D must have as many items as columns in X") 
	 .Call( "ccross", X, D, options()$cpgen.threads, PACKAGE = "cpgen" )
}
   

# cscale_inplace

cscale_inplace <- function(X,means=NULL, vars=NULL, scale=FALSE){

   if(class(X)!= "matrix") stop("X must be an object of class 'matrix'")
   if(anyNA(X)) stop("No NAs allowed in X")

   if(is.null(means))  means = ccolmv(X)
   if(!is.vector(means) | !is.numeric(means)) stop("'means' must be passed as a numeric vector")
   if(length(means)!=ncol(X)) stop("vector 'means' must have as many items as columns in X") 
   if(anyNA(means)) stop("vector 'means' has NAs") 

   if(is.null(vars)) if(scale) { vars = ccolmv(X,compute_var=T) } else { vars = rep(1,ncol(X)) } 
   if(!is.vector(vars) | !is.numeric(vars)) stop("'vars' must be passed as a numeric vector") 
   if(length(vars)!=ncol(X)) stop("vector 'vars' must have as many items as columns in X") 
   if(anyNA(vars)) stop("vector 'vars' has NAs") 
   if(any(vars==0)) stop("vector 'vars' has zeros") 

   ans <- .Call( "cscale_inplace", X, means, vars, scale, options()$cpgen.threads, PACKAGE = "cpgen" )
}
   


# Symmetrix Matrix Power operator %^% - taken from: http://stackoverflow.com/questions/16172731/how-to-compute-the-power-of-a-matrix-in-r

#`%^%` <- function(S, power) with(eigen(S), vectors %c% (values^power * t(vectors))) 
#`%^%` <- function(S, power) with(eigen(S), vectors %c% sparseMatrix(i=1:length(values),j=1:length(values),x=values^power) %c% t(vectors)) 
`%**%` <- function(X, power) with(eigen(X), ccross(vectors,values^power)) 


#####


#########################
### Internal Use only ###
#########################



# ccolmv

ccolmv <- function(X,compute_var=FALSE){

 allowed=c("matrix","dgCMatrix")
 a = class(X)

 if(!a%in%allowed) { stop("objects must match one of the following types: 'matrix', 'dgCMatrix', 'numeric'") }
 if(anyNA(X)) { stop("no NAs allowed") }
 
 if(a=="dgCMatrix") {
 
   .Call("ccolmv_sparse", X, compute_var, PACKAGE = "cpgen")
     
  } else {
 
      .Call( "ccolmv_dense", X, compute_var, PACKAGE = "cpgen" )
 
    } 

}


# ctrace

ctrace <- function(X){
 if(class(X)!="matrix"){stop("object must be of type 'matrix'")}
 .Call("ctrace",X,PACKAGE="cpgen")
}




#############################################
#### Not to be included in final package ####
#############################################


## Credit: Taken from:  http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session
# improved list of objects
#.ls.objects <- function (pos = 1, pattern, order.by,
#                        decreasing=FALSE, head=FALSE, n=5) {
#    napply <- function(names, fn) sapply(names, function(x)
#                                         fn(get(x, pos = pos)))
#    names <- ls(pos = pos, pattern = pattern)
#    obj.class <- napply(names, function(x) as.character(class(x))[1])
#    obj.mode <- napply(names, mode)
#    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
#    obj.prettysize <- napply(names, function(x) {
#                           capture.output(print(object.size(x), units = "auto")) })
#    obj.size <- napply(names, object.size)
#    obj.dim <- t(napply(names, function(x)
#                        as.numeric(dim(x))[1:2]))
#    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
#    obj.dim[vec, 1] <- napply(names, length)[vec]
#    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
#    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
#    if (!missing(order.by))
#        out <- out[order(out[[order.by]], decreasing=decreasing), ]
#    if (head)
#        out <- head(out, n)
#    out
#}
# 
## shorthand
#lsos <- function(..., n=10) {
#    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
#}


#### randmatrix

rand_data <- function(n=500,p_marker=10000,h2=0.3,prop_qtl=0.01,seed=NULL) { 

if(h2>1) h2 <- 2
if(h2<0) h2 <- 0
if(prop_qtl > 1) prop_qtl <- 1
if(prop_qtl < 0) prop_qtl <- 0
if(missing(seed)) { seed = as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) }
set.seed(seed)

### generate random data
nindividuals = n
nmarkers = p_marker

# sample allele frequencies 
p <- sample(seq(0.1,0.9,by=0.05),nmarkers,replace=T)
q <- 1-p

# random marker matrix
M<-array(NA,dim=c(nindividuals,nmarkers))

# fill matrix with marker-covariates(-1,0,1) - HWE
for (i in 1:ncol(M)) { M[,i] <- sample(c(1,0,-1),nindividuals,prob=c(p[i]**2,2*p[i]*q[i],q[i]**2),replace=T) }


# generate random phenotypes with one qtl
n_qtl <- round(p_marker*prop_qtl,digits=0)
if(n_qtl==0) n_qtl <- 1
Vp=1
h2=h2
# sample variance explained by single qtl form chisquare
#var_a = rchisq(n_qtl,Vp*h2)/n_qtl

#same variance for every marker
var_a = Vp*h2 / n_qtl
qtl<-sample(1:ncol(M),n_qtl)
alpha=sqrt((var_a)/(2*p[qtl]*q[qtl]))

y <- stats::rnorm(nindividuals,0,sqrt(Vp*(1-h2)))
y = y + (M[,qtl]%c%alpha)[,1]

# export
M <<- M
y <<- y

}





