#
# clmm.R
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

# clmm

clmm <- function(y, X = NULL , Z = NULL, ginverse = NULL, par_random = NULL, niter=10000, burnin=5000,scale_e=0,df_e=-2, beta_posterior = FALSE, verbose = TRUE, timings = FALSE, seed = NULL, use_BLAS=FALSE){

default_scale = 0
default_df = -2
h2 = 0.3
default_GWAS_window_size_proportion = 0.01
default_GWAS_threshold = 0.01

# renamed 'random' to Z in the R function
random = Z

allowed=c("numeric", "list")
a = class(y)

if(!a %in% allowed) stop("phenotypes must match one of the following types: 'numeric' 'list'") 

if(class(y) == "list") {
  p = length(y)
  if (sum(unlist(lapply(y,is.vector))) != p) stop("phenotypes must be supplied as vectors")
  n = unlist(lapply(y,length))
  if (sum(n[1]!=n) > 0) stop("phenoytpe vectors must have same length")
  n = n[1]
  if(is.null(names(y))) names(y) <- paste("Phenotype_",1:p,sep="") 
  } else {
    if(!is.vector(y)) stop("phenotype must be supplied as vector")
    n <- length(y) 
    y <- list(y)
    names(y) <- "Phenotype_1"
  }
  
y <- lapply(y,as.numeric)  
if(any(unlist(lapply(y,function(x)var(x,na.rm=TRUE)))==0)) stop("one or more phenotypes with 0 variance detected")

if(is.null(X)) {
  X = array(1,dim=c(n,1))
  par_fixed <- list(scale=default_scale,df=default_df,sparse_or_dense="dense", name="fixed_effects", method="fixed")
  } else {                   
    if(X_is_ok(X,n,"fixed")) {
      if(class(X) == "matrix") { type = "dense" } else { type = "sparse" }
      par_fixed <- list(scale=default_scale,df=default_df,sparse_or_dense=type,name="fixed_effects",method="fixed")    
    }
  }

par_fixed$GWAS=FALSE    
par_fixed$GWAS_threshold = 0.01 
par_fixed$GWAS_window_size = 1
# added 09/2015 - ginverse
par_fixed$sparse_or_dense_ginverse = "dense"

par_random_all = list()
par_temp = par_random

if(is.null(random)) {
  random = list()
  par_random_all = list()
  for(k in 1:length(y)) { par_random_all[[k]] = list(list()) }

  } else {

  if(is.null(names(random))) names(random) = paste("Effect_",1:length(random),sep="")

  for(k in 1:length(y)) {

    par_random=par_temp

    if(is.null(par_random)) {

      par_random<-list(length(random))

      for(i in 1:length(random)){

        if(X_is_ok(random[[i]],n,names(random)[i])) {
        method = "ridge"
        if(class(random[[i]]) == "matrix") type = "dense"
        if(class(random[[i]]) == "dgCMatrix") type = "sparse"
        par_random[[i]] = list(scale=default_scale,df=default_df,sparse_or_dense=type,method=method, name=as.character(names(random)[i]), GWAS=FALSE, GWAS_threshold = 0.01, GWAS_window_size = 1) }
      }

      } else {

          if(length(par_random) != length(random)) stop(" 'par_effects' must have as many items as 'random' ")

            for(i in 1:length(par_random)) {

              if(!is.list(par_random[[i]])) par_random[[i]] <- list()
              X_is_ok(random[[i]],n,names(random)[i])
              allowed_methods = c("fixed","ridge","BayesA")
              if(is.null(par_random[[i]]$method)) par_random[[i]]$method <- "ridge"
              if(is.null(par_random[[i]]$name)) par_random[[i]]$name = as.character(names(random)[i])
              if(!par_random[[i]]$method %in% allowed_methods) stop(paste("Method must be one of: ",paste(allowed_methods,collapse=" , "),sep=""))

              if(is.null(par_random[[i]]$df[k]) | !is.numeric(par_random[[i]]$df[k]) | length(par_random[[i]]$df[k]) > 1)  {

                if(par_random[[i]]$method == "BayesA") { 

                  par_random[[i]]$df = 4.0 

                } else { 

                    par_random[[i]]$df = default_df 

                 }  
        
              } else {

                  par_random[[i]]$df = par_random[[i]]$df[k]

                } 

              if(is.null(par_random[[i]]$scale[k]) | !is.numeric(par_random[[i]]$scale[k]) | length(par_random[[i]]$scale[k]) > 1) {

                if(par_random[[i]]$method == "BayesA") { 

## Fernando et al. 2012
                  dfA = par_random[[i]]$df
                  meanVar = mean(ccolmv(random[[i]],compute_var=T))
                  varG = h2*var(y[[k]],na.rm=TRUE)
                  varMarker = varG / ncol(random[[i]]) * meanVar
                  par_random[[i]]$scale = varMarker * (dfA -2) / dfA 

                } else {

                    par_random[[i]]$scale = default_scale 

                 } 
     
              } else {

                  par_random[[i]]$scale = par_random[[i]]$scale[k]

                } 

                if(class(random[[i]]) == "matrix") { type = "dense" } else { type = "sparse" }
                par_random[[i]]$sparse_or_dense = type
# GWAS
                if(is.null(par_random[[i]]$GWAS)) {

                  par_random[[i]]$GWAS=FALSE    
                  par_random[[i]]$GWAS_threshold = 0.01 
                  par_random[[i]]$GWAS_window_size = 1

                } else {

                    if(is.null(par_random[[i]]$GWAS$threshold)) { 

                      par_random[[i]]$GWAS_threshold = default_GWAS_threshold
      
                    } else { 
 
                      par_random[[i]]$GWAS_threshold = as.numeric(par_random[[i]]$GWAS$threshold) 

                    }

                  if(is.null(par_random[[i]]$GWAS$window_size)) { 

                    par_random[[i]]$GWAS_window_size = as.integer(ncol(random[[i]]) *  default_GWAS_window_size_proportion)
      
                  } else { 
 
                      par_random[[i]]$GWAS_window_size = as.integer(par_random[[i]]$GWAS$window_size) 

                    }

                  par_random[[i]]$GWAS = TRUE
 
                  }
   
            }

        }

    par_random_all[[k]] = par_random

  }

}
  



################################
### added 09/2015 - Ginverse ###
################################

# first add the parameter "sparse_or_dense_ginverse" to the par_random list
for(i in 1:length(par_random_all)) {

  for(j in 1:length(par_random_all[[i]])) par_random_all[[i]][[j]]$sparse_or_dense_ginverse = "sparse"

}


if(!is.null(ginverse)) {

  if(length(ginverse)!= length(random)) stop("If provided, ginverse must have as many items as random. Put 'NULL' for no ginverse in the list for a particular random effect")

# check dimensions - all that matters is the number of columns
  for(i in 1:length(ginverse)) {

    if(!is.null(ginverse[[i]])) {

      if(ncol(ginverse[[i]]) != ncol(random[[i]])) stop(paste("Number of columns in design matrix: '",
                                                        par_random[[i]]$name,"' dont match dimnsion of corresponding ginverse", sep=""))

      if(!class(ginverse[[i]]) %in% c("matrix","dgCMatrix")) stop(paste("Ginverse: '",par_random[[i]]$name,
                                                                       "' must be of type 'matrix' or 'dgCMatrix'",sep=""))

# set the method for that effect - Only ridge regression allowed
# also set the type of matrix for ginverse

# this is crap - but good enough for now
      for(j in 1:length(par_random_all)) {

          par_random_all[[j]][[i]]$method = "ridge_ginverse"
          par_random_all[[j]][[i]]$sparse_or_dense_ginverse = ifelse(class(ginverse[[i]])=="matrix", "dense", "sparse")

      }

    }

  }

# if no ginverse in the list, create dummy
} else {

    ginverse <- vector("list",length(random))

  }


# RNG Seed based on system time and process id
# Taken from: http://stackoverflow.com/questions/8810338/same-random-numbers-every-time
if(is.null(seed)) { seed = as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) }

par_mcmc = list()
verbose_single = verbose
if(timings | length(y) > 1) verbose_single = FALSE

if(length(y)>1) {

  if(length(scale_e)!=length(y)) scale_e = rep(scale_e[1], length(y))
  if(length(df_e)!=length(y)) df_e = rep(df_e[1], length(y))
  
}

# for CV timings is not a good thing
if(length(y) > 1) timings = FALSE

for(i in 1:length(y)) {

  par_mcmc[[i]] = list(niter=niter, burnin=burnin, full_output=beta_posterior, verbose=verbose_single,
  timings = timings, scale_e = scale_e[i], df_e = df_e[i], seed = as.character(seed), name=as.character(names(y)[i]))

}

 mod <- .Call("clmm",y, X , par_fixed ,random, par_random_all ,par_mcmc, verbose=verbose, options()$cpgen.threads, use_BLAS, ginverse, PACKAGE = "cpgen" )

 if(length(y) == 1) { return(mod[[1]]) } else { return(mod) }

#return(list(y, X , par_fixed ,random, par_random_all ,par_mcmc, verbose=verbose, options()$cpgen.threads, use_BLAS, ginverse))

}






get_pred <- function(mod) {

return(matrix(unlist(lapply(mod,function(x)x$Predicted)),ncol=length(mod),nrow=length(mod[[1]]$Predicted)))

}


get_cor <- function(predictions,cv_pheno,y) {

cv_vec <- matrix(unlist(cv_pheno),nrow=length(y),ncol=length(cv_pheno),byrow=FALSE)

mean_pred <- rep(NA,nrow(predictions))

for(i in 1:nrow(predictions)) {

mean_pred[i] <- mean(predictions[i,which(is.na(cv_vec[i,]))])

}

return(cor(mean_pred,y,use="pairwise.complete.obs"))

}

# cGBLUP


cGBLUP <- function(y,G,X=NULL, scale_a = 0, df_a = -2, scale_e = 0, df_e = -2,niter = 10000, burnin = 5000, seed = NULL, verbose=TRUE){

isy <- (1:length(y))[!is.na(y)]
if(length(y) != nrow(G)) stop("dimension of y and G dont match")

if(verbose) cat("\nComputing Eigen Decomposition\n")

if(length(isy) < length(y)) {
  UD <- eigen(G[isy,isy]) } else {
    UD <- eigen(G) }

n <- length(isy)

if(is.null(X)) X = rep(1,length(y[isy]))

Uy <- (t(UD$vectors)%c%y[isy])[,1]
UX <- t(UD$vectors)%c%X

D_sqrt <- sqrt(UD$values)
Z<-sparseMatrix(i=1:n,j=1:n,x=D_sqrt)

par_random <- list(list(scale=scale_a,df=df_a,sparse_or_dense="sparse",method="ridge"))

if(verbose) cat("Running Model\n")
if(is.null(seed)) { seed = as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) }

# set the number of threads to 1 for clmm
old_threads <- get_num_threads()
set_num_threads(1,silent=TRUE)

mod <- clmm(y = Uy, X = UX , Z = list(Z), par_random = par_random, scale_e=scale_e, df_e=df_e, verbose=verbose, niter=niter, burnin=burnin, seed=seed)

# set number of threads to old value
set_num_threads(old_threads,silent=TRUE)

u <- rep(NA,length(y))


u[isy] <- UD$vectors %c% (D_sqrt * mod[[4]]$posterior$estimates_mean)
if(length(isy) < length(y)) { u[-isy] <- G[-isy,isy] %c% csolve(G[isy,isy],u[isy]) }

e<-mod$Residual_Variance$Posterior

return(list(var_e = mod$Residual_Variance$Posterior_Mean,
	    var_a = mod[[4]]$posterior$variance_mean, 
            b = mod[[3]]$posterior$estimates_mean,
	    a = u,
	    posterior_var_e = mod$Residual_Variance$Posterior,
	    posterior_var_a = mod[[4]]$posterior$variance))


}






X_is_ok <- function(X,n,name) {

allowed=c("matrix","dgCMatrix")
a = class(X)
#if(sum(a%in%allowed)!=1) stop(paste(c("lol","rofl"))) 
if(sum(a%in%allowed)!=1) stop(paste("design matrix '",name,"' must match one of the following types: ",paste(allowed,collapse=" , "),sep="")) 

if(anyNA(X)) { stop(paste("No NAs allowed in design matrix '", name,"'", sep="")) } 
      
if(a=="matrix" | a=="dgCMatrix") { if(nrow(X) != n) stop(paste("Number of rows in design matrix '",name,"' doesnt match number of observations in y",sep="")) }


return(1) 

}



### GWAS

#cGWAS.BR <- function(mod, M, window_size, threshold, sliding_window=FALSE, verbose=TRUE) {		
#		
#  niter = mod$mcmc$niter		
#  burnin = mod$mcmc$burnin
#		
#  n_windows = ifelse(sliding_window, as.integer(ncol(M) - window_size + 1), as.integer(ncol(M) / window_size))		
# 		
#  posterior = mod[[4]]$posterior$estimates[(burnin+1):niter,]		
#		
#  genetic_values = tcrossprod(M, posterior)		
#  genetic_variance = apply(genetic_values,2,var)		
#		
#  res = array(0, dim=c(n_windows,6))		
#  colnames(res) <- c("window","mean_var","mean_var_proportion","prob_var_bigger_threshold","start","end")		
#  end=0	
#  count = 0	
#		
#  for(i in 1:n_windows) {		
#		
#    count = count + 1
#    if(verbose) print(paste("window: ", i, " out of ", n_windows,sep=""))		

#    if(sliding_window) {

#      start = count
#      end = count + window_size - 1

#    } else {

#        start = end + 1		
#        end = end + window_size		
#        if(i == n_windows) end = ncol(M) 

#      }
# 		
#    window_genetic_values = tcrossprod(M[,start:end], posterior[,start:end])		
#    window_genetic_variance = ccolmv(window_genetic_values,compute_var=TRUE)		
#    post_var_proportion = window_genetic_variance / genetic_variance		
#  			
#    res[i,"window"] = i		
#    res[i,"mean_var"] = mean(window_genetic_variance)		
#    res[i,"mean_var_proportion"] = mean(post_var_proportion)		
#    res[i,"prob_var_bigger_threshold"] = sum(post_var_proportion > threshold) / nrow(posterior)		
#    res[i,"start"] = start
#    res[i,"end"] = end	
#		
#  }		
#		
#  return(res)		
#		
#}








