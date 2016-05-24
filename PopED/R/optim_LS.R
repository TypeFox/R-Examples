#' Optimization Using a Line Search Algorithm. 
#' 
#' \code{optim_LS} performs sequential grid search optimization of an arbitrary function with respect 
#' to each of the parameters to be optimized over. 
#' The function works for both discrete and continuous optimization parameters 
#' and allows for box-constraints (by using the upper and lower function arguments) or sets of allowed values (by using the 
#' allowed_values function argument) for all parameters, or on a parameter per parameter basis. 
#' 
#' @param allowed_values A list containing allowed values for each parameter \code{list(par1=c(2,3,4,5,6),par2=c(5,6,7,8))}. 
#' A vector containing allowed values for all parameters is also allowed \code{c(2,3,4,5,6)}.
#' @param par A vector of initial values for the parameters to be optimized over.
#' @param fn A function to be minimized (or maximized), 
#' with first argument the vector of parameters over which minimization is to take place. 
#' It should return a scalar result.
#' @param upper Upper bounds on the parameters. A vector.
#' @param lower Lower bounds on the parameters. A vector.
#' @param line_length The number of different parameter values per parameter to evaluate.  The values are 
#' selected as an evenly spaced grid between the upper and lower bounds. 
#' @param trace Should the algorithm outpur results intermittently.
#' @param maximize Should the function be maximized?  Default is to minimize.
#' @param parallel Should we use parallelize the computation.
#' @param parallel_type Which type of parallelization should be used? 
#' Can be "snow" or "multicore".  "snow"  works on Linux-like systems & Windows. "multicore" works only on 
#' Linux-like systems.  By default this is chosen for you depending on your operating system. 
#' See \code{\link{start_parallel}}.
#' @param num_cores The number of cores to use in the parallelization.  By default  is set to the number 
#'  output from 
#' \code{parallel::detectCores()}. 
#' See \code{\link{start_parallel}}.
#' @param seed The random seed to use in the algorithm,
#' @param replicates_index A vector, the same length as the parameters.  
#' If two values are the same in this vector then the parameters may not assume the same value in the optimization.
#' @param ofv_initial An initial objective function value (OFV).  If not NULL then the initial design is not evaluated
#' and the OFV value is assumed to be this number.
#' @param closed_bounds Are the upper and lower limits open (boundaries not allowed) or closed (boundaries allowed) bounds?
#' @param ... Additional arguments passed to \code{fn} and \code{start_parallel}.
#' 
#' @references \enumerate{
#' \item M. Foracchia, A.C. Hooker, P. Vicini and A. Ruggeri, "PopED, a software fir optimal 
#' experimental design in population kinetics", Computer Methods and Programs in Biomedicine, 74, 2004.
#' \item J. Nyberg, S. Ueckert, E.A. Stroemberg, S. Hennig, M.O. Karlsson and A.C. Hooker, "PopED: An extended, 
#' parallelized, nonlinear mixed effects models optimal design tool",  
#' Computer Methods and Programs in Biomedicine, 108, 2012.
#' }
#' 
#' @family Optimize
#' 
#' @example tests/testthat/examples_fcn_doc/examples_optim_LS.R
#' @export
optim_LS <- function(par,
                     fn,
                     lower=NULL,
                     upper=NULL,
                     allowed_values=NULL,
                     #constraints=NULL,
                     #control=NULL,
                     #loc_fac=4,
                     #no_bounds_sd = par,
                     line_length=50,
                     trace = TRUE,
                     maximize=F, # default is to minimize
                     parallel=F, # T or F or a list of (type, n_cores)
                     parallel_type=NULL,
                     num_cores = NULL,
                     seed=round(runif(1,0,10000000)),
                     replicates_index=seq(1,length(par)), # same value, parameters can not be the same value
                     ofv_initial=NULL,
                     closed_bounds=TRUE, # are upper and lower limits allowed values?
                     ...){
  
  # parallel with # T or F or a list of (type, n_cores)
  # also change with optim_ARS and poped_optim and maybe change to optim_poped
  
  #---------------- start trace
  if((trace)){
    tic(name=".ls_savedTime")
  }
  
  #--------------- checks
  if((is.null(lower) || is.null(upper)) && is.null(allowed_values)){
    stop("At least 'lower' and 'upper' or 'allowed_values' or 'generator' should be supplied.")
  }
  if(length(lower)==1 && length(par)>1) lower <- rep(lower,length(par))
  if(length(upper)==1 && length(par)>1) upper <- rep(upper,length(par))
  if(!is.null(allowed_values)){
    if(!is.list(allowed_values)) allowed_values <- list(allowed_values)
    if(length(allowed_values) == 1 && length(par)>1) allowed_values <- rep(allowed_values,length(par))  
  }
  
  #----------------- initialization 
  itvector <- c() 
  dmfvector <- c()
  counter = 0
  best_changed = FALSE
  

  fn1 <- function(par) fn(par, ...)  
  par_opt <- par
  if(is.null(ofv_initial)){
    ofv_opt <- fn1(par_opt) 
  } else {
    ofv_opt <- ofv_initial
  }
  
  if(!is.null(seed)) set.seed(seed)
  
  # start parallel computing
  if(parallel){
    parallel <- start_parallel(parallel,seed=seed,parallel_type=parallel_type,num_cores=num_cores,...) 
    on.exit(if(parallel && (attr(parallel,"type")=="snow")) parallel::stopCluster(attr(parallel,"cluster")))
  }   
  #if(is.null(iter_chunk)) if(parallel) iter_chunk <- attr(parallel,"cores") else iter_chunk <- 1
  
  # continuous and discrete parameters
  par_type <- rep("cont",length(par))
  if(!is.null(allowed_values)){
    for(k in 1:length(par)){
      if(!is.na(allowed_values[[k]]) && length(allowed_values[[k]]>0)){
        par_type[k] <- "cat"          
      }
    }
  }
  
  #resample <- function(x,size=1) x[sample.int(length(x),size=size)]
  #if(!is.null(par)) iter = iter+1
  par_tried <- list(par) # investigated parameter vectors
  compare <-function(a,b) a < b
  which_compare <- function(x) which.min(x)
  if(maximize){
    compare <-function(a,b) a > b
    which_compare <- function(x) which.max(x)  
  }
  
  
  # ------------ generate new parameter sets
  gen_par_list <- function (par,i, ...) {

    # generate new continuous parameters
    if(par_type[i]=="cont"){
      if(!is.finite(lower[i])) lower[i]=par[i]-1000*par[i]
      if(!is.finite(upper[i])) upper[i]=par[i]+1000*par[i]
      if(closed_bounds){
        par_i_set <- seq(lower[i],upper[i],len=line_length)
      } else {
        par_i_set <- seq(lower[i],upper[i],len=(line_length+2))
        par_i_set <- par_i_set[-c(1,length(par_i_set))]
      }
    } else {
      par_i_set <- allowed_values[[i]]
    }
          
    # handle replicates
    if(length(unique(replicates_index))!=length(par)){
      par_i_set <- par_i_set[par_i_set!=par[replicates_index==replicates_index[i]]]
    }

    # create full parameter list
    par_list <- lapply(par_i_set,function(x){ par[i] <- x; par})

    # remove any design that has already been tested (works but is slow!)
    #par_list <- par_list[unlist(lapply(par_list,function(x){!(list(x) %in% par_tried)}))]
  
    return(par_list)
  } # end function
  
  # Randomize which parameter to investigate first
  par_order <- sample.int(length(par))
  
  #------------- evaluate new parameter sets
  if(trace) cat("\n   Initial parameters:",par_opt, "\n")
  if(trace) cat("   Initial OFV:",ofv_opt, "\n\n")
    
  
  for(i in par_order){
    par_list <- gen_par_list(par_opt,i)
    par_tried <- c(par_tried,par_list)
    if(trace){
      cat("   Searching parameter",i,"\n")
      #pb <- txtProgressBar(min=0, max=n, initial=1, style=3) 
      #setTxtProgressBar(pb, i)
      #close(pb)
    }
    
    if(parallel && (attr(parallel,"type")=="multicore")){
      if(!is.null(seed)) set.seed(seed+i)
      res <- parallel::mclapply(par_list,fn1,mc.cores=attr(parallel, "cores"))
    } else if(parallel && (attr(parallel,"type")=="snow")){
      res <- parallel::parLapply(attr(parallel, "cluster"),par_list,fn1)
    } else {
      res <- lapply(par_list,fn1)  
    }
    
    res <- unlist(res)
    best_index <- which_compare(res)
    ofv <- res[best_index]
    par <- par_list[[best_index]]
    
    if((compare(ofv,ofv_opt) || is.null(ofv_opt))){
      if(trace) cat(sprintf("     Changed from %g to %g",par_opt[i],par[i]))        
      if(trace) cat(sprintf(" ; OFV = %g \n",ofv))        
      par_opt <- par 
      ofv_opt <- ofv
    } else {
      if(trace) cat(paste("     No change","\n"))        
    }
      
  }
  
  #--------- Write results
  if((trace)){
    cat("\n   ")
    toc(name=".ls_savedTime")
    cat("\n   Final OFV = ", ofv_opt, "\n") 
    cat("   Parameters:",par_opt, "\n\n")
  }
  
  return(list( par=par_opt,ofv=ofv_opt)) 
}
