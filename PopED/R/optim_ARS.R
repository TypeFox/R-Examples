#' Optimization Using Adaptive Random Search. 
#' 
#' Optimize an objective function using an adaptive random search algorithm.  
#' The function works for both discrete and continuous optimization parameters 
#' and allows for box-constraints and sets of allowed values.
#'  
#' @param loc_fac Locality factor for determining the standard deviation of the sampling distribution around the current
#' position of the parameters. The initial standard deviation is normally calculated as \code{(upper - lower)/loc_fac} 
#' except in cases when  there are no upper or lower limits (e.g. when \code{upper=Inf} or \code{lower=-Inf}).
#' @param no_bounds_sd The standard deviation of the sampling distribution around the current
#' position of the parameters when there are no upper or lower limits (e.g. when \code{upper=Inf} or \code{lower=-Inf}).
#' @param adapt_scale The scale for adapting the size of the sampling distribution.  The adaptation of the 
#' standard deviation of the sampling distribution around the current
#' position of the parameters is done after \code{iter_adapt} iteration with no change in the best objective function.  
#' When adapting, the  standard deviation of the sampling distribution is calculated as 
#' \code{(upper - lower)/(loc_fac*ff*adapt_scale)} where ff starts at 1 and increases by 1 for each adaptation.
#' @param allowed_values A list containing allowed values for each parameter \code{list(par1=c(2,3,4,5,6),par2=c(5,6,7,8))}. 
#' A vector containing allowed values for all parameters is also allowed \code{c(2,3,4,5,6)}.
#' @param iter The number of iterations for the algorithm to perfrom (this is a maximum number, it could be less).
#' @param iter_adapt The number of iterations before adapting (shrinking) the parameter search space.
#' @param max_run The maximum number of iterations to run without a change in the best parameter estimates.
#' @param trace_iter How many interations between each update to the screen about the result of the search.
#' @param new_par_max_it The algorithm randomly chooses samples based on the current best set of parameters.  If when drawing 
#' these samples the new parameter set has already been tested then a new draw is performed. After \code{new_par_max_it} draws, with
#' no new parameter sets, then the algorithm stops.
#' @param allow_replicates Should the algorithm allow parameters to have the same value?
#' @param generator A user-defined function that generates new parameter sets to try in the algorithm.  See examples below.
#' 
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
#' @inheritParams optim_LS
#' 
#' @example tests/testthat/examples_fcn_doc/examples_optim_ARS.R
#' @export
optim_ARS <- function(par,
                      fn,
                      lower=NULL,
                      upper=NULL,
                      allowed_values=NULL,
                      #constraints=NULL,
                      #control=NULL,
                      loc_fac=4,
                      no_bounds_sd = par,
                      iter=400,
                      iter_adapt = 50,# number of iterations with no change before adapting search space
                      adapt_scale = 1,
                      max_run=200,
                      trace = TRUE,
                      trace_iter=5,
                      #par_grouping=NULL,
                      new_par_max_it = 200,
                      maximize=F, # default is to minimize
                      parallel=F,
                      parallel_type=NULL,
                      num_cores = NULL,
                      seed=round(runif(1,0,10000000)),
                      allow_replicates=TRUE,
                      generator=NULL,
                      ...){
  
  # constratints
  # uniform instead of normal sampling

  #---------------- start trace
  if((trace)){
    tic(name=".ars_savedTime")
    wd_iter <- nchar(iter) 
  }
  
  #--------------- checks
  if((is.null(lower) || is.null(upper)) && is.null(allowed_values) && is.null(generator)){
    stop("At least 'lower' and 'upper' or 'allowed_values' or 'generator' should be supplied.")
  }   
  if(length(lower)==1 && length(par)>1) lower <- rep(lower,length(par))
  if(length(upper)==1 && length(par)>1) upper <- rep(upper,length(par))
  if(!is.null(allowed_values)){
    if(!is.list(allowed_values)) allowed_values <- list(allowed_values)
    if(length(allowed_values) == 1 && length(par)>1) allowed_values <- rep(allowed_values,length(par))  
  }
  
  #----------------- initialization 
  nullit=1
  runs <- 1
  ff=1
  ff_scale <- 1
  par0=par
  fn1 <- function(par) fn(par, ...)  
  #npar <- max(c(length(lower),length(upper),length(allowed_values)))
  npar <- length(par)
  ofv_opt <- NULL
  par_opt <- par
  dpar=(upper-lower)/loc_fac
  dpar[is.infinite(dpar)] <- no_bounds_sd
  if(!is.null(seed)) set.seed(seed)
  
  if(parallel){
    parallel <- start_parallel(parallel,seed=seed,parallel_type=parallel_type,num_cores=num_cores,...) 
    on.exit(if(parallel && (attr(parallel,"type")=="snow")) parallel::stopCluster(attr(parallel,"cluster")))
  }  
  iter_chunk = NULL
  if(is.null(iter_chunk)) if(parallel) iter_chunk <- attr(parallel,"cores") else iter_chunk <- 1
  
  # continuous and discrete parameters
  par_type <- rep("cont",npar)
  if(!is.null(allowed_values)){
    for(k in 1:npar){
      if(!is.na(allowed_values[[k]]) && length(allowed_values[[k]]>0)){
        par_type[k] <- "cat"          
      }
    }
  }
  
  resample <- function(x,size=1) x[sample.int(length(x),size=size)]
  #if(!is.null(par)) iter = iter+1
  par_vec <- cell(iter,1) # investigated parameter vectors
  compare <-function(a,b) a < b
  if(maximize) compare <-function(a,b) a > b
  
  
  # ------------ generate and evaluate new parameter sets
  gen_par_ofv <- function (it, par, ...) {
    need_new_par <- TRUE
    new_par_it <- 0
    
    #---------------- generate new parameter set
    while(need_new_par && (new_par_it < new_par_max_it)){
      if(it==1 && !is.null(par)){
        need_new_par <- FALSE
        break
      }
      # generate new continuous parameters
      if(is.null(par[par_type=="cont"])){
        par[par_type=="cont"]  <- lower[par_type=="cont"]+(upper[par_type=="cont"] - lower[par_type=="cont"])/2
        par[is.nan(par) && par_type=="cont"] <- 0
      }
      
      generate_par <- function(par,par_type,allowed_values,...){
        par[par_type=="cont"]=par[par_type=="cont"]+dpar[par_type=="cont"]/ff_scale*rnorm(length(par[par_type=="cont"]))
        #par[par_type=="cont" && is.infinite(par)] <- par[par_type=="cont"]+1e10/ff*rnorm(length(par[par_type=="cont"]))
        
        # generate new discrete parameters
        par[par_type=="cat"] <- vapply(allowed_values[par_type=="cat"],resample,c(0)) 
        return(par)
      }
      
      if(!is.null(generator)){
        par <- do.call(generator,list(par,...)) 
      } else {
        par <- generate_par(par,par_type,allowed_values)  
      }
      
      
      # handle replicates
      if(!allow_replicates){
        while(any(duplicated(par))) {
          cur_pars <- par[!duplicated(par)]
          tmp_allowed <- lapply(allowed_values[duplicated(par)],setdiff,cur_pars)
          #if(any(sapply(tmp_allowed,length)==0)) stop("Not enough potential values for all parameters to have different values.")
          par[duplicated(par)] <- generate_par(par[duplicated(par)],par_type[duplicated(par)],tmp_allowed)
        }
      }
      
      # Group samples that should be grouped
      #       if(!is.null(par_grouping)){
      #         for(k in unique(par_grouping)){
      #           par[par_grouping==k] <- resample(par[par_grouping==k])
      #         }
      #       }
      
      # set to boundary if beyond the boundary
      if(!all(is.null(upper))) par[par>upper] <- upper[par>upper]
      #par=par-((par>upper)*(par-upper)) 
      if(!all(is.null(lower))) par[par<lower] <- lower[par<lower]
      #par=par+((par<lower)*(lower-par))
      
      # check if design has already been tested
      #need_new_par <- any(sapply(sapply(par_vec, length)!=0,function(x) all(x == par)))
      #par_vec[sapply(par_vec, length)!=0,]
      need_new_par <- any(sapply(par_vec[sapply(par_vec, length)!=0,],function(x) all(x==par)))
      #browser()
      new_par_it <- new_par_it + 1
    } # end while  
    
    #------------- evaluate new parameter sets
    if(need_new_par){
      ofv <- NULL
      par <- NULL
    } else {
      ofv <- fn1(par)
    }
    return(list(ofv=ofv,par=par,need_new_par=need_new_par))
  } # end function
  
  for(it in 1:ceiling(iter/iter_chunk)){
    start_it <- (it-1)*iter_chunk+1
    stop_it <- min(it*iter_chunk,iter)
    it_seq <- start_it:stop_it
    
    # generate new parameters and OFVs
    if(parallel && (attr(parallel,"type")=="multicore")){
      if(!is.null(seed)) set.seed(seed+it)
      res <- mapply(c,parallel::mclapply(it_seq,gen_par_ofv,par_opt,mc.cores=attr(parallel, "cores")))
    } else if(parallel && (attr(parallel,"type")=="snow")){
      res <- mapply(c,parallel::parLapply(attr(parallel, "cluster"),it_seq,gen_par_ofv,par_opt))
    } else {
      res <- mapply(c,lapply(it_seq,gen_par_ofv,par_opt))  
    }
    
    res2 <- res[,!sapply(res["ofv",],is.null),drop=F]
    
    if(maximize){
      out <- res2[,which.max(res2["ofv",])]  
    } else {
      out <- res2[,which.min(res2["ofv",])]  
    }
    
    # check if a unique new parameter vector was generated
    if(is.null(out$ofv)){
      cat(paste0("Maximum number of duplicate parameter samples reached (new_par_max_it=",new_par_max_it,"), optimization stopped.\n"))
      break
    }
    
    ofv <- out$ofv
    par <- out$par
    
    par_vec[it_seq] <- res["par",] # save new parameter vectors in a list
    
    if((compare(ofv,ofv_opt) || is.null(ofv_opt))){
      par_opt <- par 
      ofv_opt <- ofv
      nullit=1
      runs <- 1
      #ff=1
    } else {
      nullit=nullit+(length(it_seq))
      runs <- runs+(length(it_seq))
    }
    if((nullit>=iter_adapt) ){# when to make the search area smaller
      ff=ff+1
      ff_scale <- ff*adapt_scale
      nullit=1
    }
    
    if((trace && any(rem(it_seq,trace_iter)==0))){
      if(length(it_seq)==1){ 
        cat(sprintf(paste0("It. %",wd_iter,"i"),start_it))
      } else {
        cat(sprintf(paste0("It. %",wd_iter,"i to %",wd_iter,"i"),start_it,stop_it))
      }
      cat(sprintf(" | OFV = %g",ofv_opt))
      if(trace==2) cat(" | opt. par. = ",par_opt)
      #cat(" | runs = ",runs)
      #cat(" | nullit = ",nullit)
      if(trace==3) cat(" | par tried = ",par)
      cat("\n")
    }
    
    if(runs>=max_run){
      cat(paste0("Maximum number of identical optimal values reached (max_run=",max_run,"), optimization stopped.\n"))
      break
    }     
  }
  
  #--------- Write results
  if((trace)){
    cat("\nTotal iterations:",stop_it,"\n")
    toc(name=".ars_savedTime")
    cat("\nFinal OFV = ", ofv_opt, "\n") 
    cat("Parameters:",par_opt, "\n\n")
  }
  
  return(list( par=par_opt,ofv=ofv_opt)) 
}
