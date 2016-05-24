#' Optimization main module for PopED
#' 
#' Optimize the objective function. The function works for both discrete and 
#' continuous optimization variables. If more than one optimization method is 
#' specified then the methods are run in series.  If \code{loop_methods=TRUE} 
#' then the series of optimization methods will be run for \code{iter_max} 
#' iterations, or until the efficiency of the design after the current series 
#' (compared to the start of the series) is less than \code{eff_crit}.
#' 
#' This function takes information from the PopED database supplied as an 
#' argument. The PopED database supplies information about the the model, 
#' parameters, design and methods to use. Some of the arguments coming from the 
#' PopED database can be overwritten; if they are supplied then they are used 
#' instead of the arguments from the PopED database.
#' 
#' @inheritParams RS_opt
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @inheritParams Dtrace
#' @inheritParams calc_ofv_and_fim
#' @inheritParams optim_LS
#' @param ... arguments passed to other functions.
#' @param control Contains control arguments for each method specified.
#' @param method A vector of optimization methods to use in a sequential 
#'   fashion.  Options are \code{c("ARS","BFGS","LS","GA")}. \code{c("ARS")} is 
#'   for Adaptive Random Search \code{\link{optim_ARS}}.  \code{c("LS")} is for 
#'   Line Search \code{\link{optim_LS}}. \code{c("BFGS")} is for Method 
#'   "L-BFGS-B" from \code{\link[stats]{optim}}. \code{c("GA")} is for the 
#'   genetic algorithm from \code{\link[GA]{ga}}.
#' @param out_file Save output from the optimization to a file.
#' @param loop_methods Should the optimization methods be looped for
#'   \code{iter_max} iterations, or until the efficiency of the design after the
#'   current series (compared to the start of the series) is less than, or equal to,
#'   \code{eff_crit}?
#' @param eff_crit If \code{loop_methods==TRUE}, the looping will stop if the
#'   efficiency of the design after the current series (compared to the start of
#'   the series) is less than, or equal to, \code{eff_crit}.
#'   
#'   
#'   
#' @references \enumerate{ \item M. Foracchia, A.C. Hooker, P. Vicini and A. 
#'   Ruggeri, "PopED, a software fir optimal experimental design in population 
#'   kinetics", Computer Methods and Programs in Biomedicine, 74, 2004. \item J.
#'   Nyberg, S. Ueckert, E.A. Stroemberg, S. Hennig, M.O. Karlsson and A.C. 
#'   Hooker, "PopED: An extended, parallelized, nonlinear mixed effects models 
#'   optimal design tool", Computer Methods and Programs in Biomedicine, 108, 
#'   2012. }
#'   
#' @family Optimize
#'   
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_poped_optim.R
#' @export

poped_optim <- function(poped.db,
                        opt_xt=poped.db$settings$optsw[2],
                        opt_a=poped.db$settings$optsw[4],
                        opt_x=poped.db$settings$optsw[3],
                        opt_samps=poped.db$settings$optsw[1],
                        opt_inds=poped.db$settings$optsw[5],
                        method=c("ARS","BFGS","LS"),
                        control=list(),
                        trace = TRUE,
                        fim.calc.type=poped.db$settings$iFIMCalculationType,
                        ofv_calc_type=poped.db$settings$ofv_calc_type,
                        approx_type=poped.db$settings$iApproximationMethod,
                        d_switch=poped.db$settings$d_switch,
                        ED_samp_size = poped.db$settings$ED_samp_size,
                        bLHS=poped.db$settings$bLHS,
                        use_laplace=poped.db$settings$iEDCalculationType,
                        out_file="",
                        parallel=F,
                        parallel_type=NULL,
                        num_cores = NULL,
                        loop_methods=ifelse(length(method)>1,TRUE,FALSE),
                        iter_max = 10,
                        eff_crit = 1.001,
                        ...){
  
  #------------ update poped.db with options supplied in function
  called_args <- match.call()
  default_args <- formals()
  for(i in names(called_args)[-1]){
    if(length(grep("^poped\\.db\\$",capture.output(default_args[[i]])))==1) {
      #eval(parse(text=paste(capture.output(default_args[[i]]),"<-",called_args[[i]])))
      eval(parse(text=paste(capture.output(default_args[[i]]),"<-",i)))
    }
  }
  
  #----------- checks
  if((sum(poped.db$settings$optsw)==0)){
    stop('No optimization parameter is set.')
  }
  
  #------------- initialization
  fmf = 0 #The best FIM so far
  dmf = 0 #The best ofv of FIM  so far
  #output <-calc_ofv_and_fim(poped.db,...)
  output <-calc_ofv_and_fim(poped.db,d_switch=d_switch,
                   ED_samp_size=ED_samp_size,
                   bLHS=bLHS,
                   use_laplace=use_laplace,
                   ofv_calc_type=ofv_calc_type,
                   fim.calc.type=fim.calc.type,
                   ...)
  
  fmf <- output$fim
  dmf <- output$ofv
  fmf_init <- fmf
  dmf_init <- dmf
  
  #--------------------- write out info to a file
  fn=blockheader(poped.db,name="optim",e_flag=!d_switch,
                 fmf=fmf_init,dmf=dmf_init,
                 out_file=out_file,
                 trflag=trace,
                 ...)
  
  # Collect the parameters to optimize
  par <- c()
  upper <- c()
  lower <- c()
  par_grouping <- c()
  par_type <- c()
  par_dim <- list()
  allowed_values <- NULL
  build_allowed_values <- FALSE
  if(!is.null(poped.db$design_space$xt_space) ||
     !is.null(poped.db$design_space$a_space)) build_allowed_values <- TRUE
  if(opt_samps) stop('Sample number optimization is not yet implemented in the R-version of PopED.')
  if(opt_inds) stop('Optimization  of number of individuals in different groups is not yet implemented in the R-version of PopED.')
  if(opt_xt){ 
    par <- c(par,poped.db$design$xt)
    upper <- c(upper,poped.db$design_space$maxxt)
    lower <- c(lower,poped.db$design_space$minxt)
    par_grouping <- c(par_grouping,poped.db$design_space$G_xt)
    par_type <- c(par_type,rep("xt",length(poped.db$design$xt)))
    par_dim$xt <- dim(poped.db$design$xt)
    if(is.null(poped.db$design_space$xt_space) && build_allowed_values){ 
      poped.db$design_space$xt_space <- cell(par_dim$xt)
    }
    allowed_values <- c(allowed_values,poped.db$design_space$xt_space)
    
  }
  if(opt_a) { 
    par <- c(par,poped.db$design$a)
    upper <- c(upper,poped.db$design_space$maxa)
    lower <- c(lower,poped.db$design_space$mina)
    if(opt_xt){
      par_grouping <- c(par_grouping,poped.db$design_space$G_a + max(par_grouping)) 
    } else {
      par_grouping <- c(par_grouping,poped.db$design_space$G_a) 
    }
    par_type <- c(par_type,rep("a",length(poped.db$design$a)))
    par_dim$a <- dim(poped.db$design$a)
    if(is.null(poped.db$design_space$a_space) && build_allowed_values){ 
      poped.db$design_space$a_space <- cell(par_dim$a)
    }
    allowed_values <- c(allowed_values,poped.db$design_space$a_space)
    
  }
  if(opt_x) NULL # par <- c(par,poped.db$design$x)
  
  # continuous and discrete parameters
  npar <- max(c(length(lower),length(upper),length(allowed_values),length(par)))
  par_cat_cont <- rep("cont",npar)
  if(!is.null(allowed_values)){
    for(k in 1:npar){
      if(!is.na(allowed_values[[k]]) && length(allowed_values[[k]]>0)){
        par_cat_cont[k] <- "cat"          
      }
    }
  }
  
  # Parameter grouping
  par_df <- data.frame(par,par_grouping,upper,lower,par_type,par_cat_cont)
  par_df_unique <- NULL
  if(!all(!duplicated(par_df$par_grouping))){
    par_df_unique <- par_df[!duplicated(par_df$par_grouping),]
    par <- par_df_unique$par
    lower <- par_df_unique$lower
    upper <- par_df_unique$upper
    par_cat_cont <- par_df_unique$par_cat_cont
  }
  
  par_df_2 <- data.frame(par,upper,lower,par_cat_cont)
  par_fixed_index <- which(upper==lower)
  par_not_fixed_index <- which(upper!=lower)
  if(length(par_fixed_index)!=0){
    par <- par[-c(par_fixed_index)]
    lower <- lower[-c(par_fixed_index)]
    upper <- upper[-c(par_fixed_index)]
    par_cat_cont <- par_cat_cont[-c(par_fixed_index)]
  }
  
  #------- create optimization function with optimization parameters first
  ofv_fun <- function(par,only_cont=F,...){
    
    if(length(par_fixed_index)!=0){
      par_df_2[par_not_fixed_index,"par"] <- par
      par <- par_df_2$par
    }
    
    if(!is.null(par_df_unique)){
      if(only_cont){ 
        par_df_unique[par_df_unique$par_cat_cont=="cont","par"] <- par
      } else {
        par_df_unique$par <- par
      }
      for(j in par_df_unique$par_grouping){
        par_df[par_df$par_grouping==j,"par"] <- par_df_unique[par_df_unique$par_grouping==j,"par"]
      }  
      
      #par_full[par_cat_cont=="cont"] <- par 
      par <- par_df$par
    } else if (only_cont){ 
      par_df[par_df$par_cat_cont=="cont","par"] <- par
      par <- par_df$par
    }
    xt <- NULL
    if(opt_xt) xt <- matrix(par[par_type=="xt"],par_dim$xt)
    a <- NULL
    if(opt_a) a <- matrix(par[par_type=="a"],par_dim$a)
    
    if(d_switch){
      FIM <- evaluate.fim(poped.db,xt=xt,a=a,...)
      ofv <- ofv_fim(FIM,poped.db,...)
    } else{
      output <-calc_ofv_and_fim(poped.db,d_switch=d_switch,
                                ED_samp_size=ED_samp_size,
                                bLHS=bLHS,
                                use_laplace=use_laplace,
                                ofv_calc_type=ofv_calc_type,
                                fim.calc.type=fim.calc.type,
                                xt=xt,
                                a=a,
                                ...)
      
      FIM <- output$fim
      ofv <- output$ofv
    }
              
        
    #ofv <- tryCatch(ofv_fim(FIM,poped.db,...), error = function(e) e)
    if(!is.finite(ofv) && ofv_calc_type==4){
      ofv <- -Inf 
    } else {
      if(!is.finite(ofv)) ofv <- 1e-15
      #if(!is.finite(ofv)) ofv <- NA
      #if(!is.finite(ofv)) ofv <- -Inf
    }
    
    #cat(ofv,"\n")
    return(ofv)
  }
  
  #------------ optimize
  if(!(fn=="")) sink(fn, append=TRUE, split=TRUE)
  
  iter <- 0
  stop_crit <- FALSE
  output$ofv <- dmf_init
  while(stop_crit==FALSE && iter < iter_max){
    ofv_init <- output$ofv
    iter=iter+1
    method_loop <- method
    if(loop_methods){
      cat("************* Iteration",iter," for all optimization methods***********************\n\n") 
    }
    
    while(length(method_loop)>0){
      cur_meth <- method_loop[1]
      method_loop <- method_loop[-1]
      if(cur_meth=="ARS"){
        cat("*******************************************\n")
        cat("Running Adaptive Random Search Optimization\n")
        cat("*******************************************\n")
        
        # handle control arguments
        con <- list(trace = trace, 
                    parallel=parallel,
                    parallel_type=parallel_type,
                    num_cores = num_cores)
        nmsC <- names(con)
        con[(namc <- names(control$ARS))] <- control$ARS
        #if (length(noNms <- namc[!namc %in% nmsC])) warning("unknown names in control: ", paste(noNms, collapse = ", "))
        
        output <- do.call(optim_ARS,c(list(par=par,
                                           fn=ofv_fun,
                                           lower=lower,
                                           upper=upper,
                                           allowed_values = allowed_values,
                                           maximize=T
                                           #par_df_full=par_df
        ),
        #par_grouping=par_grouping),
        con,
        ...))
        
      }
      if(cur_meth=="LS"){
        cat("*******************************************\n")
        cat("Running Line Search Optimization\n")
        cat("*******************************************\n")
        
        # handle control arguments
        con <- list(trace = trace, 
                    parallel=parallel,
                    parallel_type=parallel_type,
                    num_cores = num_cores)
        nmsC <- names(con)
        con[(namc <- names(control$LS))] <- control$LS
        #if (length(noNms <- namc[!namc %in% nmsC])) warning("unknown names in control: ", paste(noNms, collapse = ", "))
        
        output <- do.call(optim_LS,c(list(par=par,
                                          fn=ofv_fun,
                                          lower=lower,
                                          upper=upper,
                                          allowed_values = allowed_values,
                                          maximize=T
                                          #par_df_full=par_df
        ),
        #par_grouping=par_grouping),
        con,
        ...))
        
      }
      if(cur_meth=="BFGS"){
        cat("*******************************************\n")
        cat("Running L-BFGS-B Optimization\n")
        cat("*******************************************\n")
        
        if(trace) trace_optim=3
        if(is.numeric(trace)) trace_optim = trace
        #if(trace==2) trace_optim = 4
        #if(trace==3) trace_optim = 5
        #if(trace==4) trace_optim = 6
        
        # handle control arguments
        con <- list(trace=trace_optim,fnscale=-1)
        nmsC <- names(con)
        con[(namc <- names(control$BFGS))] <- control$BFGS
        #if (length(noNms <- namc[!namc %in% nmsC])) warning("unknown names in control: ", paste(noNms, collapse = ", "))
        
        par_full <- par
        output <- optim(par=par[par_cat_cont=="cont"],
                        fn=ofv_fun,
                        gr=NULL,
                        #par_full=par_full,
                        only_cont=T,
                        lower=lower[par_cat_cont=="cont"],
                        upper=upper[par_cat_cont=="cont"],
                        method = "L-BFGS-B",
                        control=con)
        output$ofv <- output$value
        par_tmp <- output$par
        output$par <- par_full
        output$par[par_cat_cont=="cont"] <- par_tmp
        
        fprintf('\n')
        if(fn!="") fprintf(fn,'\n')
      }
      
      if(cur_meth=="GA"){
        if (!requireNamespace("GA", quietly = TRUE)) {
          stop("GA package needed for this function to work. Please install it.",
               call. = FALSE)
        }
        cat("*******************************************\n")
        cat("Running Genetic Algorithm (GA) Optimization\n")
        cat("*******************************************\n")
        
        # handle control arguments
        parallel_ga <- parallel
        if(!is.null(num_cores))  parallel_ga <- num_cores
        if(!is.null(parallel_type))  parallel_ga <- parallel_type
        
        con <- list(parallel=parallel_ga)
        
        nmsC <- names(con)
        con[(namc <- names(control$GA))] <- control$GA
        #if (length(noNms <- namc[!namc %in% nmsC])) warning("unknown names in control: ", paste(noNms, collapse = ", "))
        
        par_full <- par
        output_ga <- do.call(GA::ga,c(list(type = "real-valued", 
                                           fitness = ofv_fun,
                                           #par_full=par_full,
                                           only_cont=T,
                                           min=lower[par_cat_cont=="cont"],
                                           max=upper[par_cat_cont=="cont"],
                                           suggestions=par[par_cat_cont=="cont"]),
                                      #allowed_values = allowed_values),
                                      con,
                                      ...))
        
        
        output$ofv <- output_ga@fitnessValue
        
        output$par <- output_ga@solution
        
        par_tmp <- output$par
        output$par <- par_full
        output$par[par_cat_cont=="cont"] <- par_tmp
        
        fprintf('\n')
        if(fn!="") fprintf(fn,'\n')
      }
      par <- output$par
    }
    
    if(!loop_methods){
      stop_crit <- TRUE
    } else {
      cat("*******************************************\n")
      cat("Stopping criteria testing\n")
      cat("*******************************************\n")
      
      # stop based on efficiency
      p = get_fim_size(poped.db)
      eff = ofv_criterion(output$ofv,p,poped.db)/ofv_criterion(ofv_init,p,poped.db)
      cat("Efficiency of design (start of loop / end of loop) = ",eff, "\n")
      cat("Efficiency stopping criteria (lower limit) = ",eff_crit, "\n")
      if(eff<=eff_crit) stop_crit <- TRUE
      
      # stop based on change in parameters
      
      # stop based on own function
      if(stop_crit){
        cat("Stopping criteria achieved\n")
      } else {
        cat("Stopping criteria NOT achieved\n")
      }
      cat("\n")
    }
    
  } # end of total loop 
  
  if(!(fn=="")) sink()
  
  # add the results into a poped database 
  # expand results to full size 
  if(length(par_fixed_index)!=0){
    par_df_2[par_not_fixed_index,"par"] <- par
    par <- par_df_2$par
  }
  if(!is.null(par_df_unique)){
    par_df_unique$par <- par
    for(j in par_df_unique$par_grouping){
      par_df[par_df$par_grouping==j,"par"] <- par_df_unique[par_df_unique$par_grouping==j,"par"]
    }  
  } else {
    par_df$par <- par
  }  
  
  #poped.db$design$ni <- ni
  if(opt_xt) poped.db$design$xt[,]=matrix(par_df[par_type=="xt","par"],par_dim$xt)
  if(opt_a) poped.db$design$a[,]=matrix(par_df[par_type=="a","par"],par_dim$a)
  #   if((!isempty(x))){
  #     poped.db$design$x[1:size(x,1),1:size(x,2)]=x
  #     #poped.db$design$x <- x
  #   }
  #   if((!isempty(a))){
  #     poped.db$design$a[1:size(a,1),1:size(a,2)]=a
  #     #poped.db$design$a <- a
  #   }
  
  #--------- Write results
  #if((trflag)){
  #  if(footer_flag){
  #FIM <- evaluate.fim(poped.db,...)
  if(d_switch){
    FIM <- evaluate.fim(poped.db,...)
  } else{
    out <-calc_ofv_and_fim(poped.db,d_switch=d_switch,
                              ED_samp_size=ED_samp_size,
                              bLHS=bLHS,
                              use_laplace=use_laplace,
                              ofv_calc_type=ofv_calc_type,
                              fim.calc.type=fim.calc.type,
                              ...)
    
    FIM <- out$fim
  }
  blockfinal(fn=fn,fmf=FIM,
             dmf=output$ofv,
             groupsize=poped.db$design$groupsize,
             ni=poped.db$design$ni,
             xt=poped.db$design$xt,
             x=poped.db$design$x,
             a=poped.db$design$a,
             model_switch=poped.db$design$model_switch,
             poped.db$parameters$param.pt.val$bpop,
             poped.db$parameters$param.pt.val$d,
             poped.db$parameters$docc,
             poped.db$parameters$param.pt.val$sigma,
             poped.db,
             fmf_init=fmf_init,
             dmf_init=dmf_init,
             ...)
  
  
  #  }
  #}
  
  return(invisible(list( ofv= output$ofv, FIM=FIM, poped.db = poped.db ))) 
}

