#' Optimization main module for PopED
#' 
#' Optimize the objective function.
#' The function works for both discrete and continuous optimization variables.
#' This function takes information from the PopED database supplied as an argument.
#' The PopED database supplies information about the the model, parameters, design and methods to use.
#' Some of the arguments coming from the PopED database can be overwritten;  
#' if they are supplied then they are used instead of the arguments from the PopED database.
#' 
#' @inheritParams RS_opt
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @inheritParams Dtrace
#' @inheritParams calc_ofv_and_fim
#' @param ... arguments passed to other functions. See \code{\link{Doptim}}.
#' 
#' 
#' @references \enumerate{
#' \item M. Foracchia, A.C. Hooker, P. Vicini and A. Ruggeri, "PopED, a software fir optimal 
#' experimental design in population kinetics", Computer Methods and Programs in Biomedicine, 74, 2004.
#' \item J. Nyberg, S. Ueckert, E.A. Stroemberg, S. Hennig, M.O. Karlsson and A.C. Hooker, "PopED: An extended, 
#' parallelized, nonlinear mixed effects models optimal design tool",  
#' Computer Methods and Programs in Biomedicine, 108, 2012.
#' }
#' @family Optimize
#' 
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_poped_optimize.R
#' @export


poped_optimize <- function(poped.db,
                           ni=NULL, 
                           xt=NULL, 
                           model_switch=NULL, 
                           x=NULL, 
                           a=NULL, 
                           bpop=NULL, 
                           d=NULL, 
                           maxxt=NULL, 
                           minxt=NULL,
                           maxa=NULL,
                           mina=NULL,
                           fmf=0,
                           dmf=0,
                           trflag = TRUE,
                           opt_xt=poped.db$settings$optsw[2],
                           opt_a=poped.db$settings$optsw[4],
                           opt_x=poped.db$settings$optsw[3],
                           opt_samps=poped.db$settings$optsw[1],
                           opt_inds=poped.db$settings$optsw[5],
                           cfaxt=poped.db$settings$cfaxt, 
                           cfaa=poped.db$settings$cfaa,
                           rsit=poped.db$settings$rsit,
                           rsit_output=poped.db$settings$rsit_output,
                           fim.calc.type=poped.db$settings$iFIMCalculationType,
                           ofv_calc_type=poped.db$settings$ofv_calc_type,
                           approx_type=poped.db$settings$iApproximationMethod,
                           bUseExchangeAlgorithm=poped.db$settings$bUseExchangeAlgorithm,
                           iter=1,
                           d_switch=poped.db$settings$d_switch,
                           ED_samp_size = poped.db$settings$ED_samp_size,
                           bLHS=poped.db$settings$bLHS,
                           use_laplace=poped.db$settings$iEDCalculationType,
                           ...){
  
  ## update poped.db with options supplied in function
  called_args <- match.call()
  default_args <- formals()
  for(i in names(called_args)[-1]){
    if(length(grep("^poped\\.db\\$",capture.output(default_args[[i]])))==1) {
      #eval(parse(text=paste(capture.output(default_args[[i]]),"<-",called_args[[i]])))
      eval(parse(text=paste(capture.output(default_args[[i]]),"<-",i)))
    }
  }
  #=========================================
  # PERFORM OPTIMAL DESIGN
  #=========================================
  fmf = 0 #The best FIM so far
  dmf = 0 #The best ofv of FIM  so far
  if((sum(poped.db$settings$optsw)>0)){
    
    # ------------- downsizing of general design
    downsize.list <- downsizing_general_design(poped.db)
    if(is.null(ni)) ni <- downsize.list$ni
    if(is.null(xt)) xt <- downsize.list$xt
    if(is.null(model_switch)) model_switch <- downsize.list$model_switch
    if(is.null(x)) x <- downsize.list$x
    if(is.null(a)) a <- downsize.list$a    
    if(is.null(bpop)) bpop <- downsize.list$bpop
    if(is.null(d)) d <- downsize.list$d
    if(is.null(maxxt)) maxxt <- downsize.list$maxxt
    if(is.null(minxt)) minxt <- downsize.list$minxt
    if(is.null(maxa)) maxa <- downsize.list$maxa
    if(is.null(mina)) mina <- downsize.list$mina
    
    #------------- optimization routine call
    if((opt_samps==TRUE && opt_inds==TRUE)){
      stop('Cannot optimize over both number of samples and number ', 
           'of individuals in each group')
    }
    
    if((sum(poped.db$settings$optsw)==0)){
      stop(sprintf('No optimization parameter is set.'))
    }
    
    # if optimizing over number of samples then do this
    if((opt_samps==TRUE)){
      if((poped.db$settings$d_switch==TRUE) ){#D-optimal
        #         returnArgs <- Dpoptim(poped.db,model_switch,ni,xt,x,a,bpop,d,maxxt,minxt,maxa,mina,fmf,dmf) 
        #         xt <- returnArgs[[1]]
        #         x <- returnArgs[[2]]
        #         a <- returnArgs[[3]]
        #         ni <- returnArgs[[4]]
        #         fmf <- returnArgs[[5]]
        #         dmf <- returnArgs[[6]]
        #         poped.db <- returnArgs[[7]]
        stop('Sample number optimization not implemented in the R-version of PopED yet')
      } else { #ED-optimal
        #         returnArgs <- EDpoptim(poped.db,model_switch,ni,xt,x,a,bpop,d,maxxt,minxt,maxa,mina,fmf,dmf) 
        #         xt <- returnArgs[[1]]
        #         x <- returnArgs[[2]]
        #         a <- returnArgs[[3]]
        #         ni <- returnArgs[[4]]
        #         fmf <- returnArgs[[5]]
        #         dmf <- returnArgs[[6]]
        #         poped.db <- returnArgs[[7]]
        stop('Sample number optimization not implemented in the R-version of PopED yet')
      }
      
      # if optimizing over number of individuals in group then do this
    } else if (opt_inds==TRUE){
      if((poped.db$settings$d_switch==TRUE) ){#D-optimal
        #         returnArgs <- DoptIndGrp(poped.db,model_switch,ni,xt,x,a,bpop,d,maxxt,minxt,maxa,mina,fmf,dmf) 
        #         xt <- returnArgs[[1]]
        #         x <- returnArgs[[2]]
        #         a <- returnArgs[[3]]
        #         groupsize <- returnArgs[[4]]
        #         fmf <- returnArgs[[5]]
        #         dmf <- returnArgs[[6]]
        #         poped.db <- returnArgs[[7]]
        stop('Number of individuals in different groups optimization not implemented in the R-version of PopED yet')
      } else {
        #         returnArgs <- EDoptIndGrp(poped.db,model_switch,ni,xt,x,a,bpop,d,maxxt,minxt,maxa,mina,fmf,dmf) 
        #         xt <- returnArgs[[1]]
        #         x <- returnArgs[[2]]
        #         a <- returnArgs[[3]]
        #         groupsize <- returnArgs[[4]]
        #         fmf <- returnArgs[[5]]
        #         dmf <- returnArgs[[6]]
        #         poped.db <- returnArgs[[7]]
        stop('Number of individuals in different groups optimization not implemented in the R-version of PopED yet')        
      }
      #poped.db$design$groupsize = groupsize
    } else {
      #Using ea algorithm do this,
      if((bUseExchangeAlgorithm==TRUE)){
        #if(opt_xt || opt_a) stop('MFEA algorithm does not work with xt and a right now')
        returnArgs <- mfea(poped.db,model_switch,ni,xt,x,a,bpop,d,maxxt,minxt,maxa,mina,fmf,dmf,...) 
        xt <- returnArgs[[1]]
        x <- returnArgs[[2]]
        a <- returnArgs[[3]]
        fmf <- returnArgs[[4]]
        dmf <- returnArgs[[5]]
        poped.db <- returnArgs[[6]]
      } else {
        if((d_switch==TRUE) ){#D-optimal over continuous variables
          returnArgs <- Doptim(poped.db,ni, xt, model_switch, x, a, bpop, d, maxxt, minxt,maxa,mina,fmf,dmf,
                               trflag=trflag,...) 
          xt <- returnArgs[[1]]
          x <- returnArgs[[2]]
          a <- returnArgs[[3]]
          fmf <- returnArgs[[4]]
          dmf <- returnArgs[[5]]
          poped.db <- returnArgs[[6]]
        } else {
          #if((poped.db$settings$iEDCalculationType>=1)){
          returnArgs <- LEDoptim(poped.db,
                                 model_switch=model_switch,
                                 ni=ni,
                                 xt=xt,
                                 x=x,
                                 a=a,
                                 bpopdescr=bpop,
                                 ddescr=d,
                                 maxxt=maxxt,
                                 minxt=minxt,
                                 maxa=maxa,
                                 mina=mina,
                                 ofv_init=dmf,
                                 fim_init=fmf,
                                 trflag=trflag,
                                 d_switch=d_switch,
                                 use_laplace=use_laplace,
                                 ...)
          xt <- returnArgs[[1]]
          x <- returnArgs[[2]]
          a <- returnArgs[[3]]
          fmf <- returnArgs[[4]]
          dmf <- returnArgs[[5]]
          poped.db <- returnArgs[[6]]
          # stop('E-family optimization using the requested methods not implemented in the R-version of PopED yet')        
          #} else {
          #             returnArgs <- EDoptim(poped.db,model_switch,ni,xt,x,a,bpop,d,maxxt,minxt,maxa,mina,fmf,dmf) 
          #             xt <- returnArgs[[1]]
          #             x <- returnArgs[[2]]
          #             a <- returnArgs[[3]]
          #             fmf <- returnArgs[[4]]
          #             dmf <- returnArgs[[5]]
          #             poped.db <- returnArgs[[6]]
          #  stop('E-family optimization using the requested methods not implemented in the R-version of PopED yet')            
          #}
        }
      }
    }
    
    poped.db$design$ni[1:size(ni,1)]=ni
    #poped.db$design$ni <- ni
    if((!isempty(xt))){
      poped.db$design$xt[1:size(xt,1),1:size(xt,2)]=xt
      #poped.db$design$xt <- xt
    }
    if((!isempty(x))){
      poped.db$design$x[1:size(x,1),1:size(x,2)]=x
      #poped.db$design$x <- x
    }
    if((!isempty(a))){
      poped.db$design$a[1:size(a,1),1:size(a,2)]=a
      #poped.db$design$a <- a
    }
  }
  return(list( fmf= fmf, dmf= dmf, poped.db = poped.db )) 
}

