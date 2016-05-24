#' Optimize the objective function using an adaptive random search algorithm for D-family and E-family designs. 
#' 
#' Optimize the objective function using an adaptive random search algorithm.  
#' Optimization can be performed for both D-family and E-family designs.
#' The function works for both discrete and continuous optimization variables.
#' This function takes information from the PopED database supplied as an argument.
#' The PopED database supplies information about the the model, parameters, design and methods to use.
#' Some of the arguments coming from the PopED database can be overwritten;  
#' by default these arguments are \code{NULL} in the 
#' function, if they are supplied then they are used instead of the arguments from the PopED database.
#' 
#' @inheritParams evaluate.fim
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @inheritParams calc_ofv_and_fim
#' @param cfaxt First step factor for sample times 
#' @param opt_xt Should the sample times be optimized?
#' @param opt_a Should the continuous design variables be optimized?
#' @param opt_x Should the discrete design variables be optimized?
#' @param approx_type Approximation method for model, 0=FO, 1=FOCE, 2=FOCEI, 3=FOI.
#' @param iter The number of iterations entered into the \code{blockheader_2} function.
#' @param header_flag Should the header text be printed out?
#' @param footer_flag Should the footer text be printed out?
#' @param out_file Which file should the output be directed to?  A string, a file handle using 
#'        \code{\link{file}} or \code{""} will output to the screen.
#' @param compute_inv should the inverse of the FIM be used to compute expected RSE values?  Often not needed
#'        except for diagnostic purposes.
#' @param ... arguments passed to \code{\link{evaluate.fim}} and \code{\link{ofv_fim}}.
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
#' @example tests/testthat/examples_fcn_doc/warfarin_ed.R
#' @example tests/testthat/examples_fcn_doc/examples_RS_opt_gen.R
#' @export
RS_opt_gen <- function(poped.db,
                       ni=NULL, 
                       xt=NULL, 
                       model_switch=NULL, 
                       x=NULL, 
                       a=NULL, 
                       bpopdescr=NULL, 
                       ddescr=NULL, 
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
                       cfaxt=poped.db$settings$cfaxt, 
                       cfaa=poped.db$settings$cfaa,
                       rsit=poped.db$settings$rsit,
                       rsit_output=poped.db$settings$rsit_output,
                       fim.calc.type=poped.db$settings$iFIMCalculationType,
                       approx_type=poped.db$settings$iApproximationMethod,
                       iter=NULL,
                       d_switch=poped.db$settings$d_switch,
                       use_laplace=poped.db$settings$iEDCalculationType,
                       laplace.fim=FALSE,
                       header_flag=TRUE,
                       footer_flag=TRUE,
                       out_file=NULL,
                       compute_inv=TRUE,
                       ...){
  
  
  # Only get inputs that are needed, not double inputs
  # needed inputs to function: get first then run function
  # poped.db$settings$cfaxt     0.001
  # poped.db$design$m
  # max(poped.db$design_space$maxni)
  # poped.db$settings$optsw
  # poped.db$settings$cfaa
  # size(poped.db$design$a,2)
  # poped.db$parameters$covd
  # poped.db$parameters$covdocc
  # poped.db$parameters$docc
  # poped.db$design$groupsize
  # poped.db$parameters$sigma
  # poped.db$settings$bShowGraphs
  # poped.db$settings$maxrsnullit # when to make the search area smaller
  # poped.db$settings$strIterationFileName iteration file if not empty string
  
  ## update poped.db with options supplied in function
  called_args <- match.call()
  default_args <- formals()
  for(i in names(called_args)[-1]){
    if(length(grep("^poped\\.db\\$",capture.output(default_args[[i]])))==1) {
      #eval(parse(text=paste(capture.output(default_args[[i]]),"<-",called_args[[i]])))
      eval(parse(text=paste(capture.output(default_args[[i]]),"<-",i)))
    }
  }
  
  downsize.list <- downsizing_general_design(poped.db)
  if(is.null(ni)) ni <- downsize.list$ni
  if(is.null(xt)) xt <- downsize.list$xt
  if(is.null(model_switch)) model_switch <- downsize.list$model_switch
  if(is.null(x)) x <- downsize.list$x
  if(is.null(a)) a <- downsize.list$a    
  if(is.null(bpopdescr)) bpopdescr <- downsize.list$bpop
  if(is.null(ddescr)) ddescr <- downsize.list$d
  if(is.null(maxxt)) maxxt <- downsize.list$maxxt
  if(is.null(minxt)) minxt <- downsize.list$minxt
  if(is.null(maxa)) maxa <- downsize.list$maxa
  if(is.null(mina)) mina <- downsize.list$mina
  
  # ----------------- initialization of optimum variables 
  xtopt=xt
  xopt=x
  aopt=a
  
  if(sum(opt_x,opt_a,opt_xt)==0){
    cat("No optimization variables specified in input\n")
    return(invisible(list( xtopt= xtopt,xopt=xopt,aopt=aopt,fmf=fmf,dmf=dmf,poped.db =poped.db )))   
  }  
  
  loops=TRUE # for building a RS without for loops use FALSE (developmental)
  
  # ----------------- initialization of size variables
  m=size(ni,1)
  maxni=size(xt,2)
  na=size(a,2)
  
  # ----------------- type of optimization determination
  axt=opt_xt*cfaxt*matrix(1,m,maxni)
  aa=opt_a*cfaa*matrix(1,m,na)
  optxt=opt_xt
  optx=opt_x
  opta=opt_a
  
  # ----------------- initialization of model parameters
  bpop=bpopdescr[,2,drop=F]
  d=getfulld(ddescr[,2,drop=F],poped.db$parameters$covd)
  docc_full = getfulld(poped.db$parameters$docc[,2,drop=F],poped.db$parameters$covdocc)
  
  output <- calc_ofv_and_fim(poped.db,
                             ofv=dmf,
                             fim=fmf, 
                             d_switch=d_switch,  
                             bpopdescr=bpopdescr, 
                             ddescr=ddescr,
                             bpop=bpop, 
                             d=d, 
                             docc_full= docc_full, 
                             model_switch=model_switch, 
                             ni=ni, 
                             xt=xt, 
                             x=x, 
                             a=a, 
                             fim.calc.type=fim.calc.type,
                             use_laplace=use_laplace,
                             laplace.fim=laplace.fim, 
                             ...)
  fmf <- output$fim
  dmf <- output$ofv
  fmf_init <- fmf
  dmf_init <- dmf
  
  #   # ------------------ Write summary output file header
  #   param_cvs_init=NULL
  #   if((trflag && header_flag)){
  #     fn=blockheader_2(poped.db,name='RS',iter=iter,e_flag=!d_switch,
  #                      opt_xt=opt_xt,opt_a=opt_a,opt_x=opt_x,
  #                      opt_inds=F,opt_samps=F,
  #                      fmf=fmf,dmf=dmf,
  #                      bpop=bpopdescr,d=ddescr,docc=poped.db$parameters$docc,sigma=poped.db$parameters$sigma,
  #                      out_file=out_file,compute_inv=compute_inv)
  #     if(is.matrix(fmf) && compute_inv){ 
  #       param_vars_init=diag_matlab(inv(fmf))
  #       returnArgs <-  get_cv(param_vars_init,bpop=bpopdescr,d=ddescr,docc=poped.db$parameters$docc,sigma=poped.db$parameters$sigma,poped.db) 
  #       params_init <- returnArgs[[1]]
  #       param_cvs_init <- returnArgs[[2]]
  #     }
  #   } 
  #   
  #--------------------- write out info to a file
  fn=blockheader(poped.db,name="RS",e_flag=!d_switch,iter=iter,
                   opt_xt=opt_xt,opt_a=opt_a,opt_x=opt_x,
                   opt_inds=F,opt_samps=F,
                   fmf=fmf,dmf=dmf,
                   bpop=bpopdescr,d=ddescr,docc=poped.db$parameters$docc,sigma=poped.db$parameters$sigma,
                   out_file=out_file,
                   trflag=trflag,
                   header_flag=header_flag,
                   compute_inv=compute_inv,
                   ...)
  
  
  #   #=====open a file to write to if not already open
  #   if(!exists("fn")){
  #     fn <- out_file
  #     if(is.null(fn)) fn <- ''
  #     if(!any(class(fn)=="file") && !(fn=='')){
  #       fn=file(fn,'w')
  #       if(fn==-1){
  #         stop(sprintf('output file could not be opened'))
  #       }
  #     } 
  #   }
  
  # if((poped.db$settings$bShowGraphs && trflag)){
  #     figure(1)
  #     if((poped.db$settings$Engine$Type==1)){
  #         set(1,'Name','Random Search')
  #     }
  # }
  
  itvector <- c()
  dmfvector <- c()
  
  
  
  # ----------------- RS variables initialization
  dxt=(maxxt-minxt)/poped.db$settings$rslxt
  da=(maxa-mina)/poped.db$settings$rsla
  xtoptn=xtopt
  xoptn=xopt
  aoptn=aopt
  nullit=1
  ff=1
  
  
  
  # ----------------- RANDOM SEARCH BEGINS HERE
  itvector = zeros(0,ceiling(rsit/rsit_output)+1)
  dmfvector = zeros(0,ceiling(rsit/rsit_output)+1)
  dmfvector[1] = dmf
  
  if((trflag)){
    # NEEDS TO BE FIXED with right output
    Dtrace(fn,0,ni,xtopt,xopt,aopt,matrix(0,0,0),matrix(0,0,0),dmf,matrix(0,0,0),
           matrix(0,0,0),matrix(0,0,0),itvector,dmfvector,poped.db)
  }
  
  if((poped.db$settings$parallel$bParallelRS)){
    
    # Generate the input designs
    
    designsin = cell(1,0)
    for(it in 1:rsit){
      if((optxt==TRUE)){
        if((poped.db$design_space$bUseGrouped_xt)){
          xtoptn=grouped_rand(poped.db$design_space$G_xt,xtopt,dxt,ff,axt)
        } else {
          xtoptn=xtopt+dxt/ff*randn(m,maxni)*(axt>0)
        }
        xtoptn=xtoptn-((xtoptn>maxxt)*(xtoptn-maxxt))
        xtoptn=xtoptn+((xtoptn<minxt)*(minxt-xtoptn))
      }
      if((optx==TRUE)){
        xoptn=get_discrete_x(poped.db$design_space$G_x,poped.db$design_space$discrete_x,poped.db$design_space$bUseGrouped_x)
      }
      if((opta==TRUE)){
        if((poped.db$design_space$bUseGrouped_a)){
          aoptn=grouped_rand_a(poped.db$design_space$G_a,aopt,da,ff,aa)
        } else {
          aoptn=aopt+da/ff*randn(m,size(poped.db$design$a,2))*(aa>0)
        }
        aoptn=aoptn-((aoptn>maxa)*(aoptn-maxa))
        aoptn=aoptn+((aoptn<mina)*(mina-aoptn))
      }
      designsin = update_designinlist(designsin,poped.db$design$groupsize,ni,xtoptn,xoptn,aoptn,-1,0)
    }
    
    stop("Parallel execution not yet implemented in PopED for R")
    designsout = designsin
    #designsout = execute_parallel(designsin,poped.db)
    #Store the optimal design
    for(it in 1:rsit){
      if((designsout[[it]]$ofv>dmf)){
        if((optxt==TRUE)){
          xtopt=designsin[[it]]$xt
        }
        if((optx==TRUE)){
          xopt=designsin[[it]]$x
        }
        if((opta==TRUE)){
          aopt=designsin[[it]]$a
        }
        dmf=designsout[[it]]$ofv
        fmf=designsout[[it]]$FIM
      }
      
      if((trflag && (rem(it,rsit_output)==0 || it==rsit))){
        itvector[ceiling(it/rsit_output)+1]=it
        dmfvector[ceiling(it/rsit_output)+1]=dmf
        Dtrace(fn,it,ni,xtopt,xopt,aopt,matrix(0,0,0),matrix(0,0,0),dmf,matrix(0,0,0),
               matrix(0,0,0),matrix(0,0,0),itvector,dmfvector,poped.db)
      }
    }
  } else { # end paralell settings
    if(loops){
      tic()
      for(it in 1:rsit){
        if((optxt==TRUE)){
          if((poped.db$design_space$bUseGrouped_xt)){
            xtoptn=grouped_rand(poped.db$design_space$G_xt,xtopt,dxt,ff,axt)
          } else {
            xtoptn=xtopt+dxt/ff*randn(m,maxni)*(axt>0)
          }
          xtoptn=xtoptn-((xtoptn>maxxt)*(xtoptn-maxxt))
          xtoptn=xtoptn+((xtoptn<minxt)*(minxt-xtoptn))
        }
        if((optx==TRUE)){
          xoptn=get_discrete_x(poped.db$design_space$G_x,poped.db$design_space$discrete_x,poped.db$design_space$bUseGrouped_x)
        }
        if((opta==TRUE)){
          if((poped.db$design_space$bUseGrouped_a)){
            aoptn=grouped_rand_a(poped.db$design_space$G_a,aopt,da,ff,aa)
          } else {
            aoptn=aopt+da/ff*randn(m,size(poped.db$design$a,2))*(aa>0)
          }
          aoptn=aoptn-((aoptn>maxa)*(aoptn-maxa))
          aoptn=aoptn+((aoptn<mina)*(mina-aoptn))
        }
        
        output <- calc_ofv_and_fim(poped.db,
                                   d_switch=d_switch,  
                                   bpopdescr=bpopdescr, 
                                   ddescr=ddescr,
                                   bpop=bpop, 
                                   d=d, 
                                   docc_full= docc_full, 
                                   model_switch=model_switch, 
                                   ni=ni, 
                                   xt=xtoptn, 
                                   x=xoptn, 
                                   a=aoptn, 
                                   fim.calc.type=fim.calc.type,
                                   use_laplace=use_laplace,
                                   laplace.fim=FALSE, 
                                   ...)
        nfmf <- output$fim
        ndmf <- output$ofv            
        if(is.nan(ndmf)) ndmf  <- 0
        if((ndmf>dmf)){
          if((optxt==TRUE)){
            xtopt=xtoptn
          }
          if((optx==TRUE)){
            xopt=xoptn
          }
          if((opta==TRUE)){
            aopt=aoptn
          }
          dmf=ndmf
          fmf=nfmf
          nullit=1
          ff=1
        } else {
          nullit=nullit+1
        }
        if((nullit==poped.db$settings$maxrsnullit) ){# when to make the search area smaller
          ff=ff+1
          nullit=1
        }
        
        if((!isempty(poped.db$settings$strIterationFileName))){
          write_iterationfile('Random Search',it,xtopt,aopt,xopt,ni,poped.db$design$groupsize,fmf,dmf,poped.db)
        }
        
        if((trflag && (rem(it,rsit_output)==0 || trflag && it==rsit))){
          itvector[ceiling(it/rsit_output)+1]=it
          dmfvector[ceiling(it/rsit_output)+1]=dmf
          # fix so that the iterations are output to summary file.
          Dtrace(fn,it,ni,xtopt,xopt,aopt,matrix(0,0,0),matrix(0,0,0),
                 dmf,matrix(0,0,0),matrix(0,0,0),matrix(0,0,0),itvector,dmfvector,poped.db,
                 rsit=rsit,opt_xt=opt_xt,opt_a=opt_a,opt_x=opt_x)
        }
      }
      
    } else { ##no loops 
      
    }
  }
  # ----------------- RANDOM SEARCH ENDS HERE
  
  # all output should be passed here to poped.db not just these 3
  #poped.db$gxt = xtopt
  #poped.db$gx = xopt
  #poped.db$ga = aopt
  
  poped.db$design$xt <- xtopt
  poped.db$design$x <-xopt
  poped.db$design$a <-aopt
  
  
  # make sure that both the OFV and FIM exists, calculate if not there
  output <- calc_ofv_and_fim(poped.db,
                             ofv=dmf,
                             fim=fmf, 
                             d_switch=d_switch,  
                             bpopdescr=bpopdescr, 
                             ddescr=ddescr,
                             bpop=bpop, 
                             d=d, 
                             docc_full= docc_full, 
                             model_switch=model_switch, 
                             ni=ni, 
                             xt=xt, 
                             x=x, 
                             a=a, 
                             fim.calc.type=fim.calc.type,
                             use_laplace=use_laplace,
                             laplace.fim=laplace.fim, 
                             ...)
  fmf <- output$fim
  dmf <- output$ofv
  
  
  #--------- Write results
  #if((trflag)){
  #  if(footer_flag){
  blockfinal(fn,fmf,dmf,poped.db$design$groupsize,ni,xtopt,xopt,aopt,model_switch,
               bpopdescr,ddescr,poped.db$parameters$docc,poped.db$parameters$sigma,poped.db,
               opt_xt=opt_xt,opt_a=opt_a,opt_x=opt_x,fmf_init=fmf_init,
               dmf_init=dmf_init,compute_inv=compute_inv,
               out_file=out_file,trflag=trflag,
               footer_flag=footer_flag,...)
  #  }
  #}
  
  invisible(list( xtopt= xtopt,xopt=xopt,aopt=aopt,fmf=fmf,dmf=dmf,poped.db =poped.db )) 
}
