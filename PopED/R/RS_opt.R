#' Optimize the objective function using an adaptive random search algorithm for D-family designs. 
#' 
#' Optimize the objective function using an adaptive random search algorithm.
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
#' @param cfaxt First step factor for sample times 
#' @param opt_xt Should the sample times be optimized?
#' @param opt_a Should the continuous design variables be optimized?
#' @param opt_x Should the discrete design variables be optimized?
#' @param approx_type Approximation method for model, 0=FO, 1=FOCE, 2=FOCEI, 3=FOI.
#' @param iter The number of iterations entered into the \code{\link{blockheader}} function.
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
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_RS_opt.R
#' @export
## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

RS_opt <- function(poped.db,
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
                   #ofv_calc_type=poped.db$settings$ofv_calc_type,
                   #ds_index = poped.db$parameters$ds_index
                   approx_type=poped.db$settings$iApproximationMethod,
                   iter=1,
                   ...){
  # change so that only one output file and one summary file
  # enable turning off of writing output
  # allow for only global structure to be passed. or only input needed for
  # function
  # fix summary file output.
  # fix so that it works for ed as well as d...needs to be more general.
  
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
  
  if((fmf==0) ){
    #     tot.time <- 0
    #     n <- 100
    #     for(i in 1:n){
    #       times <- system.time(mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop,d,poped.db$parameters$sigma,docc_full,poped.db))
    #       tot.time <- tot.time + times
    #     }
    #     avg.time <- tot.time/n
    #     avg.time
    #     
    #     system.time(for(i in 1:100) mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop,d,poped.db$parameters$sigma,docc_full,poped.db))
    #     system.time(mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop,d,poped.db$parameters$sigma,docc_full,poped.db))
    #     
    #     system.time(for(i in 1:100){
    #     evaluate.fim(poped.db,
    #                  bpop.val=bpop,
    #                  d_full=d,
    #                  docc_full=docc_full,
    #                  sigma_full=poped.db$parameters$sigma,
    #                  model_switch=model_switch,
    #                  ni=ni,
    #                  xt=xt,
    #                  x=x,
    #                  a=a,
    #                  groupsize=poped.db$design$groupsize,
    #                  ...)
    #     }    
    #     )
    #     
    fmf <- evaluate.fim(poped.db,
                        bpop.val=bpop,
                        d_full=d,
                        docc_full=docc_full,
                        sigma_full=poped.db$parameters$sigma,
                        model_switch=model_switch,
                        ni=ni,
                        xt=xt,
                        x=x,
                        a=a,
                        groupsize=poped.db$design$groupsize,
                        fim.calc.type=fim.calc.type,
                        ...)
    
    #     returnArgs <-  mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop,d,poped.db$parameters$sigma,docc_full,poped.db) 
    #     fmf <- returnArgs[[1]]
    #     poped.db <- returnArgs[[2]]
  }
  if((dmf==0)){
    dmf=ofv_fim(fmf,poped.db,...)
  }
  fmf_init <- fmf
  dmf_init <- dmf
  
  #==========================================
  # turn warnings off
  #==========================================
  ## warning('off','MATLAB:nearlySingularMatrix');
  
  
  # ------------------ Write summary output file header
  if((trflag)){
    fn=blockheader(poped.db,name='RS',iter=iter,
                     opt_xt=opt_xt,opt_a=opt_a,opt_x=opt_x,
                     opt_inds=F,opt_samps=F,
                     fmf=fmf_init,dmf=dmf_init,
                     bpop=bpopdescr,d=ddescr,docc=poped.db$parameters$docc,sigma=poped.db$parameters$sigma)
    param_vars_init=diag_matlab(inv(fmf))
    returnArgs <-  get_cv(param_vars_init,bpop=bpopdescr,d=ddescr,docc=poped.db$parameters$docc,sigma=poped.db$parameters$sigma,poped.db) 
    params_init <- returnArgs[[1]]
    param_cvs_init <- returnArgs[[2]]
  }
  
  # if((poped.db$settings$bShowGraphs && trflag)){
  #     figure(1)
  #     if((poped.db$settings$Engine$Type==1)){
  #         set(1,'Name','Random Search')
  #     }
  # }
  
  itvector <- c()
  dmfvector <- c()
  
  # ----------------- initialization of optimum variables (and old)
  xtopt=xt
  xopt=x
  aopt=a
  
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
    Dtrace(fn,0,ni,xtopt,xopt,aopt,matrix(0,0,0),matrix(0,0,0),dmf,matrix(0,0,0),matrix(0,0,0),matrix(0,0,0),itvector,dmfvector,poped.db)
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
        Dtrace(fn,it,ni,xtopt,xopt,aopt,matrix(0,0,0),matrix(0,0,0),dmf,matrix(0,0,0),matrix(0,0,0),matrix(0,0,0),itvector,dmfvector,poped.db)
      }
    }
  } else { # } paralell settings
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

        nfmf <- evaluate.fim(poped.db,
                             bpop.val=bpop,
                             d_full=d,
                             docc_full=docc_full,
                             sigma_full=poped.db$parameters$sigma,
                             model_switch=model_switch,
                             ni=ni,
                             xt=xtoptn,
                             x=xoptn,
                             a=aoptn,
                             groupsize=poped.db$design$groupsize,
                             fim.calc.type=fim.calc.type,
                             ...)
        #         returnArgs <- mftot(model_switch,poped.db$design$groupsize,ni,xtoptn,xoptn,aoptn,bpop,d,poped.db$parameters$sigma,docc_full,poped.db) 
        #         nfmf <- returnArgs[[1]]
        #         poped.db <- returnArgs[[2]]
        ndmf=ofv_fim(nfmf,poped.db,...)
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
          Dtrace(fn,it,ni,xtopt,xopt,aopt,matrix(0,0,0),matrix(0,0,0),dmf,matrix(0,0,0),matrix(0,0,0),matrix(0,0,0),itvector,dmfvector,poped.db,
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
  
  #--------- Write results
  if((trflag)){
    blockfinal(fn,fmf,dmf,poped.db$design$groupsize,ni,xtopt,xopt,aopt,model_switch,bpopdescr,ddescr,poped.db$parameters$docc,poped.db$parameters$sigma,poped.db,
                 opt_xt=opt_xt,opt_a=opt_a,opt_x=opt_x,fmf_init=fmf_init,dmf_init=dmf_init,param_cvs_init=param_cvs_init)
    #close(fn)
  }
  
  #==========================================
  # turn warnings on
  #==========================================
  ## warning('on','MATLAB:nearlySingularMatrix');
  
  
  invisible(list( xtopt= xtopt,xopt=xopt,aopt=aopt,fmf=fmf,dmf=dmf,poped.db =poped.db )) 
}
