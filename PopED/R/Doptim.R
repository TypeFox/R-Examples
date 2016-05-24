#' D-family optimization function
#' 
#' Optimize the objective function. There are 4 different optimization
#' algorithms used in this function \enumerate{ \item Adaptive random search.
#' See \code{\link{RS_opt}}. \item Stochastic gradient. \item A Broyden Fletcher
#' Goldfarb Shanno (BFGS) method for nonlinear minimization with box
#' constraints. \item A line search. See \code{\link{a_line_search}}. } The
#' optimization algorithms run in series, taking as input the output from the
#' previous method. The stopping rule used is to test if the line search
#' algorithm fids a better optimum then its inital value. If so, then the chain
#' of algorithms is run again.  If line search is not used then the argument
#' \code{iter_tot} defines the number of times the chain of algorithms is run. 
#' This function takes information from the PopED database supplied as an
#' argument. The PopED database supplies information about the the model,
#' parameters, design and methods to use. Some of the arguments coming from the
#' PopED database can be overwritten; if they are supplied then they are used
#' instead of the arguments from the PopED database.
#' 
#' @inheritParams RS_opt
#' @inheritParams evaluate.fim
#' @inheritParams create.poped.database
#' @param bpopdescr Matrix defining the fixed effects, per row (row number =
#'   parameter_number) we should have: \itemize{ \item column 1 the type of the
#'   distribution for E-family designs (0 = Fixed, 1 = Normal, 2 = Uniform, 3 =
#'   User Defined Distribution, 4 = lognormal and 5 = truncated normal) \item
#'   column 2  defines the mean. \item column 3 defines the variance of the
#'   distribution (or length of uniform distribution). }
#' @param ddescr Matrix defining the diagnonals of the IIV (same logic as for
#'   the \code{bpopdescr}).
#' @param fmf The initial value of the FIM. If set to zero then it is computed.
#' @param dmf The inital OFV. If set to zero then it is computed.
#' @param trflag Should the optimization be output to the screen and to a file?
#' @param ls_step_size Number of grid points in the line search.
#' @param iter_tot Number of iterations to use if line search is not used. Must
#'   be less than \code{iter_max} to be used.
#' @param iter_max If line search is used then the algorithm tests if line
#'   search (always run at the end of the optimization iteration) changes the 
#'   design in any way.  If not, the algorithm stops.  If yes, then a new
#'   iteration is run unless \code{iter_max} iterations have already been run.
#'   
#' @references \enumerate{ \item M. Foracchia, A.C. Hooker, P. Vicini and A.
#'   Ruggeri, "PopED, a software for optimal experimental design in population
#'   kinetics", Computer Methods and Programs in Biomedicine, 74, 2004. \item J.
#'   Nyberg, S. Ueckert, E.A. Stroemberg, S. Hennig, M.O. Karlsson and A.C.
#'   Hooker, "PopED: An extended, parallelized, nonlinear mixed effects models
#'   optimal design tool", Computer Methods and Programs in Biomedicine, 108,
#'   2012. }
#' @family Optimize
#'   
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_Doptim.R
#' @export
#' @keywords internal
#' 
## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

Doptim <- function(poped.db,ni, xt, model_switch, x, a, bpopdescr, 
                   ddescr, maxxt, minxt,maxa,mina,fmf=0,dmf=0,
                   trflag=TRUE,
                   bUseRandomSearch=poped.db$settings$bUseRandomSearch,
                   bUseStochasticGradient=poped.db$settings$bUseStochasticGradient,
                   bUseBFGSMinimizer=poped.db$settings$bUseBFGSMinimizer,
                   bUseLineSearch=poped.db$settings$bUseLineSearch,
                   sgit=poped.db$settings$sgit,
                   ls_step_size=poped.db$settings$ls_step_size,
                   BFGSConvergenceCriteriaMinStep=poped.db$settings$BFGSConvergenceCriteriaMinStep,
                   BFGSProjectedGradientTol=poped.db$settings$BFGSProjectedGradientTol,
                   BFGSTolerancef=poped.db$settings$BFGSTolerancef,
                   BFGSToleranceg=poped.db$settings$BFGSToleranceg,
                   BFGSTolerancex=poped.db$settings$BFGSTolerancex,
                   iter_tot=poped.db$settings$iNumSearchIterationsIfNotLineSearch,
                   iter_max=10,
                   ...){
  
  
  ## update poped.db with options supplied in function
  called_args <- match.call()
  default_args <- formals()
  for(i in names(called_args)[-1]){
    if(length(grep("^poped\\.db\\$",capture.output(default_args[[i]])))==1) {
      #eval(parse(text=paste(capture.output(default_args[[i]]),"<-",called_args[[i]])))
      eval(parse(text=paste(capture.output(default_args[[i]]),"<-",i)))    }
  }
  iter=0 #Search iterations
  test_change=TRUE
  
  trflag = trflag
  
  # ----------------- type of optimization determination
  axt=poped.db$settings$optsw[2]*poped.db$settings$cfaxt*matrix(1,poped.db$design$m,max(poped.db$design_space$maxni))
  aa=poped.db$settings$optsw[4]*poped.db$settings$cfaa*matrix(1,poped.db$design$m,size(poped.db$design$a,2))
  optxt=poped.db$settings$optsw[2]
  optx=poped.db$settings$optsw[3]
  opta=poped.db$settings$optsw[4]
  bfgs_init=matrix(0,0,0)
  # ----------------- initialization of size variables
  m=size(ni,1)
  maxni=size(xt,2)
  
  iMaxSearchIterations = iter_tot
  
  # ----------------- initialization of model parameters
  bpop=bpopdescr[,2,drop=F]
  d=getfulld(ddescr[,2,drop=F],poped.db$parameters$covd)
  docc_full = getfulld(poped.db$parameters$docc[,2,drop=F],poped.db$parameters$covdocc)
  
  Engine = list(Type=1,Version=version$version.string)
  
  if((dmf==0) ){#Only first time
    returnArgs <-  mftot(model_switch,poped.db$design$groupsize,ni,xt,x,a,bpop,d,poped.db$parameters$sigma,docc_full,poped.db) 
    fmf <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    dmf=ofv_fim(fmf,poped.db)
    fmf_init <- fmf
    dmf_init <- dmf
  }
  
  if((optxt || optx || opta)){
    # ------------------ Write output file header
    #fn=blockheader(FALSE,iter,poped.db)
    #if((trflag)){
    fn=blockheader(poped.db,name='D_cont_opt',iter=1,
                     opt_xt=optxt,opt_a=opta,opt_x=optx,
                     opt_inds=F,opt_samps=F,
                     fmf=fmf_init,dmf=dmf_init,
                     bpop=bpopdescr,d=ddescr,docc=poped.db$parameters$docc,sigma=poped.db$parameters$sigma,
                     trflag=trflag,
                     ...)
    #       param_vars_init=diag_matlab(inv(fmf))
    #       returnArgs <-  get_cv(param_vars_init,bpop=bpopdescr,d=ddescr,docc=poped.db$parameters$docc,sigma=poped.db$parameters$sigma,poped.db) 
    #       params_init <- returnArgs[[1]]
    #       param_cvs_init <- returnArgs[[2]]
    #}
    while(test_change==TRUE && iMaxSearchIterations>0 && iter < iter_max){
      iter=iter+1
      
      
      
      
      if(((bUseRandomSearch || bUseStochasticGradient) && poped.db$settings$bShowGraphs)){
        #figure(1)
        #if((poped.db$settings$Engine$Type==1)){
        #    set(1,'Name','Random Search and Stochastic Gradient')
        #}
      }
      
      itvector <- c()
      dmfvector <- c()
      
      # ----------------- initialization of optimum variables (and old)
      xtopt=xt
      xopt=x
      aopt=a
      xtopto=xt
      xopto=x
      aopto=a
      lgxto=zeros(m*maxni)
      lgao=zeros(m*size(poped.db$design$a,2))
      
      # ----------------- RS variables initialization
      dxt=(maxxt-minxt)/poped.db$settings$rslxt
      da=(maxa-mina)/poped.db$settings$rsla
      xtoptn=xtopt
      xoptn=xopt
      aoptn=aopt
      nullit=1
      ff=1
      
      # ----------------- SG variables initialization
      kitxt=1
      kita=1
      mnormxt=zeros(m,maxni)
      mnorma=zeros(m,size(poped.db$design$a,2))
      normgxt=zeros(m,maxni)
      normga=zeros(m,size(poped.db$design$a,2))
      
      # ----------------- initialization of trace support variables
      odmf=0
      inversionxt=FALSE
      inversiona=FALSE
      
      if((bUseRandomSearch) ){#If we want to perform random search
        # ----------------- RANDOM SEARCH BEGINS HERE
        #save for graphical output
        if((poped.db$settings$parallel$bParallelLS == 0)){
          tic()              
        }
        itvector[1] = 0
        dmfvector[1] = dmf
        
        if((trflag)){
          Dtrace(fn,0,ni,xtopt,xopt,aopt,matrix(0,0,0),matrix(0,0,0),dmf,matrix(0,0,0),matrix(0,0,0),matrix(0,0,0),itvector,dmfvector,poped.db)
        }
        
        if((poped.db$settings$parallel$bParallelRS)){
          
          # Generate the input designs
          
          designsin = cell(1,0)
          for(it in 1:poped.db$settings$rsit){
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
          
          stop("Parallel execution not yet implemented in R")
          designsout = designsin
          #designsout = execute_parallel(designsin,poped.db)
          #Store the optimal design
          for(it in 1:poped.db$settings$rsit){
            if((designsout[[it]]$ofv>dmf) ){##ok<USENS>
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
            
            if((trflag && (rem(it,poped.db$settings$rsit_output)==0 || it==poped.db$settings$rsit))){
              itvector[ceiling(it/poped.db$settings$rsit_output)+1]=it
              dmfvector[ceiling(it/poped.db$settings$rsit_output)+1]=dmf
              Dtrace(fn,it,ni,xtopt,xopt,aopt,matrix(0,0,0),matrix(0,0,0),dmf,matrix(0,0,0),matrix(0,0,0),matrix(0,0,0),itvector,dmfvector,poped.db)
            }
          }
        } else {
          for(it in 1:poped.db$settings$rsit){
            
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
            returnArgs <- mftot(model_switch,poped.db$design$groupsize,ni,xtoptn,xoptn,aoptn,bpop,d,poped.db$parameters$sigma,docc_full,poped.db) 
            nfmf <- returnArgs[[1]]
            poped.db <- returnArgs[[2]]
            ndmf=ofv_fim(nfmf,poped.db)
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
            if((nullit==poped.db$settings$maxrsnullit)){
              ff=ff+1
              nullit=1
            }
            
            if((!isempty(poped.db$settings$strIterationFileName))){
              write_iterationfile('Random Search',it,xtopt,aopt,xopt,ni,poped.db$design$groupsize,fmf,dmf,poped.db)
            }
            
            if((trflag && (rem(it,poped.db$settings$rsit_output)==0 || it==poped.db$settings$rsit))){
              itvector[ceiling(it/poped.db$settings$rsit_output)+1]=it
              dmfvector[ceiling(it/poped.db$settings$rsit_output)+1]=dmf
              Dtrace(fn,it,ni,xtopt,xopt,aopt,matrix(0,0,0),matrix(0,0,0),dmf,matrix(0,0,0),matrix(0,0,0),matrix(0,0,0),itvector,dmfvector,poped.db)
            }
          }  
        }
        
        if((poped.db$settings$parallel$bParallelLS == 0)){
          timeLS = toc(echo=FALSE)
          fprintf('Run time for random search: %g seconds\n\n',timeLS)
          #if(trflag) fprintf(fn,'Elapsed time for Serial Random search with %g iterations: %g seconds\n',it,timeLS)
        }
        # ----------------- RANDOM SEARCH ENDS HERE 
      }
      # ----------------- initialization of best optimum variables
      
      bestxt=xtopt
      bestx=xopt
      besta=aopt
      
      diff=poped.db$settings$convergence_eps+1
      
      if((bUseBFGSMinimizer )){
        if((Engine$Type==2)){
          stop(sprintf('BFGS optimization can not be used with Freemat in this version!'))
        }
        
        f_name <- 'calc_ofv_and_grad' 
        f_options <- list(x,optxt, opta, model_switch,aa,axt,poped.db$design$groupsize,ni,xtopt,xopt,aopt,bpop,d,poped.db$parameters$sigma,docc_full,poped.db)
        x_k=matrix(0,0,0)
        lb=matrix(0,0,0)
        ub=matrix(0,0,0)
        options=list('factr'=BFGSConvergenceCriteriaMinStep,'pgtol'=BFGSProjectedGradientTol,'ftol'=BFGSTolerancef,'gtol'=BFGSToleranceg,'xtol'=BFGSTolerancex)
        if((optxt==TRUE)){
          index=t(1:numel(xtopt))
          if(poped.db$design_space$bUseGrouped_xt){
            returnArgs <- unique(poped.db$design_space$G_xt) 
            temp <- returnArgs[[1]]
            index <- returnArgs[[2]]
            temp2 <- returnArgs[[3]]
          }
          index=index[minxt!=maxxt]
          x_k=t(t(xtopt[index]))
          lb=t(t(minxt[index]))
          ub=t(t(maxxt[index]))
        }
        if((opta==TRUE)){
          index=t(1:numel(aopt))
          if(poped.db$design_space$bUseGrouped_a){
            returnArgs <- unique(poped.db$design_space$G_a) 
            temp1 <- returnArgs[[1]]
            index <- returnArgs[[2]]
            temp2 <- returnArgs[[3]]
          }
          index=index[mina!=maxa]
          x_k=t(t(c(x_k,aopt[index])))
          lb=t(t(c(lb,mina[index])))
          ub=t(t(c(ub,maxa[index])))
          
          #x_k(end+index)=aopt[index]
          #lb(end+index)=mina[index]
          #ub(end+index)=maxa[index]
        }
        if((any(x_k<lb))){
          x_k[x_k<lb]=lb[x_k<lb]
        }
        if((isempty(bfgs_init) || any(x_k!=bfgs_init))){
          bfgs_init=x_k
          fprintf('Starting BGFS minimization with OFV of %g \n', dmf)
          returnArgs <- bfgsb_min(f_name,f_options, x_k,lb,ub,options) 
          x_opt  <- returnArgs[[1]]
          f_k <- returnArgs[[2]]
          B_k <- returnArgs[[3]]
          
          if(optxt){
            notfixed=minxt!=maxxt
            if(poped.db$design_space$bUseGrouped_xt){
              xtopt[notfixed]=x_opt[poped.db$design_space$G_xt[notfixed]]
              x_opt <- x_opt[-(1:numel(unique(poped.db$design_space$G_xt[notfixed])))]
            } else {
              xtopt[notfixed]=x_opt[1:numel(xtopt[notfixed])]
              x_opt <- x_opt[-(1:numel(xtopt[notfixed]))]
            }
          }
          
          if(opta){
            notfixed=mina!=maxa
            
            if(poped.db$design_space$bUseGrouped_a){
              aopt[notfixed]=x_opt(poped.db$design_space$G_a[notfixed])
            } else {
              aopt[notfixed]=x_opt
            }
          }
          returnArgs <-  mftot(model_switch,poped.db$design$groupsize,ni,xtopt,xopt,aopt,bpop,d,poped.db$parameters$sigma,docc_full,poped.db) 
          nfmf <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
          ndmf=ofv_fim(nfmf,poped.db)
          
          fprintf('BFGS minimization finished. New OFV: %g \n', ndmf)
          if((ndmf>dmf)){
            dmf=ndmf
            fmf=nfmf
            bestxt=xtopt
            bestx=xopt
            besta=aopt
          }
        }
      }
      
      if((bUseStochasticGradient)){
        tic()
        # ----------------- SG AUTO-FOCUS BEGINS HERE
        if((optxt==TRUE)){
          returnArgs <-  gradofv_xt(model_switch,axt,poped.db$design$groupsize,ni,xtopt,xopt,aopt,bpop,d,poped.db$parameters$sigma,docc_full,poped.db) 
          gradxt <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
          normgxt=sign(gradxt)*(maxxt-minxt)
          xtopto=xtopt
          returnArgs <-  calc_autofocus(m,ni,dmf,xtopt,xtopto,maxxt,minxt,gradxt,normgxt,axt,model_switch,poped.db$design$groupsize,xtopt,xopt,aopt,ni,bpop,d,poped.db$parameters$sigma,docc_full,poped.db) 
          axt <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
        }
        
        if((optx==TRUE)){
          # No gradient optimization for discrete values
          xopto=xopt
        }
        
        if((opta==TRUE)){
          returnArgs <- gradofv_a(model_switch,aa,poped.db$design$groupsize,ni,xtopt,xopt,aopt,bpop,d,poped.db$parameters$sigma,docc_full,poped.db) 
          grada <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
          normga=grada/abs(grada)*(maxa-mina)
          aopto=aopt
          na_vector = matrix(1,1,m)*size(poped.db$design$a,2)
          returnArgs <- calc_autofocus(m,na_vector,dmf,aopt,aopto,maxa,mina,grada,normga,aa,model_switch,poped.db$design$groupsize,xtopt,xopt,aopt,ni,bpop,d,poped.db$parameters$sigma,docc_full,poped.db) 
          aa <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
        }
        
        xtopt=xtopto
        xopt=xopto
        aopt=aopto
        
        if((matrix_all(axt==0) && matrix_all(aa==0))){
          diff=0
        }
        
        # ----------------- SG AUTO-FOCUS ENDS HERE
        #
        # ----------------- STOCHASTIC GRADIENT BEGINS HERE
        it=1
        dmfvector <- c()
        itvector <- c()
        while((it<=sgit) && (abs(diff)>poped.db$settings$convergence_eps)){
          if((optxt==TRUE)){
            returnArgs <- gradofv_xt(model_switch,axt,poped.db$design$groupsize,ni,xtopto,xopto,aopto,bpop,d,poped.db$parameters$sigma,docc_full,poped.db) 
            gradxt <- returnArgs[[1]]
            poped.db <- returnArgs[[2]]
            returnArgs <- sg_search(gradxt,mnormxt,axt,maxxt,minxt,xtopto,lgxto,kitxt,it,m,maxni) 
            lgxto <- returnArgs[[1]]
            kitxt <- returnArgs[[2]]
            inversionxt <- returnArgs[[3]]
            xtopt <- returnArgs[[4]]
          }
          if((opta==TRUE)){
            returnArgs <- gradofv_a(model_switch,aa,poped.db$design$groupsize,ni,xtopto,xopto,aopto,bpop,d,poped.db$parameters$sigma,docc_full,poped.db) 
            grada <- returnArgs[[1]]
            poped.db <- returnArgs[[2]]
            returnArgs <- sg_search(grada,mnorma,aa,maxa,mina,aopto,lgao,kita,it,m,size(poped.db$design$a,2)) 
            lgao <- returnArgs[[1]]
            kita <- returnArgs[[2]]
            inversiona <- returnArgs[[3]]
            aopt <- returnArgs[[4]]
          }
          xtopto=xtopt
          xopto=xopt
          aopto=aopt
          
          
          if((sum(sum(isnan(xtopt)))>0 || sum(sum(isnan(xopt)))>0 || sum(sum(isnan(aopt)))>0)){
            break
          }
          
          
          returnArgs <-  mftot(model_switch,poped.db$design$groupsize,ni,xtopto,xopto,aopto,bpop,d,poped.db$parameters$sigma,docc_full,poped.db) 
          nfmf <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
          ndmf=ofv_fim(nfmf,poped.db)
          
          if((ndmf>dmf)){
            dmf=ndmf
            fmf=nfmf
            bestxt=xtopt
            bestx=xopt
            besta=aopt
          }
          
          if((!isempty(poped.db$settings$strIterationFileName))){
            write_iterationfile('Stochastic Gradient',it,bestxt,besta,bestx,ni,poped.db$design$groupsize,fmf,dmf,poped.db)
          }
          
          if((poped.db$settings$convergence_eps!=0 || trflag==TRUE)){
            if((it>1)){
              if((odmf!=0)){
                diff=abs((ndmf-odmf))/abs(odmf)
              }
            }
            odmf=ndmf
          }
          if((trflag==TRUE && (rem(it,poped.db$settings$sgit_output)==0 || abs(diff)<poped.db$settings$convergence_eps || it==sgit))){
            itvector[ceiling(it/poped.db$settings$sgit_output)]=it
            dmfvector[ceiling(it/poped.db$settings$sgit_output)]=dmf
            ga_tmp=0
            if(opta) ga_tmp = grada/dmf
            gxt_tmp=0
            if(optxt) gxt_tmp = gradxt/dmf
            Dtrace(fn,poped.db$settings$rsit+it,ni,bestxt,bestx,besta,gxt_tmp,ga_tmp,dmf,diff,inversionxt,inversiona,itvector,dmfvector,poped.db)
            #Dtrace(fn,poped.db$settings$rsit+it,ni,bestxt,bestx,besta,normgxt,normga,dmf,diff,inversionxt,inversiona,itvector,dmfvector,poped.db)
            inversionxt=FALSE
            inversiona=FALSE
          }
          it=it+1
        }
        # ----------------- STOCHASTIC GRADIENT ENDS HERE
        sg_time=toc(echo=FALSE)
        fprintf('Stochastic gradient run time: %g seconds\n\n',sg_time)
        #if(trflag) fprintf(fn,'Elapsed time for stochastic gradient search with %g iterations: %g seconds\n',sgit,sg_time)
        
        
      }
      xtopt=bestxt
      xopt=bestx
      aopt=besta
      xt=xtopt
      x=xopt
      a=aopt
      
      #poped.db$gxt = xtopt
      #poped.db$gx = xopt
      #poped.db$design$a = aopt
      
      poped.db$design$xt <- xtopt
      poped.db$design$x <-xopt
      poped.db$design$a <-aopt
      
      if((bUseLineSearch)){
        #------------------------------- LINE SEARCH optimization START HERE
        strLineSearchFile=sprintf('%s_LS_%g%s',poped.db$settings$strOutputFileName,iter,poped.db$settings$strOutputFileExtension)
        strLineSearchFile = fullfile(poped.db$settings$strOutputFilePath,strLineSearchFile)
        #returnArgs <- a_line_search(strLineSearchFile,FALSE,0,fmf,dmf,poped.db) 
        returnArgs <- a_line_search(poped.db,fn,FALSE,0,fmf,dmf) 
        fmf <- returnArgs[[1]]
        dmf <- returnArgs[[2]]
        test_change <- returnArgs[[3]]
        xt <- returnArgs[[4]]
        x <- returnArgs[[5]]
        a <- returnArgs[[6]]
        poped.db <- returnArgs[[7]]
        #------------------------------- LINE SEARCH optimization ENDS HERE
      } else {
        iMaxSearchIterations = iMaxSearchIterations-1
      }
      
      
    }
    #--------- Write results
    #     if((trflag)){
    blockfinal(fn,fmf,dmf,poped.db$design$groupsize,ni,xt,x,a,model_switch,
                 bpopdescr,ddescr,poped.db$parameters$docc,poped.db$parameters$sigma,poped.db,
                 opt_xt=optxt,opt_a=opta,opt_x=optx,fmf_init=fmf_init,
                 dmf_init=dmf_init,trflag=trflag,...)
    #}
    #blockfinal(fn,fmf,dmf,poped.db$design$groupsize,ni,xt,x,a,bpopdescr,ddescr,poped.db$parameters$docc,poped.db$parameters$sigma,m,poped.db)
    
    #close(fn)
    
  } else {
    fprintf('No Continuous Optimization Performed \n')
  } # matches -- if (optxt|optx|opta)
  return(list( xt= xt,x=x,a=a,fmf=fmf,dmf=dmf,poped.db =poped.db )) 
}

#' Compute the autofocus portion of the stochastic gradient routine
#' 
#' Compute the autofocus portion of the stochastic gradient routine
#'    
#' @return A list containing:
#' \item{navar}{The autofocus parameter.}
#' \item{poped.db}{PopED database.}
#' 
#' 
#' @inheritParams RS_opt
#' @inheritParams evaluate.fim
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @param ni_var The ni_var.
#' @param varopt The varopt.
#' @param varopto The varopto.
#' @param maxvar The maxvar.
#' @param minvar The minvar.
#' @param gradvar The gradvar.
#' @param normgvar The normgvar.
#' @param avar The avar.
#' @param xtopt The optimal sampling times matrix.
#' @param xopt The optimal discrete design variables matrix.
#' @param aopt The optimal continuous design variables matrix.
#' @family Optimize
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_calc_autofocus.R
#' @export
#' @keywords internal
calc_autofocus <- function(m,ni_var,dmf,varopt,varopto,maxvar,minvar,gradvar,normgvar,
                           avar,model_switch,groupsize,xtopt,xopt,aopt,ni,bpop,d,sigma,docc,poped.db){
  navar = avar
  for(i in 1:m){
    for(ct1 in 1:ni_var[i]){
      if((varopt[i,ct1]==maxvar[i,ct1] && gradvar[i,ct1]>0)){
        avar[i,ct1]=0
      }
      if((varopt[i,ct1]==minvar[i,ct1] && gradvar[i,ct1]<0)){
        avar[i,ct1]=0
      }
      if((avar[i,ct1]!=0)){
        varopt=varopto
        tavar=avar[i,ct1]
        varopt[i,ct1]=varopto[i,ct1]+tavar*normgvar[i,ct1]
        varopt[i,ct1]=limit_value(varopt[i,ct1],maxvar[i,ct1],minvar[i,ct1])
        returnArgs <-  mftot(model_switch,groupsize,ni,xtopt,xopt,aopt,bpop,d,sigma,docc,poped.db) 
        mf_tmp <- returnArgs[[1]]
        poped.db <- returnArgs[[2]]
        ndmf=ofv_fim(mf_tmp,poped.db)
        ct2=1
        while(ct2<=10 && ndmf<dmf){
          tavar=tavar/2
          varopt[i,ct1]=varopto[i,ct1]+tavar*normgvar[i,ct1]
          varopt[i,ct1]=limit_value(varopt[i,ct1],maxvar[i,ct1],minvar[i,ct1])
          returnArgs <- mftot(model_switch,groupsize,ni,xtopt,xopt,aopt,bpop,d,sigma,docc,poped.db) 
          mf_tmp <- returnArgs[[1]]
          poped.db <- returnArgs[[2]]
          ndmf=ofv_fim(mf_tmp,poped.db)
          ct2=ct2+1
        }
        if((ndmf<dmf)){
          navar[i,ct1]=0
        } else {
          navar[i,ct1]=tavar
        }
      }
    }
  }
  return(list( navar= navar,poped.db=poped.db)) 
}

sg_search <- function(graddetvar,mnormvar,avar,maxvar,minvar,varopto,lgvaro,oldkitvar,it,m,numvar){
  mnormvar=((mnormvar*(it-1))+(graddetvar^2))/it
  normgvar=(avar/oldkitvar)*graddetvar/sqrt(mnormvar)*(maxvar-minvar)
  varopt=varopto+normgvar
  varopt=varopt-((varopt>maxvar)*(varopt-maxvar))
  varopt=varopt+((varopt<minvar)*(minvar-varopt))
  lgvar=varopt-varopto
  lgvar=reshape_matlab(lgvar,m*numvar,1)
  if(any(t(lgvar)%*%lgvaro<0)){
    kitvar=oldkitvar+1
    inversionvar=TRUE
  } else {
    inversionvar=FALSE
    kitvar=oldkitvar
  }
  lgvaro=lgvar
  return(list( lgvaro= lgvaro, kitvar= kitvar, inversionvar= inversionvar, varopt= varopt)) 
}

# sign <- function(x){
# # s = zeros(size(x))
# s = x/abs(x)*(x!=0)
# return( s) 
# }


#' Compute an objective function and gradient 
#' 
#' Compute an objective function and gradient with respect to the optimization parameters.
#' This function can be passed to the Broyden Fletcher Goldfarb Shanno (BFGS) 
#' method for nonlinear minimization with box constraints implemented in \code{\link{bfgsb_min}}.
#'   
#' @return A list containing:
#' \item{f}{The objective function.}
#' \item{g}{The gradient.}
#' 
#' 
#' @inheritParams RS_opt
#' @inheritParams evaluate.fim
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @inheritParams calc_autofocus
#' @param aa The aa value
#' @param axt the axt value
#' @param xtopto the xtopto value
#' @param xopto the xopto value
#' @param optxt If sampling times are optimized
#' @param opta If continuous design variables are optimized
#' @param aopto the aopto value
#' @param only_fim Should the gradient be calculated?
#' @family Optimize
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_calc_ofv_and_grad.R
#' @export
#' @keywords internal
calc_ofv_and_grad <- function(x,optxt,opta, model_switch,aa,axt,groupsize,ni,xtopto,
                              xopto,aopto,bpop,d,sigma,docc_full,poped.db,only_fim=FALSE){
  if(optxt){
    notfixed <- poped.db$design_space$minxt!=poped.db$design_space$maxxt
    if(poped.db$design_space$bUseGrouped_xt){
      xtopto[notfixed]=x[poped.db$design_space$G_xt[notfixed]]
      ##x[1:numel(unique(poped.db$design_space$G_xt[notfixed]))]=matrix(0,0,0)
      x=x[-c(1:numel(unique(poped.db$design_space$G_xt[notfixed])))]
    } else {
      xtopto[notfixed]=x[1:numel(xtopto[notfixed])]
      x=x[-c(1:numel(xtopto[notfixed]))]
      ##x[1:numel(xtopto[notfixed])]=matrix(0,0,0)
    }
  }
  if(opta){
    notfixed <- poped.db$design_space$mina!=poped.db$design_space$maxa
    if(poped.db$design_space$bUseGrouped_a){
      aopto[notfixed]=x[poped.db$design_space$G_a[notfixed]]
    } else {
      aopto[notfixed]=x
    }
  }
  returnArgs <-  mftot(model_switch,groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc_full,poped.db) 
  nfmf <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  f=-ofv_fim(nfmf,poped.db)
  gradxt=matrix(0,0,0)
  grada=matrix(0,0,0)
  if((only_fim)){
    g=0
    return(list( f= f,g=g))
  }
  if((optxt==TRUE)){
    notfixed=poped.db$design_space$minxt!=poped.db$design_space$maxxt
    returnArgs <- gradofv_xt(model_switch,axt,groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc_full,poped.db) 
    gradxt <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    gradxt=gradxt[notfixed]
    if(poped.db$design_space$bUseGrouped_xt){
      returnArgs <- unique(poped.db$design_space$G_xt) 
      temp1 <- returnArgs[[1]]
      index <- returnArgs[[2]]
      temp2 <- returnArgs[[3]]
      gradxt=gradxt[index]
    }
  }
  if((opta==TRUE)){
    notfixed=poped.db$design_space$mina!=poped.db$design_space$maxa
    returnArgs <- gradofv_a(model_switch,aa,groupsize,ni,xtopto,xopto,aopto,bpop,d,sigma,docc_full,poped.db) 
    grada <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    grada=grada[notfixed]
    if(poped.db$design_space$bUseGrouped_a){
      returnArgs <- unique(poped.db$design_space$G_a) 
      temp1 <- returnArgs[[1]]
      index <- returnArgs[[2]]
      temp2 <- returnArgs[[3]]
      grada=grada[index]
    }
  }
  g=-matrix(c(gradxt,grada),ncol=1,byrow=T)
  return(list( f= f,g=g)) 
}

