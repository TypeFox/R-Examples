#' Summarize your experiment for optimization routines
#' 
#' Create some output to the screen and a text file that summarizes the initial design and the design space
#' you will use to optimize.
#' 
#' @inheritParams RS_opt
#' @inheritParams evaluate.fim
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @inheritParams Dtrace
#' @param fn The file handle to write to.
#' @param e_flag Shuould output be with uncertainty around parameters?
#' 
#' 
#' @family Helper
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_blockexp.R
#' @export
#' @keywords internal
# @importFrom MASS write.matrix

## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

blockexp <- function(fn,poped.db,e_flag=FALSE,
                     opt_xt=poped.db$settings$optsw[2],opt_a=poped.db$settings$optsw[4],opt_x=poped.db$settings$optsw[4],
                     opt_samps=poped.db$settings$optsw[1],opt_inds=poped.db$settings$optsw[5]){
  
  fprintf(fn,'==============================================================================\n')
  fprintf(fn,'Model description : %s \n',poped.db$settings$modtit)
  fprintf(fn,'\n')
  fprintf(fn,'Model Sizes : \n')
  fprintf(fn,'Number of individual model parameters                  g[j]    : Ng    = %g\n',poped.db$parameters$ng)
  fprintf(fn,'Number of population model fixed parameters            bpop[j] : Nbpop = %g\n',poped.db$parameters$nbpop)
  fprintf(fn,'Number of population model random effects parameters   b[j]    : Nb    = %g\n',poped.db$parameters$NumRanEff)
  fprintf(fn,'\n')
  print_params(poped.db$parameters$bpop,"bpop",fn=fn,poped.db=poped.db, 
               head_txt="Typical Population Parameters",e_flag=e_flag)
  fprintf(fn,'\n')
  if((poped.db$parameters$NumRanEff!=0)){
    fprintf(fn,"Between Subject Variability matrix D (variance units) \n")
    d=getfulld(poped.db$parameters$d[,2,drop=F],poped.db$parameters$covd)
    MASS::write.matrix(d,file=fn)
    fprintf(fn,'\n')
    print_params(poped.db$parameters$d,"D",fn=fn,poped.db=poped.db,param_sqrt=TRUE, matrix_elements=T,
                 head_txt="Diagonal Elements of D",e_flag=e_flag)
    fprintf(fn,'\n')
  }
  
  docc_full = getfulld(poped.db$parameters$docc[,2,drop=F],poped.db$parameters$covdocc)
  
  fprintf(fn,'Residual Unexplained Variability matrix SIGMA (variance units) : \n')
  sigma = poped.db$parameters$sigma
  MASS::write.matrix(sigma,file=fn)
  fprintf(fn,'\n')
  
  #sigma_d = diag_matlab(poped.db$parameters$sigma)
  sigma_d <- cbind(c(0,0),diag_matlab(poped.db$parameters$sigma),c(1,1))
  print_params(sigma_d,"SIGMA",fn=fn,poped.db=poped.db,param_sqrt=TRUE, matrix_elements=T,
               head_txt="Diagonal Elements of SIGMA",e_flag=e_flag)
  fprintf(fn,'\n')
  
  fprintf(fn,'==============================================================================\n')
  fprintf(fn,'Experiment description (design and design space)\n')
  fprintf(fn,'\n')
  
  tmp_txt <- "Numer of individuals"
  if(opt_inds) tmp_txt <- paste(tmp_txt,'(min, max)',sep=" ")
  tmp_txt <- paste(tmp_txt,': %g',sep="")
  if(opt_inds) tmp_txt <- paste(tmp_txt,'(%g, %g)',sep=" ")
  tmp_txt <- paste(tmp_txt,'\n',sep="")
  fprintf(fn,tmp_txt,sum(poped.db$design$groupsize),poped.db$design_space$mintotgroupsize,poped.db$design_space$maxtotgroupsize)
  
  fprintf(fn,'Number of groups (individuals with same design): %g\n',poped.db$design$m)
  
  tmp_txt <- "Numer of individuals per group"
  if(opt_inds) tmp_txt <- paste(tmp_txt,'(min, max)',sep=" ")
  tmp_txt <- paste(tmp_txt,':\n',sep="")
  fprintf(fn,tmp_txt)
  
  fprintf(fn," ")
  tmp_txt <- '    Group %g: %g'
  if(opt_inds) tmp_txt <- paste(tmp_txt,'(%g, %g)',sep=" ")
  tmp_txt <- paste(tmp_txt,'\n',sep="")
  fprintf(fn,tmp_txt,1:poped.db$design$m,poped.db$design$groupsize,poped.db$design_space$maxgroupsize, poped.db$design_space$maxgroupsize)
  
  tmp_txt <- "Numer of samples per group"
  if(opt_samps) tmp_txt <- paste(tmp_txt,'(min, max)',sep=" ")
  tmp_txt <- paste(tmp_txt,':\n',sep="")
  fprintf(fn,tmp_txt)
  
  fprintf(fn," ")
  tmp_txt <- '    Group %g: %g'
  if(opt_samps) tmp_txt <- paste(tmp_txt,'(%g, %g)',sep=" ")
  tmp_txt <- paste(tmp_txt,'\n',sep="")
  fprintf(fn,tmp_txt,1:poped.db$design$m,poped.db$design$ni,poped.db$design$minni, poped.db$design$maxni)
  
  fprintf(fn,'Number of discrete experimental variables: %g\n',size(poped.db$design$x,2))
  fprintf(fn,'Number of model covariates: %g\n',size(poped.db$design$a,2))
  
  fprintf(fn,'\n')
  
  print_xt(poped.db$design$xt,poped.db$design$ni,poped.db$design$model_switch,fn,
           head_txt="Initial Sampling Schedule\n")
  fprintf(fn,'\n')
  if(opt_xt){
    print_xt(poped.db$design$xt,poped.db$design$ni,poped.db$design$model_switch,fn,
          head_txt="Minimum allowed sampling values\n",xt_other=poped.db$design_space$minxt)
    fprintf(fn,'\n')
    print_xt(poped.db$design$xt,poped.db$design$ni,poped.db$design$model_switch,fn,
             head_txt="Maximum allowed sampling values\n",xt_other=poped.db$design_space$maxxt)
    fprintf(fn,'\n')
  }  
  
  
  if((size(poped.db$design$x,2)!=0)){
    tmp_txt <- "Discrete Variables"
    if(opt_x) tmp_txt <- paste(tmp_txt,' (possible vales)',sep=" ")
    tmp_txt <- paste(tmp_txt,':\n',sep="")
    fprintf(fn,tmp_txt)
    for(ct1 in 1:poped.db$design$m){
      fprintf(fn,'Group %g: ', ct1)
      for(ct2 in 1:size(poped.db$design$x,2)){
        tmp_txt <- '%g'
        if(opt_x) tmp_txt <- paste(tmp_txt,'(%s)',sep=" ")
        if(ct2<size(poped.db$design$x,2)) tmp_txt <- paste(tmp_txt,' : ',sep="")        
        discrete_val = poped.db$design_space$discrete_x[[ct1,ct2]]  
        fprintf(fn,tmp_txt,poped.db$design$x[ct1,ct2],get_vector_str(discrete_val))
      }
      fprintf(fn,'\n')
    }
    fprintf(fn,'\n')
  }
  
  
  if((size(poped.db$design$a,2)!=0)){   
    tmp_txt <- "Covariates"
    if(opt_a) tmp_txt <- paste(tmp_txt,' (min, max)',sep=" ")
    tmp_txt <- paste(tmp_txt,':\n',sep="")
    fprintf(fn,tmp_txt)
    for(ct1 in 1:poped.db$design$m){
      fprintf(fn,'Group %g: ', ct1)
      for(ct2 in 1:size(poped.db$design$a,2)){
        tmp_txt <- '%g'
        if(opt_a) tmp_txt <- paste(tmp_txt,'(%g, %g)',sep=" ")
        if(ct2<size(poped.db$design$a,2)) tmp_txt <- paste(tmp_txt,' : ',sep="")
        fprintf(fn,tmp_txt,poped.db$design$a[ct1,ct2],poped.db$design_space$mina[ct1,ct2],poped.db$design_space$maxa[ct1,ct2])
      }
      fprintf(fn,'\n')
    }
    fprintf(fn,'\n')
  }
  
  return( ) 
}

print_params <- function (params,name_str, fn, poped.db, param_sqrt=FALSE,head_txt=NULL,matrix_elements=F,e_flag=FALSE) {
  if(is.null(head_txt)) head_txt <- "Parameter Values"
  uncer_txt <- ""
  if(e_flag) uncer_txt <- " (Uncertainty Distribution)"
  sqrt_txt <- ""
  if(param_sqrt) sqrt_txt <- " [sqrt(param)]"
  fprintf(fn,paste(head_txt,sqrt_txt,uncer_txt,":\n",sep=""))
  for(ct in 1:size(params,1)){
    
    par_val <- params[ct,2]
    
    uncer_val <- params[ct,3]
    
    cv_uncer <- sqrt(uncer_val)/par_val*100
    
    cv_str <- ", %CV="
    uncer_str <- ", Var="
    
    par_val_sqrt <- ""
    if(param_sqrt) par_val_sqrt =sqrt(par_val)
    #if(param_sqrt) cv_uncer =cv_uncer/2
    
    if((params[ct,1]==0)){ 
      dist_str <- "Point value" 
      uncer_val <- ""
      cv_uncer <- ""
      cv_str <- ""
      uncer_str <- ""
      
    }
    
    if((params[ct,1]==2)){
      dist_str <- "Uniform" 
      cv_uncer <- ""
      cv_str <- ""
      uncer_str <- ", Max-Min" 
    }
    
    if((params[ct,1]==4)){
      dist_str <- "Log-Normal" 
    }
    if(params[ct,1]==1){
      dist_str <- "Normal" 
    }
    
    if((params[ct,1]==3)){
      dist_str <- "User Defined"
      cv_uncer <- ""
      cv_str <- ""
    }
    
    if((params[ct,1]==5)){
      dist_str <- "Zero-Truncated Normal" 
    }
    
    if(!is.character(uncer_val)) uncer_val <- sprintf("%5.4g",uncer_val)
    if(!is.character(cv_uncer)) cv_uncer <- sprintf("%5.4g",cv_uncer)
    if(!is.character(par_val_sqrt)) par_val_sqrt <- sprintf("[%5.4g] ",par_val_sqrt)
    mat_str <- ""
    if(matrix_elements) mat_str <- sprintf(",%g",ct)

    
    if(e_flag){ 
      fprintf(fn,'%s[%g%s]: %5.4g %s(%s%s%s%s%s)\n', name_str,ct,mat_str,par_val, par_val_sqrt,dist_str, uncer_str, uncer_val, cv_str, cv_uncer)
    } else {
      fprintf(fn,'%s[%g%s]: %5.4g %s\n', name_str,ct,mat_str,par_val, par_val_sqrt,dist_str, uncer_str, uncer_val, cv_str, cv_uncer)
    }
  }
}
