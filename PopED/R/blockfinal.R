#' Result function for optimization routines
#' 
#' Create some output to the screen and a text file that summarizes the problem you solved.
#' 
#' @inheritParams RS_opt
#' @inheritParams evaluate.fim
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @inheritParams blockexp
#' @inheritParams blockheader
#' @inheritParams RS_opt_gen
#' @param fmf_init Initial FIM.
#' @param dmf_init Initial OFV.
#' @param param_cvs_init The inital design parameter RSE values.
#' 
#' @family Helper
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_blockfinal.R
#' @export
#' @keywords internal
#' 
#' 
# @importFrom MASS write.matrix
## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

blockfinal <- function(fn,fmf,dmf,groupsize,ni,xt,x,a,model_switch,bpop,d,docc,sigma,poped.db,
                         opt_xt=poped.db$settings$optsw[2],opt_a=poped.db$settings$optsw[4],opt_x=poped.db$settings$optsw[3],
                         fmf_init=NULL,dmf_init=NULL,param_cvs_init=NULL,
                         compute_inv=TRUE,out_file=NULL,trflag=TRUE,footer_flag=TRUE,...){
  
  if(!trflag) return(invisible() ) 
  if(footer_flag){
    if(is.null(fmf)) compute_inv <- FALSE
    if(!is.matrix(fmf)) compute_inv <- FALSE
    
    
    
    fprintf(fn,'===============================================================================\nFINAL RESULTS\n')
    if(fn!="") fprintf('===============================================================================\nFINAL RESULTS\n')
    
    time_value <- NULL
    if(exists(".poped_total_time", envir=.PopedNamespaceEnv)) time_value = toc(echo=FALSE,name=".poped_total_time")
    if((opt_xt==TRUE)){
      print_xt(xt,ni,model_switch,fn,head_txt="Optimized Sampling Schedule\n")
      if(fn!="") print_xt(xt,ni,model_switch,head_txt="\nOptimized Sampling Schedule\n")
    }
    if((opt_x==TRUE)){
      #     fprintf(fn,'x :\n')
      #     fprintf(fn,'%g \n',x)
      #     cat("Optimized x values:\n")
      #     print(x)
      tmp_txt <- "\nOptimized Discrete Variables"
      tmp_txt <- paste(tmp_txt,':\n',sep="")
      fprintf(fn,tmp_txt)
      if(fn!="") fprintf(tmp_txt)
      for(ct1 in 1:poped.db$design$m){
        fprintf(fn,'Group %g: ', ct1)
        if(fn!="") fprintf('Group %g: ', ct1)
        for(ct2 in 1:size(poped.db$design$x,2)){
          tmp_txt <- '%g'
          if(ct2<size(poped.db$design$x,2)) tmp_txt <- paste(tmp_txt,' : ',sep="")        
          fprintf(fn,tmp_txt,x[ct1,ct2])
          if(fn!="") fprintf(tmp_txt,x[ct1,ct2])
        }
        fprintf(fn,'\n')
        if(fn!="") fprintf('\n')
      }
      #     fprintf(fn,'\n')
      #     fprintf('\n')
    }
    if((opt_a==TRUE)){
      tmp_txt <- "\nOptimized Covariates"
      tmp_txt <- paste(tmp_txt,':\n',sep="")
      fprintf(fn,tmp_txt)
      if(fn!="") fprintf(tmp_txt)
      for(ct1 in 1:poped.db$design$m){
        fprintf(fn,'Group %g: ', ct1)
        if(fn!="") fprintf('Group %g: ', ct1)
        for(ct2 in 1:size(poped.db$design$a,2)){
          tmp_txt <- '%g'
          if(ct2<size(poped.db$design$a,2)) tmp_txt <- paste(tmp_txt,' : ',sep="")
          fprintf(fn,tmp_txt,a[ct1,ct2])
          if(fn!="") fprintf(tmp_txt,a[ct1,ct2])
        }
        fprintf(fn,'\n')
        if(fn!="") fprintf('\n')
      }
      #     fprintf(fn,'\n')
      #     fprintf('\n')
      #     
      #     fprintf(fn,'a :\n')
      #     fprintf(fn,'%g \n',a)
      #     cat("Optimized a values:\n")
      #     print(a)
    }
    if((poped.db$settings$d_switch==TRUE && (fn!="" || trflag>1))){
      fprintf(fn,'\n FIM: \n')
      #write_matrix(fn,fmf)
      MASS::write.matrix(fmf,file=fn)
      fprintf(fn,'\n\nInverse(FIM):\n')
      #write_matrix(fn,inv(fmf))
      if(compute_inv) MASS::write.matrix(inv(fmf),file=fn)
    }
    fprintf(fn,'\nOFV = %g\n',dmf)
    if(fn!="") fprintf('\nOFV = %g\n',dmf)
    
    if(compute_inv){
      param_vars=diag_matlab(inv(fmf))
      returnArgs <-  get_cv(param_vars,bpop,d,docc,sigma,poped.db) 
      params <- returnArgs[[1]]
      param_cvs <- returnArgs[[2]]
    }
    
    output <- get_unfixed_params(poped.db)
    npar <- length(output$all)
    
    if(fn!="" || trflag>1) fprintf(fn,'\nEfficiency criterion [usually defined as OFV^(1/npar)]  = %g\n',
            ofv_criterion(dmf,npar,poped.db))
    
    fprintf(fn,'\nEfficiency [typically: (OFV_final/OFV_initial)^(1/npar)]: %g\n',
            ofv_criterion(dmf,npar,poped.db)/ofv_criterion(dmf_init,npar,poped.db))
    if(fn!=""){
      fprintf('\nEfficiency [typically: (OFV_final/OFV_initial)^(1/npar)]: %g\n',
              ofv_criterion(dmf,npar,poped.db)/ofv_criterion(dmf_init,npar,poped.db))
    }
    #fprintf(fn,'\nEfficiency criterion: det(FIM)^(1/npar) = %g\n',dmf^(1/length(params)))
    #fprintf(fn,'\nEfficiency (final_design/initial_design): %g\n',(dmf^(1/length(params)))/(dmf_init^(1/length(params))))
    #if(fn!="") fprintf('\nEfficiency (final_design/initial_design): %g\n',(dmf^(1/length(params)))/(dmf_init^(1/length(params))))
    
    if(is.null(param_cvs_init) && !is.null(fmf_init) && is.matrix(fmf_init) && compute_inv){
      param_vars_init=diag_matlab(inv(fmf_init))
      returnArgs <-  get_cv(param_vars_init,bpop,d,docc,sigma,poped.db) 
      params_init <- returnArgs[[1]]
      param_cvs_init <- returnArgs[[2]]
    }
    
    if(compute_inv){
      parnam <- get_parnam(poped.db)
      fprintf(fn,'\nExpected parameter \nrelative standard error (%sRSE):\n','%')
      if(fn!="") fprintf('\nExpected parameter \nrelative standard error (%sRSE):\n','%')
      df <- data.frame("Parameter"=parnam,"Values"=params, #"Variance"=param_vars, 
                       "RSE_0"=t(param_cvs_init*100),"RSE"=t(param_cvs*100))
      print(df,digits=3, print.gap=3,row.names=F)
      if(fn!="") capture.output(print(df,digits=3, print.gap=3,row.names=F),file=fn)
    }
    
    #if(!is.null(time_value)){
    fprintf(fn,'\nTotal running time: %g seconds\n',time_value)
    if(fn!="") fprintf('\nTotal running time: %g seconds\n',time_value)
    #}
  } # end footer_flag
  if(!any(class(out_file)=="file") && (fn != '')) close(fn)
  
  return(invisible() ) 
}

print_xt <- function (xtopt, ni, model_switch,fn="",head_txt="Optimized xt values:\n",xt_other=NULL) {
  cat(head_txt,file=fn)
  for(j in 1:size(xtopt,1)){
    xtopt_i = xtopt[j,1:ni[j]]
    model_switch_i = model_switch[j,1:ni[j]]
    if(!is.null(xt_other)) xt_other_i = xt_other[j,1:ni[j]]
    for(i in unique(as.vector(model_switch_i))){
      xtopt_i_sort = sort(xtopt_i[model_switch_i==i])
      if(!is.null(xt_other)) xt_other_i_sort = xt_other_i[order(xtopt_i[model_switch_i==i])]
      if(length(unique(as.vector(model_switch_i)))>1) cat(sprintf("Model %g : ", i),file=fn)
      if(size(xtopt,1)>1) cat(sprintf("Group %g : ", j),file=fn)
      if(!is.null(xt_other)) {
        cat(sprintf("%6.4g", xt_other_i_sort),file=fn)
      } else {
        cat(sprintf("%6.4g", xtopt_i_sort),file=fn)
      }
      cat("\n",file=fn)
    }
  }
  invisible()
}

get_parnam <- function (poped.db) {
  nbpop = length(poped.db$parameters$notfixed_bpop)
  nd = length(poped.db$parameters$notfixed_d)
  ncovd = length(poped.db$parameters$notfixed_covd)
  ndocc = length(poped.db$parameters$notfixed_docc)
  ncovdocc = length(poped.db$parameters$notfixed_covdocc)
  nsigma = length(poped.db$parameters$notfixed_sigma)
  ncovsigma = length(poped.db$parameters$notfixed_covsigma)
  
  not_fixed <- list("bpop"=poped.db$parameters$notfixed_bpop,
                    "D"=poped.db$parameters$notfixed_d,
                    "D_cov"=poped.db$parameters$notfixed_covd,
                    "D.occ"=poped.db$parameters$notfixed_docc,
                    "D.occ_cov"=poped.db$parameters$notfixed_covdocc,
                    "SIGMA"=poped.db$parameters$notfixed_sigma,
                    "SIGMA_cov"=poped.db$parameters$notfixed_covsigma)
  parnam <- c()
  for(i in 1:size(not_fixed,2)){
    if(length(not_fixed[[i]])==0) next
    for(j in 1:length(not_fixed[[i]])){
      #     if(grep("_cov",names(not_fixed[i]))){
      #       k <- c()
      #       l <- 1
      #       while(length(k)<length(not_fixed[[i]])){
      #         k <- c(k,rep(l,l))
      #         l <- l+1
      #       }
      #     }
      if(not_fixed[[i]][j]==1){ 
        if(names(not_fixed[i])=="bpop") parnam <- c(parnam,paste(names(not_fixed[i]),"[",j,"]",sep=""))  
        if(any(names(not_fixed[i])==c("D","SIGMA"))) parnam <- c(parnam,paste(names(not_fixed[i]),"[",j,",",j,"]",sep=""))
        if(length(grep("_cov",names(not_fixed[i])))!=0) parnam <- c(parnam,paste(names(not_fixed[i]),"[",j,"]",sep="")) 
      }
    }
  }
  return(parnam)
}
