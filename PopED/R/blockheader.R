#' Header function for optimization routines
#' 
#' Create some output to the screen and a text file that summarizes the problem you are tying to solve.
#' 
#' @inheritParams RS_opt
#' @inheritParams evaluate.fim
#' @inheritParams blockexp
#' @inheritParams Doptim
#' @inheritParams create.poped.database
#' @inheritParams RS_opt_gen
#' @param name The name used for the output file. Combined with \code{name_header} and \code{iter}. 
#' If \code{""} then output is to the screen.
#' @param iter The last number in the name printed to the output file, combined with \code{name}.
#' @param name_header The initial portion of the file name.
#' @param file_path The path to where the file should be created.
#' @param header_flag Should the header text be printed out?
#' @param ... Additional arguments passed to further functions.
#' 
#' @family Helper
#' @return fn A file handle (or \code{''} if \code{name=''})
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_blockheader.R
#' @keywords internal
#' @export
## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

blockheader <- function(poped.db,name="Default",iter=NULL,
                          e_flag=!(poped.db$settings$d_switch),opt_xt=poped.db$settings$optsw[2],
                          opt_a=poped.db$settings$optsw[4],opt_x=poped.db$settings$optsw[3],
                          opt_samps=poped.db$settings$optsw[1],opt_inds=poped.db$settings$optsw[5],
                          fmf=0,dmf=0,bpop=NULL,d=NULL,docc=NULL,sigma=NULL,
                          name_header=poped.db$settings$strOutputFileName,
                          file_path=poped.db$settings$strOutputFilePath,
                          out_file=NULL,compute_inv=TRUE,
                          trflag=TRUE,
                          header_flag=TRUE,
                          ...)
{
  # BLOCKHEADER_2
  #   filename to write to is 
  #   poped.db$settings$strOutputFilePath,poped.db$settings$strOutputFileName,NAME,iter,poped.db$settings$strOutputFileExtension
  
  #   if((bDiscreteOpt)){
  #     tmpfile=sprintf('%s_Discrete_%g%s',poped.db$settings$strOutputFileName,iter,poped.db$settings$strOutputFileExtension)
  #   } else {
  #     tmpfile=sprintf('%s_RS_SG_%g%s',poped.db$settings$strOutputFileName,iter,poped.db$settings$strOutputFileExtension)
  #   }
  
  #tmpfile=sprintf('%s_%s_%g%s',poped.db$settings$strOutputFileName,name,iter,poped.db$settings$strOutputFileExtension)

  if(!trflag) return('')
  
  if(!is.null(out_file)){
    fn <- out_file
    if(!any(class(fn)=="file") && (fn!='')){
      fn=file(fn,'w')
      if(fn==-1){
        stop(sprintf('output file could not be opened'))
      }
    }
  } else if(name!=""){
    tmpfile <- name_header
    if(name!="Default") tmpfile=paste(tmpfile,"_",name,sep="")
    if(!is.null(iter))  tmpfile=paste(tmpfile,"_",iter,sep="")
    tmpfile=paste(tmpfile,".txt",sep="")
    #tmpfile=sprintf('%s_%s.txt',name_header,name)
    #if(!is.null(iter)) tmpfile=sprintf('%s_%s_%g.txt',name_header,name,iter)
    tmpfile = fullfile(poped.db$settings$strOutputFilePath,tmpfile)
    fn=file(tmpfile,'w')
    if((fn==-1)){
      stop(sprintf('output file could not be opened'))
    }
  } else {
    fn <- ''
    #     filename=readline("File to open for output: ")
    #     fn = file(filename, 'w')
    #     if((fn == -1)){
    #       stop(sprintf('output file could not be opened'))
    #     }
  }
  
  if(!header_flag) return(fn)
  
  
  #tic()
  tic(name=".poped_total_time")
  
  # -------------- LOG FILE: initial status
  if(name=="RS"){
    alg_name <- "Adaptive Random Search"
    if(fn!="") fprintf(fn,'PopED Optimization Results for the %s Algorithm \n\n',alg_name)
  }  else {
    if(fn!="") fprintf(fn,'PopED Results \n\n')
  }
  if(fn!="") fprintf(fn,'        ')
  if(fn!="") fprintf(fn,datestr_poped(poped.db$settings$Engine$Type))
  if(fn!="") fprintf(fn,'\n\n')
  
  if(fn!="" || trflag>1) blockexp(fn,poped.db,
                                   e_flag=e_flag,opt_xt=opt_xt,
                                   opt_a=opt_a,opt_x=opt_x,
                                   opt_samps=opt_samps,opt_inds=opt_inds)
  
  if(dmf!=0 || fmf != 0){ 
    fprintf(fn,paste0("===============================================================================\n",
                      "Initial design evaluation\n"))
    if(fn!="") fprintf(paste0("===============================================================================\n",
                          "Initial design evaluation\n"))
  }
  if(dmf!=0) fprintf(fn,'\nInitial OFV = %g\n',dmf)
  if(dmf!=0 && fn!="") fprintf('\nInitial OFV = %g\n',dmf)
  
  if(dmf!=0 && (fn!="" || trflag>1)){
    output <- get_unfixed_params(poped.db)
    npar <- length(output$all)
    
    fprintf(fn,'\nEfficiency criterion [usually defined as OFV^(1/npar)]  = %g\n',
            ofv_criterion(dmf,npar,poped.db))
    if(fn!=""){
      fprintf('\nEfficiency criterion [usually defined as OFV^(1/npar)]  = %g\n',
              ofv_criterion(dmf,npar,poped.db))
    }
  }
  
  if(is.matrix(fmf) && compute_inv){
    param_vars=diag_matlab(inv(fmf))
    returnArgs <-  get_cv(param_vars,bpop,d,docc,sigma,poped.db) 
    params <- returnArgs[[1]]
    param_cvs <- returnArgs[[2]]
    
    
      
    #fprintf(fn,'\nEfficiency criterion [usually defined as OFV^(1/npar)]  = %g\n',dmf^(1/length(params)))
    #fprintf(fn,'\nEfficiency criterion [usually defined as OFV^(1/npar)]  = %g\n',
    #        ofv_criterion(dmf,length(params),poped.db))
    
    parnam <- get_parnam(poped.db)
    fprintf(fn,'\nInitial design expected parameter \nrelative standard error (%sRSE)\n','%')
    if(fn!="") fprintf('\nInitial design expected parameter \nrelative standard error (%sRSE)\n','%')
    df <- data.frame("Parameter"=parnam,"Values"=params,"RSE_0"=t(param_cvs*100))
    print(df,digits=3, print.gap=3,row.names=F)
    if(fn!="") capture.output(print(df,digits=3, print.gap=3,row.names=F),file=fn)
    fprintf('\n')
    if(fn!="") fprintf(fn,'\n')
    
  }
  
  if(fn!="" || trflag>1) blockopt(fn,poped.db,opt_method=name)
  if(fn!="" || trflag>1) blockother(fn,poped.db,d_switch=!e_flag)
  
  if(fn!="" || trflag) blockoptwrt(fn,poped.db$settings$optsw, opt_xt=opt_xt,
                                     opt_a=opt_a,opt_x=opt_x,
                                     opt_samps=opt_samps,opt_inds=opt_inds)
  
  #fprintf('\n')
  #if(fn!="") fprintf(fn,'\n')
  
  return( fn) 
}
