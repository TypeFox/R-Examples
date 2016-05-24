utils::globalVariables(c("dmp.txt"))

#' @export
phxnlme <-
  function(inst.path=NULL, path, model.file, cols.file, data, method=5, iterlimit=200){
    
    setwd(path)
    if(is.null(inst.path)){
      inst.path = "C:/Program Files (x86)/Pharsight/Phoenix"
    }
    file.copy(paste(inst.path,"/application/Examples/NLME Command Line/Model 1/RunNLME.bat",sep=""), path)
    batdir2 = paste("\"",path,"/RunNLME.bat\"",sep="")
    cmd.line2 = paste(batdir2,method,iterlimit,model.file,cols.file,data)
    system(cmd.line2,invisible=FALSE,wait=TRUE)
    if(file.exists("out.txt")) file.rename("out.txt","out0001.txt")
    
    modl = list(model.file=model.file, 
                cols.file=cols.file ,
                data=data ,
                method=method,
                iterlimit=iterlimit,bat.file="RunNLME.bat") 
    
    write.csv(modl,file="model.spec.csv")
    
    source("dmp.txt")
      
    ## Summary of results
    theta = data.frame(dmp.txt$coefficients$fixed)
    names(theta) = "Estimates"    
      
    omega = data.frame(diag(dmp.txt$omega))
    names(omega) = "Estimates" 
      
    summary.results = rbind(theta,omega)
    write.csv(summary.results,"summary.results.csv")
      
    ### Create results folder and copy files in
    if(file.exists("Results")){
      res = paste(path,"summary.results.csv",sep="\\")
      file.copy(res,paste(path,"Results",sep="\\"),overwrite=TRUE)
    }else{
      dir.create("Results")
      res = paste(path,"summary.results.csv",sep="\\")
      file.copy(res,paste(path,"Results",sep="\\"))
    }    
  }
