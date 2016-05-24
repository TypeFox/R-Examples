"print.spa" <-
function(x,...){
  if(class(x)!="spa"){
    stop("Error:  x is not of type spa")
  } 
  if(!is.null(cl <- x$call)) {
    names(cl)[2] <- ""
    cat("Call:\n")
    dput(cl)
  }

  tab=x$model$conf
  mea=x$model$measure
  type=x$type
  n=x$mode$dims[1]
  m=x$model$dims[2]

  gcv<-x$control$gcv
  con=switch(gcv,tGCV="transductive",lGCV="labeled",aGCV="approximate transductive",fGCV="full transductive")
  
  cat("\nSequential Predictions Algorithm (SPA) with training=",round(m/n*100,0),"%\n\nType=",type," Parameter=",x$model$parm.est$cvlam," GCV type: ",con)

  if(!is.null(x$model$xstr)){
    cat("\n\nCoefficients:\n")
    print(x$model$xstr$coefs[,1])
  }
  
  if(type=="hard"){
    errm=1-sum(diag(tab))/sum(tab)
    cat("\n\nFit Measures:\n\nFinal Confusion Matrix for Trainnig Data:\n")
    print(tab)
    cat("\nTrain Error:", round(errm, digits=3),"\n\n")
  }
  if(type=="soft"){
    cat("\n\nFit Measures:\n\nRMSE=",round(tab[1],3),"  tGCV=", round(tab[2],3),"  tDF= ",round(tab[3],3),"\n\n")##,"  tSIG=", round(tab[4],3),"\n\n")
  }
  cat("Transductive Parameters:\n\nRegions: ",x$model$dat[1],"  Regularization:",x$model$dat[2],"\n\n")
  if(is.null(x$model$xstr)){
    cat("Unlabeled Data Measures:\n\nTraining:  ",round(mea[1],3),"<=1 (Labeled/Supervised ~1)",
        "\nUnlabeled: ",round(mea[2],3),"<=1 (No Supervised Equivalent)\n\n")
  }
  invisible(x)
}

