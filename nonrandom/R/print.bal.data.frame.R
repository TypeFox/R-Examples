print.bal.data.frame <- function(x,
                                 ...){

  object <- x
  
  if(!is.null(object$match.index)){
    label <- "mat"

    if (is.null(object$bal.test$p.value)){
      meth <- "diff"
    }else{
      meth <- "pval"
    }
  }else{
    label <- "str"
    if (is.null(object$bal.test$p.value)){
      meth <- "diff"
    }else{
      meth <- "pval"
    }
  }

  if (meth=="pval"){
    if (dim(object$bal.test$p.value)[2]==1){
      str.val <-
        t(t(object$bal.test$p.value[2:dim(object$bal.test$p.value)[1],]))
    }else{
      str.val <-
        object$bal.test$p.value[2:dim(object$bal.test$p.value)[1],]
    }
  }else{
    if (dim(object$bal.test$Stand.diff)[2]==1){
      str.val <-
        t(t(round(object$bal.test$Stand.diff[2:dim(object$bal.test$Stand.diff)[1],],3)))
    }else{
      str.val <-
        round(object$bal.test$Stand.diff[2:dim(object$bal.test$Stand.diff)[1],],3)
    }
  }

  
  cat("\n Summary of balance check: \n\n")
  print(object$bal.test$balance.table.summary)

  if ( length(object$bal.test$covariates.NA)==0 ){
    cat("\n\n Covariates not completely tested: ---\n")
  }else{
    cat("\n\n Covariates not completely tested:\n")
    cat(object$bal.test$covariates.NA, "\n")
  }

  cat("\n\n Detailed balance check (overall): \n\n")
  print(object$bal.test$balance.table)


  if (label=="str"){
    if (meth=="pval"){
      cat(paste("\n\n Detailed balance check (per stratum):\n [p-values from tests (significance level: ",
                object$bal.test$alpha/100, ")]\n\n", sep=""))

      print(format(data.frame(rbind(object$bal.test$p.value[1,],
                                    rep("-----", times=dim(object$bal.test$p.value)[2]),
                                    str.val,
                                    rep("",times=dim(object$bal.test$p.value)[2]),
                                    rep("----", times=dim(object$bal.test$p.value)[2]),
                                    object$bal.test$method),
                              row.names=c("Before",
                                "------",
                                paste("Stratum", seq(1:(dim(object$bal.test$p.value)[1]-1)), sep=" "),
                                "",
                                "---------",
                                "Test"))))
      cat("\n")     
    }else{
      cat(paste("\n\n Detailed balance check (per stratum):\n [standardized differences (cut point: ",
                object$bal.test$alpha, ")]\n\n", sep=""))

      print(format(data.frame(rbind(round(object$bal.test$Stand.diff[1,],3),
                                    rep("-----", times=dim(object$bal.test$Stand.diff)[2]),
                                    str.val,
                                    rep("",times=dim(object$bal.test$Stand.diff)[2]),
                                    rep("----", times=dim(object$bal.test$Stand.diff)[2]),
                                    object$bal.test$method),
                              row.names=c("Before",
                                "------",
                                paste("Stratum", seq(1:(dim(object$bal.test$Stand.diff)[1]-1)), sep=" "),
                                "",
                                "---------",
                                "Scale"))))
      cat("\n")   
    }
  }else{ ## label="mat"
    if (meth=="pval"){      
      cat(paste("\n\n Detailed balance check:\n [p-values from tests (significance level: ",
                object$bal.test$alpha/100, ")]\n\n", sep=""))
      
      print(format(data.frame(rbind(object$bal.test$p.value[1,],
                                    rep("-----", times=dim(object$bal.test$p.value)[2]),
                                    str.val,
                                    rep("",times=dim(object$bal.test$p.value)[2]),
                                    rep("----", times=dim(object$bal.test$p.value)[2]),
                                    object$bal.test$method),
                              row.names=c("Before",
                                "------",
                                "After",
                                "",
                                "-------",
                                "Test"))))
      cat("\n")      
    }else{     
      cat(paste("\n\n Detailed balance check:\n [standardized differences (cut point: ",
                object$bal.test$alpha, ")]\n\n", sep=""))
      
      print(format(data.frame(rbind(round(object$bal.test$Stand.diff[1,],3),
                                    rep("-----", times=dim(object$bal.test$Stand.diff)[2]),
                                    str.val,
                                    rep("",times=dim(object$bal.test$Stand.diff)[2]),
                                    rep("----", times=dim(object$bal.test$Stand.diff)[2]),
                                    object$bal.test$method),
                              row.names=c("Before",
                                "------",
                                "After",
                                "",
                                "-------",
                                "Scale"))))
      cat("\n")    
    }
  } 
}
