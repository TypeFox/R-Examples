summary.bal.data.frame <- function(object,
                                    ...){

  bal.tab.sum   <- object$bal.test$balance.table.summary
  bal.tab.short <- object$bal.test$balance.table


  if(!is.null(object$match.index)){
    label <- "mat"
    if (is.null(object$bal.test$p.value))
      meth <- "diff"    
    else
      meth <- "pval"
  }else{label <- "str"
        if (is.null(object$bal.test$p.value))
          meth <- "diff"
        else
          meth <- "pval"
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
    if (dim(object$bal.test$Standardized.differences)[2]==1){
      str.val <-
        t(t(round(object$bal.test$Standardized.differences[2:dim(object$bal.test$Standardized.differences)[1],],3)))
    }else{
      str.val <-
        round(object$bal.test$Standardized.differences[2:dim(object$bal.test$Standardized.differences)[1],],3)
    }
  }


  if (label=="str"){
    if (meth=="pval"){
      
      bal.tab.det <- format(data.frame(rbind(object$bal.test$p.value[1,],
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
                                     "Test")))
    }else{
    
      bal.tab.det <- format(data.frame(rbind(round(object$bal.test$Standardized.differences[1,],3),
                                         rep("-----", times=dim(object$bal.test$Standardized.differences)[2]),
                                         str.val,
                                         rep("",times=dim(object$bal.test$Standardized.differences)[2]),
                                         rep("----", times=dim(object$bal.test$Standardized.differences)[2]),
                                         object$bal.test$method),
                                   row.names=c("Before",
                                     "------",
                                     paste("Stratum", seq(1:(dim(object$bal.test$Standardized.differences)[1]-1)), sep=" "),
                                     "",
                                     "---------",
                                     "Scale")))
    }
  }else{ ## label="mat"
    if (meth=="pval"){      
    
      bal.tab.det <- format(data.frame(rbind(object$bal.test$p.value[1,],
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
                                     "Test"))) 
    }else{     
    
      bal.tab.det <- format(data.frame(rbind(round(object$bal.test$Standardized.differences[1,],3),
                                         rep("-----", times=dim(object$bal.test$Standardized.differences)[2]),
                                         str.val,
                                         rep("",times=dim(object$bal.test$Standardized.differences)[2]),
                                         rep("----", times=dim(object$bal.test$Standardized.differences)[2]),
                                         object$bal.test$method),
                                   row.names=c("Before",
                                     "------",
                                     "After",
                                     "",
                                     "-------",
                                     "Scale")))
    }
  } 

  cov.not <- object$bal.test$covariates.NA

  bal.all <- list(bal.sum = bal.tab.sum,
                  bal.sho = bal.tab.short,
                  bal.det = bal.tab.det,
                  cov.not = cov.not,
                  sig.lev = object$bal.test$alpha,
                  method  = meth,
                  label   = label)

  class(bal.all) <- "summary.bal.data.frame"

  bal.all

}


  



