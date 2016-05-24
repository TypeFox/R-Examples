
summary.bal.stratified.data.frame <- function(object,
                                              ...){

  bal.tab.sum   <- object$bal.test$balance.table.summary
  bal.tab.short <- object$bal.test$balance.table

  
  if (is.null(object$bal.test$p.value)){

    meth <- "diff"
    
    if (dim(object$bal.test$Standardized.difference)[2]==1){
      str.val <-
        t(t(round(object$bal.test$Standardized.difference[2:dim(object$bal.test$Standardized.difference)[1],],3)))
    }else{
      str.val <-
        round(object$bal.test$Standardized.difference[2:dim(object$bal.test$Standardized.difference)[1],],3)
    }   
    bal.tab.detail <- format(data.frame(rbind(round(object$bal.test$Standardized.difference[1,],3),
                                              rep("-----", times=dim(object$bal.test$Standardized.difference)[2]),
                                              str.val,
                                              rep("",times=dim(object$bal.test$Standardized.difference)[2]),
                                              rep("----", times=dim(object$bal.test$Standardized.difference)[2]),
                                              object$bal.test$method),
                                        row.names=c("Before",
                                          "------",
                                          paste("Stratum", seq(1:(dim(object$bal.test$Standardized.difference)[1]-1)), sep=" "),
                                          "",
                                          "---------",
                                          "Scale")))
    cat("\n")    
  }else{

    meth <- "test"
   
    if (dim(object$bal.test$p.value)[2]==1){
      str.val <- t(t(object$bal.test$p.value[2:dim(object$bal.test$p.value)[1],]))
    }else{
      str.val <- object$bal.test$p.value[2:dim(object$bal.test$p.value)[1],]
    }
    
    bal.tab.detail <- format(data.frame(rbind(object$bal.test$p.value[1,],
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
    cat("\n")  
  }

  cov.not <- object$bal.test$covariates.NA
  
  bal.all <- list(bal.sum = bal.tab.sum,
                  bal.sho = bal.tab.short,
                  bal.det = bal.tab.detail,
                  cov.not = object$bal.test$covariates.NA,
                  sig.lev = object$bal.test$alpha,
                  method  = meth)

  class(bal.all) <- "summary.bal.stratified.data.frame"

  bal.all

}


  



