test.Null <-
function( trtsel, alpha){
 
 # 6/6/12 changed from test based on the bootstrap to the likelihood ratio test. 
 # below is the old code

 #n.boot <- ncol(boot.data)
 #potential.pvals <- (1:n.boot)/n.boot

 #a3.b <- boot.data[4,]


 #if(!bounded){ 
   
   # first check to see if the max and min cover 0,
   # if this is true the pvalue is less than our lowest possible detectable value

 #  if(!cover(min(a3.b), max(a3.b), 0) ){ 
     
 #    p.val<- paste("<", 1/n.boot)
 #    reject <- TRUE

 #  }else{

 #    reject.all <- unname( mapply( cover, 
 #                                  quantile(a3.b, potential.pvals/2, , study.design = 1, na.rm = TRUE),
 #                                  quantile(a3.b, 1 - potential.pvals/2, study.design = 1, na.rm = TRUE), 
 #                                  rep(0, n.boot))  )
 #    tmp <- which(reject.all==FALSE)[1] 
 #    p.val <- potential.pvals[ifelse(tmp==1, 1, tmp - 1)]
 #    reject <- p.val <= alpha

 #   }

 #  }

  if(is.null(trtsel$model.fit$coefficients)) p.val <- NA
  else p.val <- trtsel$model.fit$coefficients[4,4]
  reject <- p.val <= alpha
  z.value <- trtsel$model.fit$coefficients[4,3] 

  list( reject = reject, p.value = p.val, z.statistic = z.value, alpha = alpha)

}
