refitDIFlasso <-
function(dif.obj){
    
    P <- dif.obj$P
    I <- dif.obj$I
    l <- dif.obj$m
    
  names.y <- dif.obj$names.y
  names.x <- dif.obj$names.x
    
    refit.data <- as.data.frame(dif.obj$refit.matrix)
    names(refit.data) <- c(paste("V",1:(ncol(refit.data)-1),sep=""),"XP1")
    
  
    form1 <- as.formula("~ 0")
    form2 <- as.formula(paste("~ ",paste("V",1:(ncol(refit.data)-1),sep="",collapse="+"),sep=""))
    XP1 <- refit.data$XP1
    
    unres2<-penalized(response=XP1,unpenalized=form1,penalized=form2,lambda1=0,lambda2=0.0001,data=refit.data,model="logistic")
    
    theta <- head(coef(unres2),P)
    beta <- coef(unres2)[(P+1):(P+I-1)]
    beta <- append(beta,0,dif.obj$ref.item-1)
    names(beta)<-names.y
    gamma <- matrix(tail(coef(unres2),length(dif.obj$dif.items)*l),nrow=l)
    dimnames(gamma) <- list(names.x,names.y[dif.obj$dif.items])
    
  returns <- list(theta = theta, beta = beta, gamma = gamma, P = P, I = I, m = l, 
                  ref.item = dif.obj$ref.item, dif.items = dif.obj$dif.items, 
                  names.y = names.y, names.x = names.x)
  
  class(returns) <- "DIFlasso.refit"
  
  return(returns)
}
