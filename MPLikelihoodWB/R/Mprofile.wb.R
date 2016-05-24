Mprofile.wb <-
function(formula,censor,data,initial=1){
              Y <- model.frame(formula,data=data)[,1]  
              X <- model.matrix(formula,data=data)
          delta <- data[[censor]]            
     design.len <- ncol(X)
                  initial.val <- NULL
          if(length(initial)==1){
         initial.val <- rep(initial,design.len+1)
             }else{
         initial.val <- initial
            }
parm.beta <- NULL
for(i in 1:design.len){
opt.bi <- optim(par=initial.val, mplik.wb.bi, Y=Y,X=X,delta=delta,whc=i,method="BFGS")
parm.beta[i] <- opt.bi$par[i+1]
}
names(parm.beta) <- colnames(X)
opt.sigma <- optim(par=initial.val,mplik.wb.s,Y=Y,X=X,delta=delta,method="BFGS")
scaleP <- exp(opt.sigma$par[1])
names(scaleP) <- "Scale"
results <- list(formula,parm.beta,scaleP)
names(results) <- c("Formula","Coefficients","Scale")
return(results)
}
