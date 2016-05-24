dimdesc.MCA <- function(resmca,ncp=3,proba=0.05) {
   res <- list()
   X <- resmca$call$X
   classe <- class(resmca)[1] # new
   if(classe=='stMCA') classe=resmca$call$input.mca # new
   for(i in 1:ncol(X)) levels(X[,i]) <- paste(names(X)[i],levels(X[,i]),sep='.')
   for(i in 1:ncp) {
      if(classe %in% c('MCA','speMCA')) temp <- condes(data.frame(resmca$ind$coord[,i],X),1,weights=resmca$call$row.w,proba=proba) # new
      if(classe == 'csMCA') temp <- condes(data.frame(resmca$ind$coord[,i],X[resmca$call$subcloud,]),1,weights=resmca$call$row.w[resmca$call$subcloud],proba=proba) # new
      res[[i]] <- temp
      }
   names(res) <- paste('axe',1:ncp,sep='.')
   return(res)
   }
