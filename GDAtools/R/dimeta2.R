dimeta2 <- function(resmca,l,n=names(l),dim=1:resmca$call$ncp) {
  eta2 <- matrix(nrow=length(l),ncol=length(dim))
  sub <- rep(TRUE,times=nrow(resmca$ind$coord))
  if(class(resmca)[1]=='csMCA') sub <- resmca$call$subcloud
  if(class(resmca)[1]=='stMCA') if(resmca$call$input.mca=='csMCA') sub <- resmca$call$subcloud #new
  ww <- resmca$call$row.w[sub] #new
  if(class(resmca)[1]=='stMCA') ww <- resmca$call$fit$weights #new
  for(i in 1:length(dim)) {
    for(j in 1:length(l)) eta2[j,dim[i]] <- summary(lm(resmca$ind$coord[,dim[i]]~l[[j]][sub],weights=ww))$r.squared
   }
  rownames(eta2) <- n
  colnames(eta2) <- colnames(resmca$var$contrib)[dim]
  return(round(100*eta2,1))
  }