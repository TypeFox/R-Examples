stMCA <- function(resmca,control) {
  if(attr(resmca,'class')[1] %in% c('MCA','speMCA','csMCA')) temp <- resmca
  if(attr(resmca,'class')[1]=='multiMCA') temp <- resmca$my.mca[[1]]
  if(attr(temp,'class')[1] %in% c('MCA','speMCA')) {
     wt <- temp$call$row.w
     X <- temp$call$X
     covariate <- control
     }
  if(attr(temp,'class')[1] == 'csMCA') {
     wt <- temp$call$row.w[temp$call$subcloud]
     X <- temp$call$X[temp$call$subcloud,]
     covariate <- lapply(control,function(x) x[temp$call$subcloud])
     }
  f <- 'resmca$ind$coord ~ covariate[[1]]'
  if(length(covariate)>1) {
     for(i in 2:length(covariate)) f <- paste(f,paste('covariate[[',i,']]',sep=""),sep="+")
     }
  .env <- environment() ## identify the environment of stMCA
  f <- as.formula(f,env=.env)
  fit <- lm(f,weights=wt)
  res <- residuals(fit)
  z <- data.frame(res,X)
  acp <- PCA(z,quali.sup=(ncol(res)+1):ncol(z),scale.unit=FALSE,graph=FALSE,row.w=wt)
  acp$call$fit <- fit
  acp$call$row.w <- temp$call$row.w
  if(attr(resmca,'class')[1] %in% c('MCA','speMCA','csMCA')) {
    if(attr(resmca,'class')[1] %in% c('speMCA','csMCA')) {
       acp$quali.sup$coord <- acp$quali.sup$coord[-resmca$call$excl,] 
       acp$quali.sup$cos2 <- acp$quali.sup$cos2[-resmca$call$excl,] 
       acp$quali.sup$v.test <- acp$quali.sup$v.test[-resmca$call$excl,]
       acp$call$excl <- resmca$call$excl
       }
    if(attr(resmca,'class')[1] == 'csMCA') acp$call$subcloud <- resmca$call$subcloud
    rownames(acp$quali.sup$coord) <- rownames(resmca$var$coord)
    rownames(acp$quali.sup$cos2) <- rownames(resmca$var$coord)
    rownames(acp$quali.sup$v.test) <- rownames(resmca$var$coord)
    acp$quali.sup$eta2 <- matrix(nrow=ncol(resmca$call$X),ncol=ncol(acp$ind$coord))
    for(i in 1:ncol(X)) acp$quali.sup$eta2[i,] <- apply(acp$ind$coord,2,function(x) summary(lm(x~X[,i],weights=wt))$r.squared)
    dimnames(acp$quali.sup$eta2) <- list(colnames(X),colnames(acp$ind$coord))
    acp$quali.sup$weight <- resmca$var$weight
    acp$var <- acp$quali.sup
    acp$quali.sup <- acp$call$quali.sup <- NULL
    #acp$call$X <- acp$call$X[,-(1:ncp)]
    acp$call$X <- acp$call$X[,-(1:ncol(res))]
    class(acp) <- c('stMCA','list') # new
    acp$call$input.mca <- attr(resmca,'class')[1] # new
    }
  if(attr(resmca,'class')[1]=='multiMCA') {
    class(acp) <- c('stMCA','list') # new
    acp$call$input.mca <- 'multiMCA' # new
    VAR <- list()
    for(i in 1:resmca$call$ngroups) {
      if(attr(resmca$my.mca[[i]],'class')[1] %in% c('MCA','speMCA')) DATA <- resmca$my.mca[[i]]$call$X
      if(attr(resmca$my.mca[[i]],'class')[1] == 'csMCA') DATA <- resmca$my.mca[[i]]$call$X[resmca$my.mca[[i]]$call$subcloud,]
      cond1 <- colSums(apply(dichotom(DATA),2,as.numeric),na.rm=TRUE)>0
      cond2 <- !((1:ncol(dichotom(DATA))) %in% resmca$my.mca[[i]]$call$excl)
      coord <- do.call('rbind',lapply(as.list(colnames(DATA)), function(x) varsup(acp,DATA[,x])$coord))[cond2[cond1],]
      rownames(coord) <- colnames(dichotom(DATA))[cond1 & cond2]
      cos2 <- do.call('rbind',lapply(as.list(colnames(DATA)), function(x) varsup(acp,DATA[,x])$cos2))[cond2[cond1],]
      rownames(cos2) <- rownames(coord)
      vrc <- list()
      for(j in 1:ncol(DATA)) vrc[[colnames(DATA)[j]]] <- varsup(acp,DATA[,j])$var
      #long <- do.call('c',lapply(as.list(colnames(DATA)),function(x) rep(length(DATA[,x]),times=nlevels(DATA[,x]))))[-resmca$my.mca[[i]]$call$excl]
      #v.test <- sqrt(cos2)*sqrt(long-1)
      #v.test <- (((abs(coord)+coord)/coord)-1)*v.test
      v.test <- do.call('rbind',lapply(as.list(colnames(DATA)), function(x) varsup(acp,DATA[,x])$v.test))[cond2[cond1],]
      rownames(v.test) <- rownames(coord)
      VAR[[paste('mca',i,sep='')]] <- list(weight=resmca$my.mca[[i]]$var$weight,coord=round(coord,6),cos2=round(cos2,6),v.test=round(v.test,6),var=vrc)
      }
    acp$VAR <- VAR
    acp$call$quali.sup <- acp$quali.sup <- NULL
    acp$call$X <- acp$call$X[,1:ncol(resmca$ind$coord)]
    acp$call$ngroups <- resmca$call$ngroups
    }
  #acp$eig$mrate <- round(acp$eig[[2]],1)
  #acp$eig$cum.mrate <- cumsum(acp$eig$mrate)    
  RES <- acp
  return(RES)
  }