multiMCA <- function(l_mca,ncp=5,compute.rv=FALSE) {
  co <- data.frame(lapply(l_mca,function(x) x$ind$coord/x$eig[[1]][1]))
  if(attr(l_mca[[1]],'class')[1] %in% c('MCA','speMCA')) wt <- l_mca[[1]]$call$row.w
  if(attr(l_mca[[1]],'class')[1] == 'csMCA') wt <- l_mca[[1]]$call$row.w[l_mca[[1]]$call$subcloud]
  afm <- PCA(co,scale.unit=FALSE,row.w=wt,ncp=ncp,graph=FALSE)
  afm$call$row.w <- wt
  attr(afm,'class') <- c('multiMCA','list')
  ngroups <- length(l_mca)
  VAR <- list()
  for(i in 1:ngroups) {
    if(attr(l_mca[[i]],'class')[1] %in% c('MCA','speMCA')) DATA <- l_mca[[i]]$call$X
    if(attr(l_mca[[i]],'class')[1] == 'csMCA') DATA <- l_mca[[i]]$call$X[l_mca[[i]]$call$subcloud,]
    cond1 <- colSums(apply(dichotom(DATA),2,as.numeric),na.rm=TRUE)>0
    cond2 <- !((1:ncol(dichotom(DATA))) %in% l_mca[[i]]$call$excl)
    coord <- do.call('rbind',lapply(as.list(colnames(DATA)), function(x) varsup(afm,DATA[,x])$coord))[cond2[cond1],]
    rownames(coord) <- colnames(dichotom(DATA))[cond1 & cond2]
    cos2 <- do.call('rbind',lapply(as.list(colnames(DATA)), function(x) varsup(afm,DATA[,x])$cos2))[cond2[cond1],]
    rownames(cos2) <- rownames(coord)
    vrc <- list()
    for(j in 1:ncol(DATA)) vrc[[colnames(DATA)[j]]] <- varsup(afm,DATA[,j])$var
    long <- do.call('c',lapply(as.list(colnames(DATA)),function(x) rep(length(DATA[,x]),times=nlevels(DATA[,x]))))[-l_mca[[i]]$call$excl]
    v.test <- sqrt(cos2)*sqrt(long-1)
    v.test <- (((abs(coord)+coord)/coord)-1)*v.test
    rownames(v.test) <- rownames(coord)
    VAR[[paste('mca',i,sep='')]] <- list(weight=l_mca[[i]]$var$weight,coord=round(coord,6),cos2=round(cos2,6),v.test=round(v.test,6),var=vrc)
    }
  afm$VAR <- VAR
  agg <- factor()
  for(i in 1:ngroups) agg <- c(agg,rep(i,times=l_mca[[i]]$call$ncp))
  contrib <- do.call('rbind',by(afm$var$contrib,agg,colSums))
  rownames(contrib) <- paste('mca',1:length(l_mca),sep='')
  correl <- do.call('rbind',(lapply(l_mca,function(x) diag(cor(x$ind$coord[,1:ncp],afm$ind$coord[,1:ncp])))))
  rownames(correl) <- paste('mca',1:length(l_mca),sep='')
  colnames(correl) <- paste('Dim',1:ncp,sep='.')
  afm$group <- list(contrib=round(contrib,2),correl=round(correl,3))
  afm$my.mca <- l_mca
  afm$call$ngroups <- ngroups
  afm$eig <- as.list(afm$eig)
  l <- lapply(l_mca,function(x) x$ind$coord)
  l[[length(l)+1]] <- afm$ind$coord
  if(compute.rv==TRUE) { #new
      rv <- matrix(0,nrow=length(l),ncol=length(l))
      for(i in 2:length(l)) {
	for(j in 1:(i-1)) rv[i,j] <- coeffRV(l[[i]],l[[j]])$rv
	}
      rv <- rv+t(rv)
      diag(rv) <- 1
      rownames(rv) <- c(paste('mca',1:length(l_mca),sep=''),'mfa')
      colnames(rv) <- rownames(rv)
      afm$RV <- rv
      } # new
  return(afm)
  }
