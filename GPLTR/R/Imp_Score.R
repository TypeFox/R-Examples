## Compute the depth of any variable within a tree
.depthvar <- function(treee){
  varrs <- as.character(treee$frame$var)
  nodes <- as.numeric(row.names(treee$frame))
  varrs <- varrs[order(nodes)]
  nodes <- sort(nodes)
  depth = rep(0, max(nodes))
  for(k in nodes[-1]){
    depth[k] <- depth[k/2] +1
  }
  depth <- c(0, depth[depth!=0])
  leafs <- which(varrs == '<leaf>')
  varrs <- varrs[-leafs]
  depth <- depth[-leafs]
  names(depth) <- varrs
  depth2 <- sort(tapply(depth, names(depth), min))
  return(list(depth = depth, depth2 = depth2))  
} 

.chenIMP <- function(treee){
  improvee <- treee$splits[,'improve']
  names(improvee) <- row.names(treee$splits)
  DEPTH1 <- .depthvar(treee)$depth
  DEPTH1 <- DEPTH1[names(improvee)]
  Chen_IMP <- 2^(-DEPTH1)*improvee
#  DEV_IMP <- improvee
  Chen_IMP <- tapply(Chen_IMP, names(Chen_IMP), sum)
  DEV_IMP <- tapply(improvee, names(improvee), sum)
  return(list(Chen_IMP = Chen_IMP, DEV_IMP = DEV_IMP))
}
## Compute the variable importance on a tree by noising up the oob data
.VARIMPPERM <- function(fit,tree, dataoob, eoobt, Y.name, thres){
  Varlist <- sort(unique(as.character(tree$frame$var)))[-1]
  VOOB <- matrix(rep(0, length(Varlist)* length(thres)), ncol = length(thres))
  for(ii in 1: length(Varlist))
  {
    varperm <- dataoob[,Varlist[ii]][sample(1:nrow(dataoob))]
    datpermoob <- dataoob
    datpermoob[,Varlist[ii]] <- varperm
    pred <- predict(fit, newdata = datpermoob, type = 'response')
    pred <- sapply(thres, function(vv) as.numeric(pred >= vv))
    VOOB[ii, ] <- apply(pred, 2, function(wz) mean(dataoob[,Y.name] != wz))
  }
  row.names(VOOB) <- Varlist
  VIMPP <- t(t(VOOB) - eoobt)
  return(VIMPP)
}

## compute the variable importance of the bagging proc using
## random permutation, goodness of split,  depth criterion or depth deviance criterion.
VIMPBAG <- function(BAGGRES, data, Y.name){
  dimtrees <- sapply(BAGGRES$Tree_BAG, function(ttt) sum(ttt$frame$var == '<leaf>') )
  if(length(which(dimtrees == 1))>0){
    BAGGRES$IND_OOB[which(dimtrees == 1)] <- NULL
    BAGGRES$Tree_BAG[which(dimtrees == 1)] <- NULL
    BAGGRES$Glm_BAG[which(dimtrees == 1)] <- NULL
    BAGGRES$OOB_ERRORS_PBP <- BAGGRES$OOB_ERRORS_PBP[, -which(dimtrees == 1)]
  }
  NEWBAG <- length(BAGGRES$IND_OOB)
  DATAOOB <- lapply(BAGGRES$IND_OOB, function(w) return(data[w,])) 
  FITTREEOOBEOOB <- lapply(1:NEWBAG, function(ww) 
  {
    return(list(BAGGRES$Glm_BAG[[ww]], BAGGRES$Tree_BAG[[ww]], DATAOOB[[ww]], BAGGRES$OOB_ERRORS_PBP[, ww]))
  })
  thres <- BAGGRES$CUT
  VIMPPBP <- lapply(FITTREEOOBEOOB, function(uu)
  {
    .VARIMPPERM(uu[[1]], uu[[2]], uu[[3]], uu[[4]], Y.name, thres = thres) 
  })
  
  NAMESPBP <- lapply(VIMPPBP, function(ww) row.names(ww))
  NAMESVAR <- unique(unlist(NAMESPBP))
  #VARR <- unlist(VIMPPBP)
  IMPPER = matrix(rep(0, length(NAMESVAR)* length(thres)), ncol = length(thres))
  ESP2 = matrix(rep(0, length(NAMESVAR)* length(thres)), ncol = length(thres))
  Stderr = matrix(rep(0, length(NAMESVAR)* length(thres)), ncol = length(thres))
  OCCUR = c()
  for(j in 1 : length(NAMESVAR))
  {
    INDJ <- sapply(NAMESPBP, function(uu) is.element(NAMESVAR[j], uu))
    VARRJ <- VIMPPBP[INDJ]
    VARRJJ <- sapply(VARRJ, function(ww) ww[NAMESVAR[j], ])
    IMPPER[j, ] <- apply(VARRJJ, 1, function(wzu) sum(wzu)/NEWBAG)
    ESP2[j, ] <- apply(VARRJJ, 1, function(wzu) sum(wzu^2)/NEWBAG)
    Stderr[j, ] <- sqrt(ESP2[j, ] - (IMPPER[j, ])^2)
    OCCUR[j] <- length(VARRJ)
  }
  rm(ESP2)
  names(OCCUR) <- NAMESVAR
  colnames(IMPPER) <- paste('CUT', 1 : length(thres), sep = '')
  IMPPER <- as.list(data.frame(IMPPER))
  IMPPER <- lapply(IMPPER, function(vv){names(vv) = NAMESVAR ; return(vv)})
  IMPPER <- lapply(IMPPER, function(uu) sort(uu, decreasing = TRUE))
  colnames(Stderr) <- paste('CUT', 1 : length(thres), sep = '')
  Stderr <- as.list(data.frame(Stderr))
  Stderr <- lapply(Stderr, function(vv){names(vv) = NAMESVAR ; return(vv)}) 
  OCCUR <- sort(OCCUR, decreasing = TRUE)

  depthBAG <- lapply(BAGGRES$Tree_BAG, function(treee) .depthvar(treee)$depth2)
  depth <- unlist(depthBAG)
  DEPTH <- sort(tapply(depth, names(depth), mean))
  
  CHEN_BAG <- lapply(BAGGRES$Tree_BAG, function(treez) .chenIMP(treez)$Chen_IMP)
  chen <- unlist(CHEN_BAG)
  #STD_CHEN <- tapply(chen, names(chen), sd)
  CHEN_IMP <- sort(tapply(chen, names(chen), sum)/NEWBAG, decreasing = T)
  CHEN_IMP <- CHEN_IMP/sum(CHEN_IMP)*100
  #STD_CHEN <- STD_CHEN[names(CHEN_IMP)]
  
  DEV_IMP_BAG <- lapply(BAGGRES$Tree_BAG, function(treez) .chenIMP(treez)$DEV_IMP)
  devimp <- unlist(DEV_IMP_BAG)
  #STD_DEV <- tapply(devimp, names(devimp), function(ww) sd(c(rep(0, NEWBAG - length(ww)), ww)))
  DEV_IMP <- sort(tapply(devimp, names(devimp), sum), decreasing = T)
  DEV_IMP <- DEV_IMP/sum(DEV_IMP)*100
  #STD_DEV <- STD_DEV[names(DEV_IMP)]
  return(list(PIS = IMPPER, StdPIS = Stderr, OCCUR = OCCUR, DIS = DEV_IMP, MinDepth = DEPTH, DDIS = CHEN_IMP, dimtrees = dimtrees, EOOB = BAGGRES$EOOB, Bagfinal = NEWBAG))
}

