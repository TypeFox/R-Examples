wtable <- function(var1,var2=NULL,w=rep.int(1,length(var1)),digits=0,mar=TRUE,na=TRUE) {
  #v1 <- factor(var1)
  if(na==TRUE) {
    ww <- w
    v1 <- factor(var1) #new
    levels(v1) <- c(levels(v1),'NA')
    v1[is.na(v1)] <- 'NA'
    v2 <- var2
    if(!is.null(v2)) {
      v2 <- factor(v2)
      levels(v2) <- c(levels(v2),'NA')
      v2[is.na(v2)] <- 'NA'
      }
    }
  if(na==FALSE & is.null(var2)) {
    ww <- w[!is.na(var1)]
    v1 <- factor(var1[!is.na(var1)]) #new
    }
  if(na==FALSE & !is.null(var2)) {
    ww <- w[!is.na(var1) & !is.na(var2)]
    v1 <- factor(var1[!is.na(var1) & !is.na(var2)]) #new
    v2 <- factor(var2[!is.na(var1) & !is.na(var2)]) #new
    }    
  x <- as.matrix(dichotom(v1,out='numeric'))
  if(is.null(var2)) {
    wtab <- t(x)%*%ww
    wtab <- rbind(wtab,sum(wtab))
    rownames(wtab) <- c(levels(v1),'tot')
    if(mar==FALSE) wtab <- as.matrix(wtab[-length(wtab),])
  } else {
    y <- as.matrix(dichotom(v2,out='numeric'))
    wtab <- t(x)%*%diag(ww)%*%y
    wtab <- cbind(wtab,rowSums(wtab))
    wtab <- rbind(wtab,colSums(wtab))
    dimnames(wtab) <- list(c(levels(v1),'tot'),c(levels(v2),'tot'))
    if(mar==FALSE) wtab <- wtab[-nrow(wtab),-ncol(wtab)]
    }
  wtab <- round(wtab,digits)
  return(wtab)
  }
