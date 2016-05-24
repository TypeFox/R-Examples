decat <- function(donnee,formul,firstvar,lastvar=length(colnames(donnee)),proba = 0.05,graph=TRUE, col.lower = "mistyrose", col.upper = "lightblue", nbrow = NULL, nbcol = NULL, random = TRUE){

    old.contr = options()$contrasts
    options(contrasts=c("contr.sum", "contr.sum"))
    for (j in 1 :(firstvar-1)) donnee[,j] <- as.factor(donnee[,j])
    level.lower = -qnorm(proba/2)
    formul = as.formula(formul)
    lab.sauv <- lab <- colnames(donnee)
    for (i in 1:length(lab)) lab[i]=gsub(" ",".",lab[i])
    colnames(donnee) = lab

    equation <- as.character(formul)

  Terms=attr(terms(as.formula(equation)),"term.labels")
  equation = paste("~",Terms[1])
  if (length(Terms) > 1) for (i in 2:length(Terms)) equation <- paste(equation,"+",Terms[i])
  equation <- as.character(as.formula(equation))

    dim.donnee <- dim(donnee)[2]

    if (length(strsplit(equation,split="+",fixed=TRUE)[[2]]) == 1) random = FALSE # if there is 1 effect, there is not random effect
    for (i in 1:dim.donnee) {
      if (gsub(" ","",strsplit(equation,split="+",fixed=TRUE)[[2]][1])==lab[i]) col.p <- i
      if (random){
        if (gsub(" ","",strsplit(equation,split="+",fixed=TRUE)[[2]][2])==lab[i]) col.j <- i
      }
    }
    nb.modalite <- nlevels(donnee[,col.p])
    don.aux <- cbind.data.frame(donnee,fac=ordered(donnee[,col.p],rev(levels(donnee[,col.p]))))
    dim.don.aux <- dim(don.aux)[2]
    don.aux[,col.p] <- as.factor(don.aux[,dim.don.aux])
    tabF <- matrix(0,lastvar+1-firstvar,2)
    adjmean <- coeff <- tabT <- matrix(0,lastvar+1-firstvar,nb.modalite)
    lab2 <- labels(don.aux)[[2]]
  for (varendo in firstvar:lastvar) {
    formule <- paste(lab[varendo],"~",equation[2])
    formule <- as.formula(formule)
    res <- summary(aov( formule , data = donnee, na.action =na.exclude))[[1]]

    nrow.facteur=nrow(res)
    if (random) {
      panelist=colnames(donnee)[col.j]
      product = colnames(donnee)[col.p]
      for (i in 3:length(Terms)){                ## 1 is product, 2 is panelist
        if ((any(grep(product,Terms[i])))&(any(grep(":",Terms[i])))&(any(grep(panelist,Terms[i])))) nrow.facteur = i
      }
    }

##    tabF[varendo-firstvar+1,1] <- -qnorm(pf(res[1,4],res[1,1],res[dim(res)[1],1],lower.tail=FALSE))
##    tabF[varendo-firstvar+1,2] <-        pf(res[1,4],res[1,1],res[dim(res)[1],1],lower.tail=FALSE)
    tabF[varendo-firstvar+1,1] <- -qnorm(pf(res[1,3]/res[nrow.facteur,3],res[1,1],res[nrow.facteur,1],lower.tail=FALSE))
    tabF[varendo-firstvar+1,2] <-        pf(res[1,3]/res[nrow.facteur,3],res[1,1],res[nrow.facteur,1],lower.tail=FALSE)
    res2 <- summary.lm(aov( formule , data = donnee, na.action =na.exclude))$coef[1:nb.modalite,]
    moy <- res2[1,1]
    res2 <- res2[-1,]
    if (nb.modalite >2){
##      tabT[varendo-firstvar+1,1:(nb.modalite-1)] <-  -qnorm(( pf(res2[,3]^2,1,res[(dim(res)[[1]]),1],lower.tail=FALSE) )/2)*(res2[,1]/abs(res2[,1]))
      tabT[varendo-firstvar+1,1:(nb.modalite-1)] <- -qnorm((pf(res2[,3]^2*(res[nrow(res),3]/res[nrow.facteur,3]) ,1,res[nrow.facteur,1],lower.tail=FALSE) )/2)*sign(res2[,1])
      coeff[varendo-firstvar+1,1:(nb.modalite-1)] <-  res2[,1]
    }
    if (nb.modalite ==2){
##      tabT[varendo-firstvar+1,1:(nb.modalite-1)] <-  -qnorm(( pf(res2[3]^2,1,res[(dim(res)[[1]]),1],lower.tail=FALSE) )/2)*(res2[1]/abs(res2[1]))
      tabT[varendo-firstvar+1,1:(nb.modalite-1)] <- -qnorm((pf(res2[3]^2*(res[nrow(res),3]/res[nrow.facteur,3]) ,1,res[nrow.facteur,1],lower.tail=FALSE) )/2)*sign(res2[1])
      coeff[varendo-firstvar+1,1:(nb.modalite-1)] <-  res2[1]
    }
    res2 <- summary.lm(aov( formule , data = don.aux, na.action =na.exclude))$coef[2,]
##    tabT[varendo-firstvar+1,nb.modalite] <-  -qnorm(( pf(res2[3]^2,1,res[(dim(res)[[1]]),1],lower.tail=FALSE) )/2)*(res2[1]/abs(res2[1]))
      tabT[varendo-firstvar+1,nb.modalite] <- -qnorm((pf(res2[3]^2*(res[nrow(res),3]/res[nrow.facteur,3]) ,1,res[nrow.facteur,1],lower.tail=FALSE) )/2)*sign(res2[1])
    coeff[varendo-firstvar+1,nb.modalite] <-  res2[1]
    adjmean[varendo-firstvar+1,] <- moy+coeff[varendo-firstvar+1,]
  }
  nomdescripteur <- lab.sauv[firstvar:lastvar]
  dimnames(tabF) <- list(nomdescripteur,c("Vtest","P-value"))
  dimnames(adjmean) <-   dimnames(coeff) <- dimnames(tabT) <- list(nomdescripteur,levels(donnee[,col.p]))
  resF <- vector("list",length=1)
  select1 <- (1:nrow(tabF))[tabF[order(tabF[,2]),2]< proba]
  if (length(select1) >0 ){
    resF <- cbind.data.frame(qnorm(tabF[order(tabF[,2]),2],lower.tail=FALSE)[select1],tabF[order(tabF[,2]),2][select1])
    dimnames(resF)[[2]]=c("Vtest","P-value")
    resT <- vector("list",length=nb.modalite)
    for (i in 1:nb.modalite) {
      select <- (1:nrow(tabT))[abs(tabT[rev(order(tabT[,i])),i])>=level.lower]
      resT[[i]] <- cbind.data.frame(coeff[rev(order(tabT[,i])),i][select],adjmean[rev(order(tabT[,i])),i][select],2*(pnorm(-abs(tabT[rev(order(tabT[,i])),i][select]))),tabT[rev(order(tabT[,i])),i][select])
      dimnames(resT[[i]])[[2]]=c("Coeff","Adjust mean","P-value","Vtest")
    }
    names(resT) = c(levels(donnee[,col.p]))
  }
  if (graph){
    par(las=3)
    barplot(tabF[order(tabF[,2]),2],ylim=c(0,1),names.arg=rownames(tabF[order(tabF[,2]),]),ylab="P-value",main="P-value associated with the F-test of the product effet for each descriptor",cex.main=0.9,cex.names=0.8)
    par(las=0)

  }
 
  result = list() 
  result$tabF = tabF
  result$tabT = t(tabT)
  result$coeff = t(coeff)
  result$adjmean = t(adjmean)
  if (length(select1) > 0){
    result$resF = resF
    result$resT = resT
  }
  if (graph) {
    if (is.null(nbrow)) nbrow = nrow(result$adjmean)
    if (is.null(nbcol)) nbcol = ncol(result$adjmean)
    if ((nrow(result$tabT)>2)&(ncol(result$tabT)>1)){
      aux.sort = result$adjmean[rownames(magicsort(result$tabT)),
          colnames(magicsort(result$tabT))]
      aux.col = magicsort(result$tabT)
    }
    else {
      if (nrow(result$tabT)==2){
        aux.sort = result$adjmean[, names(sort(result$tabT[1,]))]
        aux.col = result$tabT[, names(sort(result$tabT[1,]))]
      }
      if (ncol(result$tabT)==1){
        aux.sort = as.matrix(result$adjmean[names(sort(result$tabT[,1])),])
        aux.col = as.matrix(result$tabT[names(sort(result$tabT[,1])),])
      }
    }  
    coltable(aux.sort, aux.col, col.lower = col.lower,
        col.upper = col.upper, level.lower = qnorm(proba/2),
        level.upper = -qnorm(proba/2), nbrow = nbrow, nbcol = nbcol,
        main.title = "Ajusted mean")
  }

  if (length(select1) == 0) print("Warning: No variables are discriminant")
  return(result)
  options(contrasts=old.contr)
}
