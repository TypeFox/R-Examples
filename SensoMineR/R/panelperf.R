"panelperf" <- function(donnee,formul,subset=NULL,firstvar,lastvar=ncol(donnee),random=TRUE){

old.contr = options()$contrasts
options(contrasts=c("contr.sum", "contr.sum"))

  for (j in 1 :(firstvar-1))  donnee[,j] <- as.factor(donnee[,j])
  formul = as.formula(formul)
  lab.sauv <- lab <- colnames(donnee)
  for (i in 1:length(lab)) lab[i]=gsub(" ",".",lab[i])
  colnames(donnee) = lab
  equation <- as.character(formul)

  Terms=attr(terms(as.formula(equation)),"term.labels")
  equation = paste("~",Terms[1])
  if (length(Terms) > 1) for (i in 2:length(Terms)) equation <- paste(equation,"+",Terms[i])
  equation <- as.character(as.formula(equation))

  dim.donnee <- ncol(donnee)
  for (i in 1:dim.donnee) {
    if (gsub(" ","",strsplit(equation,split="+",fixed=TRUE)[[2]][1])==lab[i]) col.p <- i
    if (gsub(" ","",strsplit(equation,split="+",fixed=TRUE)[[2]][2])==lab[i]) col.j <- i
  }
  res <- matrix(0,lastvar+1-firstvar,1)
  r2 <- matrix(0,lastvar+1-firstvar,1)
  variab <- perf <- matrix(0,lastvar+1-firstvar,length(strsplit(equation[2],"+",fixed=TRUE)[[1]]))

  for (varendo in firstvar:lastvar) {
      formule <- paste(lab[varendo],"~",equation[2])
    formule <- as.formula(formule)
    aux1 <- aov( formule , data = donnee, subset=subset,na.action =na.exclude)
    aux <- summary(aux1)[[1]]
    perf[varendo-firstvar+1,] <- aux[-nrow(aux),5]
    variab[varendo-firstvar+1,] <- aux[-nrow(aux),2]/sum(aux[,2])
    res[varendo-firstvar+1,] <- sqrt(aux[nrow(aux),3])
    r2[varendo-firstvar+1,] <- summary.lm(aux1)$r.squared
 
    if (random) {
      panelist=colnames(donnee)[col.j]
      for (i in 1:length(Terms)){
        if (any(grep(panelist,Terms[i]))){ 
          if (any(grep(":",Terms[i]))){
            facteur = gsub(":","",Terms[i])
            facteur = gsub(panelist,"",facteur)
            for (k in 1:nrow(aux)) if(gsub(" ","",rownames(aux)[k])==facteur) nrow.facteur = k
            perf[varendo-firstvar+1,nrow.facteur] <- pf(aux[nrow.facteur,3]/aux[i,3],aux[nrow.facteur,1],aux[i,1],lower.tail=FALSE)
          }
        }
      }
    }
  }
aa <- strsplit(as.character(formule),split="~",fixed=TRUE)[[3]]
dimnames(variab) <- dimnames(perf) <- list(lab.sauv[firstvar:lastvar],rownames(aux)[-nrow(aux)])
dimnames(res) <- list(lab.sauv[firstvar:lastvar],"stdev residual")
dimnames(r2) <- list(lab.sauv[firstvar:lastvar],"r2")

panelperf = list() 
panelperf$p.value = perf
panelperf$variability = variab
panelperf$res = res
panelperf$r2 = r2
return(panelperf)
options(contrasts=old.contr)
}
