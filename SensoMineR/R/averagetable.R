"averagetable" <- function(donnee,formul,subset=NULL,method="coeff",firstvar,lastvar=ncol(donnee),file=NULL){

if ((method!="coeff")&(method!="mean")) stop(paste("The method",method,"is unknown. Use coeff or mean"))

for (j in 1:(firstvar-1))  donnee[,j] <- as.factor(donnee[,j])

formul = as.formula(formul)
lab.sauv <- lab <- colnames(donnee)
for (i in 1:length(lab)) lab[i] = chartr("(), ?;/:'!$=+\n;{}<>[]-","......................", lab[i])
#for (i in 1:length(lab)) lab[i] = chartr("(), ?;/:'!$§%=+\n;{}<>[]@-",".........................", lab[i])
## for (i in 1:length(lab)) lab[i] = chartr(" '", "..", lab[i])
colnames(donnee) = lab
equation <- as.character(formul)
Terms=attr(terms(as.formula(equation)),"term.labels")
equation = paste("~",Terms[1])
if (length(Terms) > 1) for (i in 2:length(Terms)) equation <- paste(equation,"+",Terms[i])
equation <- as.character(as.formula(equation))

dim.donnee <- ncol(donnee)
for (i in 1:dim.donnee) {
  if (gsub(" ","",strsplit(equation,split="+",fixed=TRUE)[[2]][1])==lab[i]) col.p <- i
  if (gsub(" ","",strsplit(equation,split="*",fixed=TRUE)[[2]][1])==lab[i]) col.p <- i
}
nbprod <- length(levels(donnee[,col.p]))
tab<-matrix(NA,nbprod,lastvar-firstvar+1)
row.names(tab) = levels(donnee[,col.p])

if (method =="mean"){
#  for (j in firstvar:lastvar){
    for (i in 1:nbprod){             
      if (length(subset)==0) tab[i,(firstvar:lastvar)-firstvar+1]<-colMeans(donnee[donnee[,col.p]==levels(donnee[,col.p])[i],firstvar:lastvar],na.rm=TRUE)
      if (length(subset)!=0) tab[i,(firstvar:lastvar)-firstvar+1]<-colMeans(donnee[subset&(donnee[,col.p]==levels(donnee[,col.p])[i]),firstvar:lastvar],na.rm=TRUE)
    }
#  }
}

old.contr = options()$contrasts
if (method =="coeff"){
  options(contrasts=c("contr.sum", "contr.sum"))
  for (varendo in firstvar:lastvar) {
    formule <- as.formula(paste(lab[varendo],"~",equation[2]))
    if (length(subset)!=0) aa=tapply(donnee[subset,varendo],donnee[subset,col.p],mean,na.rm=TRUE)
    if (length(subset)==0) aa=tapply(donnee[,varendo],donnee[,col.p],mean,na.rm=TRUE)
    if (!any(is.na(aa))) {
      aux <- summary.lm(aov( formule, data = donnee, subset=subset,na.action =na.exclude))$coef  
      tab[-nbprod,varendo-firstvar+1] <- aux[2:nbprod,1]
      tab[nbprod,varendo-firstvar+1] <-  - sum(tab[,varendo-firstvar+1],na.rm=TRUE)
      tab[,varendo-firstvar+1] <-  tab[,varendo-firstvar+1]+aux[1,1]
    }
    if (any(is.na(aa))) {
      bb = sum(is.na(aa))
      aux <- summary.lm(aov( formule, data = donnee, subset=subset,na.action =na.exclude))$coef 
      tab[!is.na(aa),varendo-firstvar+1][-nrow(aa[!is.na(aa)])] <- aux[2:(nbprod-bb),1]
      tab[!is.na(aa),varendo-firstvar+1][nrow(aa[!is.na(aa)])] <-  - sum(tab[,varendo-firstvar+1],na.rm=TRUE)
      tab[,varendo-firstvar+1] <-  tab[,varendo-firstvar+1]+aux[1,1]
    }
  }
}

dimnames(tab) = list(levels(donnee[,col.p]),lab.sauv[firstvar:lastvar])
tab=as.data.frame(tab)
if (length(file)!=0) write.csv2(tab,file=file,sep=";",dec=",")
return(tab)
options(contrasts=old.contr)
}
