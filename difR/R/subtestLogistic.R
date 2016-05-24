subtestLogistic<-function(x,items,groups,alpha=0.05){
if (class(x)!="genLogistic") stop("'x' must be an output of the 'genLogistic' function",call.=FALSE)
if (length(groups)<=1) stop("at least two groups must be specified in 'groups'",call.=FALSE)
if (length(groups)>length(x$focal.names)+1) stop("there are more groups specified than possible",call.=FALSE)
nGroup<-length(x$focal.names)
if (is.character(items)==TRUE){
ITEMS<-NULL
for (i in 1:length(items)) ITEMS[i]<-(1:length(x$genLogistik))[x$names==items[i]]
}
else {
if (is.numeric(items)==TRUE) ITEMS<-items
else stop("'items' must be either numeric or character",call.=FALSE)
}
if (is.character(groups)==TRUE){
GROUPS<-NULL
for (i in 1:length(groups)){
if (sum(groups[i]==x$focal.names)==0) GROUPS[i]<-0
else GROUPS[i]<-(1:nGroup)[groups[i]==x$focal.names]
}
}
else{
if (is.numeric(groups)==TRUE) GROUPS<-groups
else stop("'groups' must be either numeric or character",call.=FALSE)
}
GROUPS<-sort(GROUPS)
nCol<-ifelse(x$type=="udif",2+nGroup,2+2*nGroup)
nCont<-length(groups)-1
nRow<-ifelse(x$type=="both",2*nCont,nCont)
C<-matrix(0,nRow,nCol)
if (min(GROUPS)==0){
if (x$type=="nudif"){
for (i in 1:nCont) C[i,2+nGroup+GROUPS[i+1]]<-1
}
else{
if (x$type=="udif"){
for (i in 1:nCont) C[i,2+GROUPS[i+1]]<-1
}
else{
for (i in 1:nCont) C[i,2+GROUPS[i+1]]<-1
for (i in 1:nCont) C[i+nCont,2+nGroup+GROUPS[i+1]]<-1
}

}
}
else{
if (x$type=="nudif"){
for (i in 1:nCont){
C[i,2+nGroup+GROUPS[1]]<-1
C[i,2+nGroup+GROUPS[1+i]]<-(-1)
}
}
else{
if (x$type=="udif"){
for (i in 1:nCont){
C[i,2+GROUPS[1]]<-1
C[i,2+GROUPS[1+i]]<-(-1)
}
}
else{
for (i in 1:nCont){
C[i,2+GROUPS[1]]<-1
C[i,2+GROUPS[1+i]]<-(-1)
}
for (i in 1:nCont){
C[i+nCont,2+nGroup+GROUPS[1]]<-1
C[i+nCont,2+nGroup+GROUPS[1+i]]<-(-1)
}
}
}
}

stats<-NULL
for(it in 1:length(ITEMS)){
coeff<-as.numeric(x$logitPar[ITEMS[it],])
if (x$type=="udif") coeff<-coeff[1:(2+nGroup)]
Sig<-x$covMat[,,ITEMS[it]]
stats[it]<-t(C %*% coeff) %*% solve(C %*% Sig %*% t(C)) %*% C %*% coeff
}
tab<-cbind(ITEMS,round(stats,3),rep(nRow,length(ITEMS)))
tab<-cbind(tab,round(1-pchisq(tab[,2],tab[,3]),3))
colnames(tab)<-c("Item","Q stat","df","P-value")
RES<-list(stats=tab,contrastMatrix=C,items=items,groups=groups,type=x$type,purification=x$purification,alpha=alpha)
class(RES)<-"subLogistic"
return(RES)}


# METHODS
print.subLogistic<-function (x, ...) 
{
    res <- x
    cat("\n")
    mess1 <- switch(res$type, both = " both types of ", nudif = " nonuniform ", 
        udif = " uniform ")
    cat("Detection of", mess1, "Differential Item Functioning", 
        "\n", "using Generalized logistic regression method,", "\n",sep = "")
    if (res$purification == TRUE) 
        pur <- "with "
    else pur <- "without "
    cat(pur, "item purification and with a subset of ", length(res$groups), " groups ","\n", sep = "")
    cat("of examinees", 
        "\n", "\n", sep = "")
    if (is.character(res$groups) == TRUE | is.factor(res$groups) == 
        TRUE) {
        cat("Groups compared:", "\n")
        nagr <- cbind(res$groups)
        rownames(nagr) <- rep("", nrow(nagr))
        colnames(nagr) <- ""
        print(nagr, quote = FALSE)
        cat("\n")
    }
cat("DIF flagging criterion: Wald test","\n","\n")

    cat("Generalized Logistic regression statistic:", "\n", "\n")

m1<-as.numeric(res$stat[,2:4])
if (length(res$items)==1) m1<-rbind(m1)
else m1<-matrix(m1,length(res$items),length(m1)/length(res$items))

    symb <- symnum(m1[,3], c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
        "**", "*", ".", ""))
    m1 <- noquote(cbind(format(m1, justify = "right"), symb))
if (is.numeric(res$items)==FALSE) rownames(m1)<-res$items
    else {
        rn <- NULL
        for (i in 1:nrow(m1)) rn[i] <- paste("Item", res$items[i], sep = "")
        rownames(m1) <- rn
    }
    colnames(m1) <- c("Wald", "df", "P-value", "")
    print(m1)
    cat("\n")
    cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ", 
        "\n")
    cat("\n", "Detection threshold: ", round(qchisq(1-res$alpha,res$stat[1,3]), 4), " (significance level: ", 
        res$alpha, ")", "\n", "\n", sep = "")
 }