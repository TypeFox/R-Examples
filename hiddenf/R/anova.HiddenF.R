anova.HiddenF <-
function(object,warncat=TRUE,method="HiddenF",return=FALSE,print=TRUE,stars=FALSE,...)
{
tall <- object$tall
y <- tall$y
row <- tall$row
col <- tall$col
a <- nlevels(row)
b <- nlevels(col)
if(method %in% c("HiddenF","ACMIF")){
group <- as.factor(object$config.vector)
lm.out <- lm(y~group*col+row/group)
aov.obj <- anova(lm.out)
aov.obj <- rbind(aov.obj,c(sum(aov.obj$Df),sum(aov.obj[,"Sum Sq"]),NA,NA,NA))
row.names(aov.obj)[dim(aov.obj)[1]] <- "C.Total"
# cat("The ACMIF test for the hidden additivity form of interaction\n")
#if(print){stats::print.anova(aov.obj,signif.stars=stars,...)}
if(print){
cat("The ACMIF test for the hidden additivity form of interaction\n")
anovaprint(aov.obj,signif.stars=stars,...)}
#if(print){print(aov.obj)}
if(warncat){cat("(Pvalues in ANOVA table are NOT corrected for multiplicity.) \n")}
#if(return){return(anova(lm.out))}
if(return){return(aov.obj)}
}
else if(method %in% c("TUKEY","Tukey","tukey")){
non_additivity <- fitted(lm(y ~ row + col))^2
#return(anova(lm(y~row+col+psq)))
lm.out <- lm(y~row+col+non_additivity)
aov.obj <- anova(lm.out)
aov.obj <- rbind(aov.obj,c(sum(aov.obj$Df),sum(aov.obj[,"Sum Sq"]),NA,NA,NA))
row.names(aov.obj)[dim(aov.obj)[1]] <- "C.Total"
cat("Tukey's one degree-of-freedom test for multiplicative interaction \n") 
# print(anova(lm(y~row+col+non_additivity)))
#if(print){stats::print.anova(aov.obj,signif.stars=stars,...)}
if(print){anovaprint(aov.obj,signif.stars=stars,...)}
#if(print){print(aov.obj)}
#if(return){return(anova(lm(y~row+col+non_additivity)))}
if(return){return(aov.obj)}
}
else if(method %in% c("Mandel","MANDEL","mandel")){
ymtx <- matrix(y,nrow=a,ncol=b,byrow=T)
coldevs <- apply(ymtx,2,mean)-mean(ymtx)
rowdevs <- apply(ymtx,1,mean)-mean(ymtx)
SSRow <- b*sum(rowdevs^2)
SSCol <- a*sum(coldevs^2)
SSTot <- (a*b-1)*var(y)
slopes <- ymtx %*% coldevs/sum(coldevs^2)
SSMandel <- sum((slopes-1)^2) * sum(coldevs^2)
SSE <- SSTot-SSMandel-SSRow-SSCol
dfE <- (a-1)*(b-2)
MSE <- SSE/dfE
SumSq <- format(c(SSRow,SSCol,SSMandel,SSE,SSTot),digits=4)
MeanSq <- rep(NA,5)
MeanSq[1:4] <- format(c(SSRow,SSCol,SSMandel,SSE)/c(a-1,b-1,a-1,(a-1)*(b-2)),digits=4)
Frow <- (SSRow/(a-1))/MSE
Fcol <- (SSCol/(a-1))/MSE
Fratio <- (SSMandel/(a-1))/MSE
FValue <- rep(NA,5)
FValue[1:3] <- format(c(Frow,Fcol,Fratio),digits=4)
pvalue <- 1-pf(Fratio,(a-1),(a-1)*(b-2))
pvalues <- rep(NA,5)
pvalues[1:3] <- format(c(1-pf(Frow,a-1,dfE),1-pf(Fcol,a-1,dfE),pvalue),digits=4)
#pvalues <- format(c(1-pf(Frow,a-1,dfE),1-pf(Fcol,a-1,dfE),pvalue,NA),digits=4)
Df <- c(a-1,b-1,a-1,(a-1)*(b-2),a*b-1)
table1 <- cbind(Df,"Sum Sq"=SumSq,"Mean Sq"=MeanSq,"F Value"=FValue,"Pr(>F)"=pvalues)
table.df <- as.data.frame(table1,row.names=c("Rows","Columns","Slopes","Residual","C.Total"))
#table.df <- as.table(table1,row.names=c("Rows","Columns","Slopes","Residual"))
if(print){print(table.df,na.print="")}
if(return){return(return(table.df))}
}
else if(method %in% c("KKSA","kksa")){
KKSA.pvalue <- KKSAPvalue(object)
df1 <- KKSA.pvalue$NumDf
df2 <- KKSA.pvalue$DenomDf

l1 <- (1:a)[KKSA.pvalue$grp.vector==1]
l2 <- (1:a)[KKSA.pvalue$grp.vector==2]

index1 <- row %in% l1
index2 <- row %in% l2

y1 <- y[index1]
rows <- row[index1]
cols <- col[index1]
#anova1 <- anova(lm(y[index1] ~ row[index1]+col[index1]))
#anova2 <- anova(lm(y[index2] ~ row[index2]+col[index2]))
aov.obj1 <- anova(lm(y1 ~ rows+cols))
aov.obj1 <- rbind(aov.obj1,c(sum(aov.obj1$Df),sum(aov.obj1[,"Sum Sq"]),NA,NA,NA))
row.names(aov.obj1)[dim(aov.obj1)[1]] <- "C.Total"
#anova1 <- anova(lm(y1 ~ rows+cols))

y2 <- y[index2]
rows <- row[index2]
cols <- col[index2]
#anova2 <- anova(lm(y2 ~ rows + cols))
aov.obj2 <- anova(lm(y2 ~ rows+cols))
aov.obj2 <- rbind(aov.obj2,c(sum(aov.obj2$Df),sum(aov.obj2[,"Sum Sq"]),NA,NA,NA))
row.names(aov.obj2)[dim(aov.obj2)[1]] <- "C.Total"

#fstat <- anova1$M[3]/anova2$M[3]
fstat <- aov.obj1$M[3]/aov.obj2$M[3]
fmax <- max(fstat,1/fstat)

pvalue1 <- 1-pf(fmax,df1,df2)+pf(1/fmax,df1,df2)
cc <- 2^(b-1)-1-b
pvalue <- cc*pvalue1
pvalue <- min(1,pvalue)
cat("Kharrati-Kopaei & Sadooghi-Alvandi's test for interaction \n")
cat("Rows in the first group are: ")
cat(l1)
cat("\n")
cat("Rows in the second group are: ")
cat(l2)
cat("\n")
#print(anova1)
#if(print){stats::print.anova(aov.obj1,signif.stars=stars,...)}
if(print){anovaprint(aov.obj1,signif.stars=stars,...)}
#if(print){print(aov.obj1)}
cat("-----------------------------------------------------\n")
#print(anova2)
#if(print){stats:::print.anova(aov.obj2,signif.stars=stars,...)}
if(print){anovaprint(aov.obj2,signif.stars=stars,...)}
cat("-----------------------------------------------------\n")
# KKSA.pvalue <- round(KKSAPvalue(ymtx.out)$pvalue,digits=4)
cat(paste("The multiplicity-adjusted pvalue from KKSA's ratio of \n"))
cat(paste("Residual Mean Squares is ",round(KKSA.pvalue$pvalue,digits=4),"\n"))
# if(warncat){cat("(Pvalues in ANOVA table are NOT corrected for multiplicity.) \n")}
if(return){return(list(anova1=aov.obj1,anova2=aov.obj2))}
}
else if(method %in% c("Malik","MALIK","malik")){
cat(paste("The anova() function is not supported for method=",method,"\n"))
}
else {
# cat(paste("Error: ",method,"is not a supported method \n"))
stop(paste(method,"is not a supported method \n"))
}
}
