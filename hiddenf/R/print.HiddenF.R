print.HiddenF <-
function(x,method="ACMIF",...)
{
if(method %in% c("HiddenF","ACMIF")){
# cat("The ACMIF-Test for Hidden Additivity \n")
hfanova <- anova.HiddenF(x,warncat=FALSE,print=FALSE,return=TRUE)
Fratio <- prettyNum(hfanova$F[4],digits=4)
pvalue <- prettyNum(x$adjpvalue,digits=4)
cat(paste("F=",Fratio," p-value =",pvalue," df=",hfanova$Df[4],",",hfanova$Df[5],"\n",sep=""))
cat(paste("(Bonferroni-adjusted for all",x$cc,"possible configurations) \n"))
}
else if(method %in% c("Tukey","TUKEY","tukey")){
tukey.out <- TukeyPvalue(x)
t.anova <- anova(tukey.out$singledf.out)
pvalue <- prettyNum(tukey.out$pvalue,digits=4)
Fstat <- prettyNum(t.anova$F[3],digits=4)
cat("Tukey's one degree of freedom test for non-additivity:\n")
cat(paste("Test statistic: F=",Fstat," with pvalue p=",pvalue," on df=1,",t.anova$D[4],"\n",sep=""))
#list(pvalue=pvalue,singledf.out=singledf.out)
#list(pvalue=pvalue,Fstat=anova(singledf.out)$F[3],partial.r2=partial.r2)
}
else if(method %in% c("Mandel","MANDEL","mandel")){
m.out <- MandelPvalue(x)
Fstat <- prettyNum(m.out$Fratio,digits=4)
pvalue <- prettyNum(m.out$pvalue,digits=4)
cat("Mandel's rows-linear test for non-additivity \n")
cat(paste("Test statistic F=",Fstat," with pvalue=",pvalue," on DF=",m.out$df[1],",",m.out$df[2],"\n",sep="")) 
}


else if(method %in% c("KKSA","kksa")){
kksa.out <- KKSAPvalue(x)
Fstat <- prettyNum(kksa.out$fmax,digits=4)
pvalue <- prettyNum(kksa.out$pvalue,digits=4)
cat("The Kharrati-Kopaei and Sadoogi-Alvandi test for non-additivity \n")
cat(paste("Test statistic Fmax=",Fstat," with pvalue=",pvalue," on DF=",kksa.out$NumDf,",",kksa.out$DenomDf,"\n",sep="")) 
}
else if(method %in% c("Malik","malik")){
# still working on this one HERE, sept 16, 2015
m.out <- MalikPvalue(x,pnote=FALSE,...)
pvalue <- prettyNum(m.out$pvalue,digits=4)
Tc <- prettyNum(m.out$Tc,digits=4)
cat("Malik's test for non-additivity \n")
cat(paste("Test statistic Tc=",Tc," with pvalue=",pvalue," on ",m.out$N," Monte Carlo simulations \n",sep="")) 
}
else {
stop(paste(method,"is not a supported method \n"))
}
#output <- list(Fratio=Fratio,adjpvalue=pvalue,method="ACMIF")
#class(output) <- "print.HiddenF"
#return(output)

}
