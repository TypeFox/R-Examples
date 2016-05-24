summary.CAvariants <-
function(object,printdims=3,digits=3,...) {
cat("\n    SUMMARY",object$catype,  "Correspondence Analysis\n")
cat("\n Names of output objects\n")
print(names(object))
d <- object$r
d <- min(printdims, object$r)
#---------------------------------------------------------------------------
if ((object$catype=="CA")|(object$catype=="NSCA") ){
cat("\n Total inertia ", round(object$inertiasum,digits=digits), "\n\n")
cat("Inertias, percent inertias and cumulative percent inertias of the row and column space\n\n")
print(data.frame(object$inertias))
}

#----------------------------------------------------------------------------------------------
if ((object$catype=="DONSCA")|(object$catype=="DOCA") ){
cat("\n Total inertia ", round(object$inertiasum,digits=digits), "\n\n")
cat("Inertias, percent inertias and cumulative percent inertias of the row space\n\n")
print(data.frame(object$inertias))
cat("Inertias, percent inertias and  cumulative percent inertias of the column space \n\n")
print(data.frame(object$inertias2))
cat("\n Polynomial Components of Inertia \n
** Row Components ** \n")
print(object$comps$compsR)
cat("\n** Column Components ** \n")
print(object$comps$compsC)

}
#-----------------------------------------------------------------------------------------------
if ((object$catype=="SONSCA")|(object$catype=="SOCA") ){
cat("\n Total inertia ", round(object$inertiasum,digits=digits), "\n\n")
cat("Inertias, percent inertias and cumulative percent inertias of the row space\n\n")
print(data.frame(object$inertias))
cat("Inertias, percent inertias and  cumulative percent inertias of the column space \n\n")
print(data.frame(object$inertias2))
cat("\n Polynomial Components of Inertia \n
** Column Components ** \n")
print(object$comps)
}

#############################################################
if ((object$catype=="NSCA")||(object$catype=="DONSCA")||(object$catype=="SONSCA")){
cat("\n    Predictability Index for Variants of Non symmetrical Correspondence Analysis:\n")
cat("\nTau Index predicting from column \n\n")
print(object$tau)
Cstatistic<-(sum(object$DataMatrix)-1)*(nrow(object$DataMatrix)-1)*object$tau
#browser()
pvalueC<-1 - pchisq(Cstatistic, (nrow(object$DataMatrix)-1)*(ncol(object$DataMatrix)-1))
cat("\n C-statistic", Cstatistic, "and p-value", pvalueC, "\n")
}
if ((object$catype=="DOCA")|(object$catype=="DONSCA")){
cat("\n Column standard polynomial coordinates \n")
print(data.frame(object$Cstdcoord[,1:d], row.names=object$collabels), digits=digits)
cat("\n Row standard polynomial coordinates \n")
print(data.frame(object$Rstdcoord[,1:d], row.names=object$rowlabels), digits=digits)
cat("\n Column principal polynomial coordinates \n")
print(data.frame(object$Cprinccoord[,1:d], row.names=object$collabels), digits=digits)
cat("\n Row principal polynomial coordinates \n")
print(data.frame(object$Rprinccoord[,1:d], row.names=object$rowlabels), digits=digits)
}

if ((object$catype=="SOCA")|(object$catype=="SONSCA")){
cat("\n Column standard polynomial coordinates \n")
print(data.frame(object$Cstdcoord[,1:d], row.names=object$collabels), digits=digits)
cat("\n Row standard  coordinates \n")
print(data.frame(object$Rstdcoord[,1:d], row.names=object$rowlabels), digits=digits)
cat("\n Column principal  coordinates \n")
print(data.frame(object$Cprinccoord[,1:d], row.names=object$collabels), digits=digits)
cat("\n Row principal polynomial coordinates \n")
print(data.frame(object$Rprinccoord[,1:d], row.names=object$rowlabels), digits=digits)
}
else{
cat("\n Column standard coordinates \n")
print(data.frame(object$Cstdcoord[,1:d], row.names=object$collabels), digits=digits)
cat("\n Row standard coordinates \n")
print(data.frame(object$Rstdcoord[,1:d], row.names=object$rowlabels), digits=digits)
cat("\n Column principal  coordinates \n")
print(data.frame(object$Cprinccoord[,1:d], row.names=object$collabels), digits=digits)
cat("\n Row principal coordinates \n")
print(data.frame(object$Rprinccoord[,1:d], row.names=object$rowlabels), digits=digits)
}
#cat("\n Inner product of coordinates (first two axes)   \n")
#print(round(object$Trend,digits=digits))
}
