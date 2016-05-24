ict <-
function(crosstab){
crosstab=as.matrix(crosstab)
n1=1:nrow(crosstab)
n2=1:ncol(crosstab)
d1=c()
d2=c()
for(i in n1)
for(j in n2){
count=crosstab[i,j]
if(is.na(count)) warning(paste("NA in row",i,", column",j,"exist!"),immediate. = TRUE) 
else if(count>0) for(k in 1:count) {
d1=c(d1,n1[i])
d2=c(d2,n2[j])
} # end if
} # end for j
d1=factor(d1,labels=rownames(crosstab))
d2=factor(d2,labels=colnames(crosstab))
list(d1=d1,d2=d2)
}
