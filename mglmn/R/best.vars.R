best.vars <-
function(x){
temp<-ifelse(x[[1]][1,-1:-5]==1,colnames(x[[1]][,-1:-5]),NA)
temp2<-NULL
j<-1
for (i in 1:length(temp)){
if(is.na(temp[i])==F)
{temp2[j]<-temp[i]
j <- j+1}
}
temp2
}
