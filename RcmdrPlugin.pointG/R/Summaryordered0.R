Summaryordered0<-
function (x, maxsum = 7, pourcent=1) 
{
    Table <- summary(na.omit(x), maxsum = maxsum)
    if(pourcent==0){
return(Table)}
else{
Table<-round(100*Table/sum(Table))
return(Table)
}
}