Summaryfactor0<-
function (x, maxsum = 7, pourcent=1) 
{
    k <- length(levels(x))
    Table <- summary(na.omit(x), maxsum = maxsum)
    if (k > maxsum) {
        Table[1:(maxsum - 1)] <- rev(sort(Table[1:(maxsum - 1)]))
    }
    else {
        Table <- rev(sort(Table))
    }
    if(pourcent==0){
return(Table)}
else{
Table<-round(100*Table/sum(Table))
return(Table)
}
}