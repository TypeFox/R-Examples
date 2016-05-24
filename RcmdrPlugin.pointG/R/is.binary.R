is.binary<-
function(x){
if(length(table(x))==2){binary<-TRUE}
else {binary<-FALSE}
binary
}
