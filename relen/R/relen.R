relen <-
function (x) {z<-(-1/log(length(levels(x))))*(sum(prop.table(table(x))*log(prop.table(table(x))))); if(length(levels(x))==1) {z<- 0}; return(z) }
