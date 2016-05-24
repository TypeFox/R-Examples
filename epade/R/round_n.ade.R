round_n.ade <-
function(x, digits=0){
x<-unlist(x)
if(!is.vector(x) & !is.matrix(x)) stop('x must be vector or matrix')

#########################################################
# split in single numbers
if(length(x)>1){

if(is.vector(x)){
y<-x
for(i in 1:length(x)){
y[i]<- round_n.ade(x[i], digits)
}
return(y)
}

if(is.matrix(x)){
M<-x
for(i in 1:dim(x)[1]){
for(j in 1:dim(x)[2]){
M[i, j]<- round_n.ade(x[i, j], digits)
}}
return(M)
}

}
#########################################################


if(length(x)==1){
if(is.na(x))  return(x)
if(x==0)      return(as.character(x))
y<-base::format(base::round(x, digits = digits), nsmall=digits, scientific = FALSE)
return(y)
}

}
