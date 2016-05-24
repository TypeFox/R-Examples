"as.column" <-function(x){

y<-NULL

# function: y<-as.column(x)
# input is row or column vector
# output is column vector
# modified from ascolumn.m, copyright
# Maarten Jansen 7 november 1997

x<-as.matrix(x)     #MAN:gives x dimensions for ncol,nrow

if (nrow(x)<ncol(x)){
	y<-t(x) 
}
else {
	y<-x
}

y
}
