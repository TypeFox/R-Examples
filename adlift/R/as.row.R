"as.row" <-
function(x){

y<-NULL

# function y=as.row(x)
# input is row or column vector
# output is row vector
# modified from asrow.m, copyright
# Maarten Jansen 14 november 1997

x<-as.matrix(x)        #MAN:gives x dimensions for ncol,nrow

if (nrow(x)> ncol(x)){
	y<-t(x)
}
else{
	y<-x
}

y

}
