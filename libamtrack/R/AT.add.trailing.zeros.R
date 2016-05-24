AT.add.trailing.zeros <- function(x, digits = 5){
# Convert 99.1 into "99.10000" etc.
# Library: clan
# Created March 12, 2009, S. Greilich
# corrects for possible round-off error when 
# "digits" < number of significant digits in fractional part of "x"

mx 	<- paste("0.",paste(rep("0",digits),sep="",collapse=""),"1",sep="")
fx	<- substring(as.character(x), grep(".", x) + 2)
if (nchar(fx) >= digits){
	nx	<-	round(x, digits = digits)
	res	<- as.character(nx)
}else{
	nx	<- as.numeric(mx)+x
	res	<- substring(nx, 1, nchar(as.character(nx))-1)	
}
return(res)
}

