
# transform base
#   v = value of base 10 to be transformed
#   b = new base
#   l = minimal length of returned array (default is 1)
# return value: array of factors, highest exponent first
baseOf<-function(v,b,l=1) {
	remainder<-v
	i<-l
	ret<-NULL
	while(remainder>0 || i>0) {
		#print(paste("i=",i," remainder=",remainder))
		m<-remainder%%b
		if (is.null(ret)) {
			ret<-m
		}
		else {
			ret<-c(m,ret)
		}
		remainder <- remainder %/% b
		i<-i-1
	}
	return(ret)
}
                                   
