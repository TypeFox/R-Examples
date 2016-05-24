power.f2 <-
function(df1,df2,delta, sig.level=.05){
	n <- df2/(df1+1)+1
	fc <- qf(p=sig.level,df1=df1,df2=(df1+1)*(n-1),lower.tail=FALSE)
	lamda <- (delta^2)*(n*(df1+1))
	v <- (df1+1)*(n-1)
	z1b <- (sqrt(2*(df1+lamda)-((df1+2*lamda)/(df1+lamda)))-
	sqrt((2*v-1)*((df1*fc)/v)))/
	sqrt(((df1*fc)/v)+((df1+2*lamda)/(df1+lamda)))
	pnorm(z1b)
}

