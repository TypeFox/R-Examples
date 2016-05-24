wald.cc=function(b,c,n,alpha){
		tetadach=(b-c)/n
		z=qnorm(1-alpha/2)
		se_w=sqrt(b+c-(b-c)^2/n)/n
		waldiu=tetadach-(z*se_w+1/n)
		waldio=tetadach+(z*se_w+1/n)
		cint=c(waldiu,waldio)
		attr(cint, "conf.level") <- 1-alpha
		rval <- list(conf.int = cint, estimate = tetadach)
		class(rval) <- "htest"
		return(rval)
		}
		