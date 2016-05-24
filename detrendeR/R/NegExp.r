NegExp = function(Y,  pos.slope = FALSE){

	
y<-Y[!is.na(Y)] 			#
replace(y, y==0, 0.001)->y		#

	  fit.linear.model = FALSE

	  x = 1:length(y)
 	  a = mean(y[1:floor(length(y) * 0.05)])
        b = -0.01
        k = mean(y[floor(length(y) * 0.95):length(y)])

        try(nls(y ~ Const + A * exp(B * x), start= list(A=a,B=b,Const=k ), trace=F), silent = T)->fits
 
       if (class(fits) == "try-error") {fit.linear.model = TRUE} else fits = predict(fits)
	
       if (fits[1] < fits[length(fits)]) fit.linear.model = TRUE
       if (fits[length(fits)] < 0)  fit.linear.model = TRUE

	 if( fit.linear.model )
       {
       linear.model = lm (y ~ x)
       fits = predict(linear.model)
       if (coef(linear.model)[2] > 0 & !pos.slope) fits = rep(mean(y), length(x))
       }
Y[!is.na(Y)]<-fits #
       #return(fits)
return(Y) #
}