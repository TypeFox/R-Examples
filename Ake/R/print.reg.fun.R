print.reg.fun <-
function(x, digits=NULL, ...)
           {
   	     class(x) <- "reg.fun"
		  if(x$kernel=="GA") 		kerf="gamma"
	 	  else if(x$kernel=="BE") 	kerf= " extended beta "
	    	  else if(x$kernel=="LN") 	kerf= " lognormal"
	     	  else if(x$kernel=="RIG")	kerf= " reciprocal inverse Gaussian"
		  else if(x$kernel=="bino") 	kerf=" Binomial"
	 	  else if(x$kernel=="triang") 	kerf= " Triangular "
	    	  else if(x$kernel=="dirDU") 	kerf= "DiracDU"
	cat("\nBandwidth h:",formatC(x$h,digits=digits), "\tCoef_det = ",x$Coef_det,"\n",
	   	 "\nNumber of points: ",x$n,";","\tKernel = ",kerf, "\n\n",sep="")
        print(summary(as.data.frame(x[c("data","y")])), digits=digits, ...)
        print(summary(as.data.frame(x[c("eval.points","m_n")])), digits=digits, ...)
	 #invisible(x)
	}
