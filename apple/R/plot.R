#plot function
#x: the path object
#add: add or not
#col: color

plot.apple=function(x,col='black',add=FALSE,main='apple', type='l', lty=1, ...){
	
	class=class(x)[1]
	
	if(class=='apple'){
		matplot(x$lambda,t(x$beta),type=type, col=col, add=add,xlab=expression(lambda), ylab=expression(hat(beta)), main=main, xlim=rev(range(x$lambda)), lty=lty)
		min=min(x$beta)
		max=max(x$beta)
		
		#add the cv selected line	
		if(is.null(x$cv.loc)==FALSE){
		    lam.cv=x$lambda[x$cv.loc]
            lines(x=rep(lam.cv, 1000), y=seq(from=min,to=max, length.out=1000), col=col)
  
            #add the ebic selected line
            lam.ebic=x$lambda[x$ebic.loc]
            lines(x=rep(lam.ebic,1000), y=seq(from=min, to=max, length.out=1000), col=col,lty=2)
        	}

		    }		              
        
    }
