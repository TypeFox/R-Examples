scaleBathy = function (mat, deg=1, x="bottomleft", y=NULL, inset=10, angle=90, ...) {
			
	usr=par("usr")
	
	if(is.numeric(x) == TRUE & is.null(y)==FALSE) {
		x->X
		y->Y		
		abs(Y) -> lat
	}
	
	if(is.numeric(x) == FALSE & is.null(y)==TRUE) {
		insetx = abs((usr[2]-usr[1])*inset/100)
		insety = abs((usr[4]-usr[3])*inset/100)
	    X <- switch(x,  bottomright = (usr[2]-insetx-deg), 
	    				topright = (usr[2]-insetx-deg), 
	    				bottomleft = (usr[1]+insetx),
	    				topleft = (usr[1]+insetx) )
	    Y <- switch(x,  bottomright = (usr[3]+insety), 
	    				topright = (usr[4]-insety),
	    				bottomleft = (usr[3]+insety),
	    				topleft = (usr[4]-insety) )
	    lat<- switch(x, bottomright = abs(min(as.numeric(colnames(mat)))), 
	    				topright =    abs(max(as.numeric(colnames(mat)))), 
	    				bottomleft =  abs(min(as.numeric(colnames(mat)))),
	    				topleft =     abs(max(as.numeric(colnames(mat)))) )
	}
     
	cos.lat <- cos((2 * pi * lat)/360)
	perdeg <- (2 * pi * (6372.798 + 21.38 * cos.lat) * cos.lat)/360
          	
	arrows(X, Y, X+(deg),Y, code=3, length=0.05, angle=angle)
	text((X + X+(deg))/2, Y, 
			adj=c(0.5,-0.5),
			labels=paste(round(perdeg*deg,0),"km"), ...)	
	
}