print.indicators <- function (x, At = 0, Bt = 0, sqrtIVt = 0, selection = NULL, confint=FALSE,...) {
	if(is.null(selection)) selection = rep(TRUE, nrow(x$C))
	if(length(dim(x$A))==2) {
		A = x$A[selection,]
		B = x$B[selection,]
		sqrtIV = x$sqrtIV[selection,]
	} else {
		A = x$A[selection]
		B = x$B[selection]
		sqrtIV = x$sqrtIV[selection]
	}
	C = subset(x$C,selection)
	spnames = names(C)	
    if(length(dim(A))==2) {
    	sel = A$stat>=At & B$stat>=Bt & sqrtIV$stat>=sqrtIVt 
    	if(confint) {
    		nc= 9
	    	A = A[sel,]
    		B = B[sel,]
    		sqrtIV = sqrtIV[sel,]
    	} else {
	    	nc= 3	
	    	A = A$stat[sel]
    		B = B$stat[sel]
    		sqrtIV = sqrtIV$stat[sel]
    	}
    } else {
    	sel = A>=At & B>=Bt & sqrtIV >=sqrtIVt
    	A = A[sel]
    	B = B[sel]
    	sqrtIV = sqrtIV[sel]
    	nc = 3
    }
    sel[is.na(sel)]=FALSE
    CM = subset(C,sel)
    m = data.frame(matrix(0,nrow = sum(sel), ncol=nc))
    if(nc==3) names(m) = c("A","B","sqrtIV")
    else if(nc==9) names(m) = c("A","LA","UA","B","LB","UB","sqrtIV", "LsqrtIV","UsqrtIV")
    if(sum(sel)>0) {
    	if(nc==3) {
    		m[,1] = A
    		m[,2] = B
    		m[,3] = sqrtIV
    	} else if (nc==9){
    		m[,1:3] = A
    		m[,4:6] = B
    		m[,7:9] = sqrtIV
	    }
    	for(r in 1:sum(sel)) {
    		row.names(m)[r] <- paste(spnames[CM[r,]==1], collapse = "+")
    	}
    }
    m = m[order(m$sqrtIV,decreasing=TRUE),]
	print(m,...)
}
