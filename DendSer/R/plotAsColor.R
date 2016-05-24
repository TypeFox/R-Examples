plotAsColor <-
function(m,order=NULL,order.col=order,order.row=order,rank=FALSE, 
    border.color = "grey70", labels=FALSE,x=1:ncol(d),y=1:nrow(d),...)	{
    d <- as.matrix(m) 
	if (rank)  {
		dr <- rank(d)
		dr <- matrix(dr, nrow=nrow(d))
		} 
		else dr <- d    
	 n <- nrow(dr)
     p <- ncol(dr)
    rownames(dr)<- rownames(d)
    colnames(dr)<- colnames(d)
	 if (is.null(order.row) | length(order.row) !=n) order.row <- 1:n
	 if (is.null(order.col) | length(order.col) !=p) order.col <- 1:p
	 
    
	 dr <- t(dr[rev(order.row),order.col])
	  image(x,y,dr,axes = FALSE, xlab = "", ylab = "",...)
     if (! is.null (border.color)) box(col = border.color)
     if (labels){
     if (! is.null(rownames(dr)))
       axis(3, at = x, labels = rownames(dr))
     if ( ! is.null(colnames(dr)))
       axis(2, at = y, labels = colnames(dr),las=2)
       }
}



