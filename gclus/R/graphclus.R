														
lower2upper.tri.inds <-
#copied from lower.to.upper.tri.inds from cluster library
function (n) 
{
    n1 <- as.integer(n - 1)
    if (n1 < 1) 
        stop("`n' must be >= 2")
    else if (n1 == 1) 
        1:1
    else rep(1:n1, 1:n1) + c(0, unlist(lapply(2:n1, function(k) cumsum(c(0, 
        (n - 2):(n - k))))))
}

vec2distm <- function(vec){
    #convert from a vector to a distance matrix
    m <- length(vec)
    n <- (1+sqrt(1+8*m))/2
    ans<- matrix(0,n,n)
    ans[lower.tri(ans)] <- vec
    ans[upper.tri(ans)] <- vec[lower2upper.tri.inds(n)]
    ans
}



vec2dist <- function(vec){
#convert from a vector to a "dis"
     as.dist(vec2distm(vec))
}

# Returns a vector of off-diagonal elements in m.
# The off parameter specifies the distance above the main (0) diagonal.

diag.off <- function(m,off=1)
	m[col(m)==row(m)+off]

#-----------------------------------------------------------

# Accepts a dissimilarity matrix or "dist" m, and
# returns a  matrix of colors.
# M values are cut into categories using breaks (ranked distances if 
# byrank is true) and categories  are assigned the values in colors.

default.dmat.color <- c("#FDFFDA", "#D2F4F2", "#F4BBDD")

dmat.color <-
function(m, colors = default.dmat.color,byrank=NULL, breaks=length(colors) ){ 
   if (is.matrix(m)) m <- as.dist(m)
   if (is.null(byrank))
   byrank <- length(breaks) == 1
   if (byrank ==TRUE)
       m1 <- rank(as.vector(m))
   else
       m1 <- as.vector(m)
   fac <- cut(m1,breaks,include.lowest=TRUE)
   ans <- colors[as.numeric(fac)]
   ans <- vec2distm(ans)
   diag(ans) <- NA
   attr(ans,"Levels") <- levels(fac)
   if (length(labels(m)) == nrow(ans)){
       rownames(ans) <- labels(m)
       colnames(ans) <- labels(m)}
   ans
	   
}


#-----------------------------------------------------------
#

# Extracts information from a matrix of colors suitable for use by
# image.
# 
imageinfo <- function(cmat) {
    n <- nrow(cmat)    
    p <- ncol(cmat) 
    levels <- sort(unique(as.vector(cmat)))
    z <- unclass(factor(cmat,levels= levels, labels=1:length(levels)))
    z <- matrix(z,nrow=n,p)
    list(x=1:p,y=1:n, z =t(z),col=levels)
}
 

# This draws the color matrix cmat.


plotcolors <- function(cmat,  na.color="white", dlabels = NULL, rlabels = FALSE, clabels = FALSE, 
    ptype ="image", border.color = "grey70", pch=15,cex=3,label.cex = .6,...) {
        
    n <- nrow(cmat)    
    p <- ncol(cmat) 
    cmat[is.na(cmat)] <- na.color
    if (ptype=="image") {
	info <- imageinfo(cmat)
	image(info$x, info$y, info$z[, n:1], col = info$col, 
	    axes = FALSE, xlab = "", ylab = "", ...)}
    else {
	y <- rep(n:1,p)
	x <- rep(1:p,rep(n,p)) 
	cmat <- as.vector(cmat)
	plot(x,y,col=cmat,cex=cex,pch=pch,axes=FALSE,xlab="",ylab="",
	    xlim=c(.5,p+.5),ylim=c(.5,n+.5),...)
	
    }
    axis(3, at = 1:p, tick=FALSE,labels = clabels, 
	las = 2, cex.axis = label.cex)
    axis(2, at = n:1, tick=FALSE,labels = rlabels, 
	las = 2, cex.axis =label.cex)
    if (is.vector(dlabels)){
	nl <- length(dlabels)
	text(1:nl,nl:1,dlabels,cex=label.cex)}
    box(col = border.color)
}

    




#-----------------------------------------------------------
# This function draws a scatterplot matrix of data.
# Order, if present, specifies the order of the variables and
# panel.colors, if present should be a matrix of panel colors.
# (...) are graphical parameters.

cpairs <-
function(data,order=NULL,panel.colors=NULL,border.color="grey70",show.points=TRUE,...) {
    textPanelbg <- function(x = 0.5, y = 0.5, txt, cex, font) {
	box(col= border.color)
	text(x, y, txt, cex = cex, font = font)
    }
    
    if (!is.null(order)) {
	data <- data[,order]
	if (!(is.null(panel.colors)))
	   panel.colors <- panel.colors[order,order]}
    
    if (!is.null(panel.colors)) {
	if (ncol(data) != nrow(panel.colors) || ncol(data) != ncol(panel.colors))
	   stop("dimensions do not match")	
	diag(panel.colors) <- NA
	panel.colors <- t(panel.colors)[!is.na(panel.colors)]}
    
    env<- new.env()
    assign("j",1,envir=env) 
    pairs.default(data,...,text.panel = textPanelbg,
	panel = function(x,y,...){
	    j <- get("j",envir=env)
	    reg <- par("usr")
	    if (!(is.null(panel.colors)))
	       rect(reg[1],reg[3],reg[2],reg[4],col=panel.colors[j])
	    box(col=border.color)
	    j <- j+1 
	    assign("j",j,envir=env)
	    if (show.points == TRUE) points(x,y,...)	
	})
    
}




# This function draws a parallel coordinate plot  of the data.
# Order, if present, specifies the order of the variables and
# panel.colors, if present should either be a vector of panel colors,
# or a matrix whose i,j the element gives the color for the panel
# showing columns i and j of data. (...) are graphical parameters.
# This function is adapted from parcoord(MASS).


cparcoord <-
function (data, order=NULL,panel.colors=NULL,col=1,lty=1,horizontal=FALSE,mar=NULL,...) {
    if (is.null(mar))
       if (horizontal==TRUE)
          mar <- c(5, 2, 2, 2) + 0.1
        else mar <- c(2, 8, 2, 2) + 0.1
    if (!is.null(order)) {
	data <- data[,order]
	if (is.matrix(panel.colors))
	   panel.colors <- panel.colors[order,order]}
    
    if (is.matrix(panel.colors))
        panel.colors <- diag.off(panel.colors)
      
    if (is.vector(panel.colors))
       if (ncol(data) -1 != length(panel.colors))
          stop("dimensions do not match")
      
    oldpar <- par(mar=mar)
    x <- apply(data, 2, function(x) (x - min(x))/(max(x) - min(x)))
    p <- ncol(x)
    if (horizontal==TRUE){
	matplot(1:p, t(x), 
	    xlab = "", ylab = "", axes = FALSE, type="n",...)
	axis(1, at = 1:p, labels = colnames(x))
	if (!(is.null(panel.colors)))
	for (i in 1:(p-1)) rect(i,0,i+1,1,  lty=0,col =panel.colors[i])
	for (i in 1:p) lines(c(i, i), c(0, 1), col = "grey70")
	matpoints(1:p, t(x), type = "l",col=col,lty = lty,...)  
    }
    else {  
	matplot(t(x), p:1, 
	    xlab = "", ylab = "", axes = FALSE, type="n",...)
	axis(2, at = p:1, labels = colnames(x),las=2)
	if (!(is.null(panel.colors)))
	for (i in 1:(p-1)) rect(0,i,1,i+1,  lty=0,col =panel.colors[p-i])
	for (i in 1:p) lines(c(0, 1),c(i, i), col = "grey70")
	matpoints(t(x), p:1, type = "l",col=col,lty = lty,...)  
    }
    on.exit(par(oldpar))
    invisible()
}


