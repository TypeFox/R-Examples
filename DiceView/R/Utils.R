DiceView.env <- new.env()


##=========================================================
## make colors semi-transparent
##
## alpha = 0 OPAQUE
## alpha = 1 FULLY TRANSPARENT
##
##=========================================================

translude <- function(colors, alpha = 0.6) {
    
    alpha <- rep(alpha, length.out = length(colors))
    rgb <- as.matrix(col2rgb(colors)/255)
    colors2 <- rgb(red = rgb["red", ],
                   green = rgb["green", ],
                   blue = rgb["blue", ],
                   alpha = alpha)
    
    
}

##========================================================
## level color: one base color incermented in hsv to provide a palette
## 
##========================================================

col.levels <- function(color,nlevels){
        col.rgb=col2rgb(color)
        col.hsv=rgb2hsv(r=col.rgb[1],g=col.rgb[2],b=col.rgb[3])
        col = hsv(h=col.hsv[1],s=seq(f=0,t=col.hsv[2],l=nlevels),v=col.hsv[3])
        return(col)
}

##========================================================
## fade color: one color and several fading values
## between 0 and 1
##
## alpha ~ 0 nearly equal to the input 'color'
## alpha ~ 1 nearly invisible grayed/transparent 
##
## The returned value is vector. It is build
## from a matrix with four rows as
## those produced by col2rgb(x, alpha = TRUE)
## 
##========================================================

fade <- function(color = "red",
                 alpha =  seq(from = 0, to = 1, length.out = 5),
                 plot = FALSE) {
    
    if (any(alpha < 0) || any(alpha > 1)) stop("'alpha' values must be >=0 and <= 1")
    if (length(color) > 1) stop("'color' must be of length 1")  
    
    ## a matrix with  1 col
    rgbcol <- col2rgb(color)/255
    mat <- matrix(1-alpha, nrow = 3, ncol = length(alpha), byrow = TRUE)
    
    mat <- mat + rgbcol %*% alpha 
    
    colors2 <- rgb(red = mat[1, ],
                   green = mat[2, ],
                   blue = mat[3, ],
                   alpha = alpha)
    
    if (plot) {
        x <- seq(from = 0, to = 1, length.out = length(alpha))
        plot.new( )
        for ( i in 1:(length(alpha)) ){
            rect(xleft = x[i],
                 xright = x[i+1],
                 ybottom = 0,
                 ytop =  1,
                 border = NA,
                 col = colors2[i])
        }
        
    }
    
    colors2
    
}

##========================================================
## try to find a good formatted value for a numeric vector
## x using a vector of diff range drx
##
## For instance, if drx is 1000, no decimal or very few
## decimals should be used.
## 
##========================================================

tryFormat <- function(x, drx) {
    
    d <- length(x)
    ldx <- log(drx, base = 10)
    ff <- rep(1, d)
    fd <- rep(1, d)
    ff[ldx > 0] <- ceiling(ldx[ldx > 0])
    fd[ldx < 0] <- ceiling(-ldx[ldx < 0]) + 2
    ff <- ff + fd +1
    
    formats <- paste("%", ff, ".", fd, "f", sep = "")
    fx <- sprintf(formats, x)
    
}