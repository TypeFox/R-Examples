##' extended bar plot
##'
##' @param x data
##' @param names names
##' @param names_height height of names
##' @param col color
##' @param ... passed on
##' @return NULL
##' @export
superbarplot <- function(x,
                         names = 1:dim(x)[2],
                         names_height=NULL,
                         col = gray(seq(.8,.5,length=dim(x)[1]/2)), ...) {

  plot.bar <- function(x,min,max,width=1,...) {
    alpha <- (1-width)/2
    polygon(x + c(alpha,alpha,1-alpha,1-alpha,alpha),
            c(min,max,max,min,min),
            ...)                        # pass in col
  }
    
  ## x is a matrix with rows alternating of High, Min.
  n = dim(x)[2]
  m = dim(x)[1]
  no.bars = dim(x)[1]/2
  
  y.range = c(min(x),max(x))
  x.range = c(1,n+1)

  
  ## setup plot
  plot.new()
  plot.window(xlim=x.range,ylim=y.range,
              xaxt="n",
              bty="n",ann=FALSE)
  title(...)
  
  for(i in 1:no.bars) {
    for(j in 1:n) {
      plot.bar(j,x[2*i-1,j],x[2*i,j],width=1 - i/(3*no.bars),col=col[i])
    }
  }

  ## names
  if(!is.null(names)) {
    ## height
    if(is.null(names_height)) {
      f = par("yaxp")
      names_height= f[1] + (f[2]-f[1])*(f[3]-1)/f[3]
    }
    text(0.5 + 1:n,rep(names_height,n),format(names))
  }
    

}

  
