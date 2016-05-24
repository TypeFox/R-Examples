## The functions imagematrix', 'rgb2grey', 'clipping' and 'read.jpeg' are adopted from the ReadImages package, which has been orphaned and archived by CRAN.

graph.extract<-function(MT, refX, refY, save="no", image=read.jpeg(file.choose())){

imagematrix <- function(mat, type=NULL, ncol=dim(mat)[1], nrow=dim(mat)[2],
                        noclipping=FALSE) {
  if (is.null(dim(mat)) && is.null(type)) stop("Type should be specified.")
  if (length(dim(mat)) == 2 && is.null(type)) type <- "grey"
  if (length(dim(mat)) == 3 && is.null(type)) type <- "rgb"
  if (type != "rgb" && type != "grey") stop("Type is incorrect.")
  if (is.null(ncol) || is.null(nrow)) stop("Dimension is uncertain.")
  imgdim <- c(ncol, nrow, if (type == "rgb") 3 else NULL)
  if (length(imgdim) == 3 && type == "grey") {
    # force to convert grey image
    mat <- rgb2grey(mat)
  }
  if (noclipping == FALSE && ((min(mat) < 0) || (1 < max(mat)))) {
    warning("Pixel values were automatically clipped because of range over.") 
    mat <- clipping(mat)
  }
  mat <- array(mat, dim=imgdim)
  attr(mat, "type") <- type
  class(mat) <- c("imagematrix", class(mat))
  mat
}

rgb2grey <- function(img, coefs=c(0.30, 0.59, 0.11)) {
  if (is.null(dim(img))) stop("image matrix isn't correct.")
  if (length(dim(img))<3) stop("image matrix isn't rgb image.")
  imagematrix(coefs[1] * img[,,1] + coefs[2] * img[,,2] + coefs[3] * img[,,3],
              type="grey")
}
clipping <- function(img, low=0, high=1) {
  img[img < low] <- low
  img[img > high] <- high
  img
}
  
read.jpeg <- function(filename) {
  res <- .C("get_imagesize_of_JPEG_file", as.character(filename),
            width=integer(1), height=integer(1), depth=integer(1),
            ret=integer(1), PACKAGE="SCVA")
  if (res$ret < 0)
    stop(if (res$ret==-1) "Can't open file." else "Internal error")
  imgtype <- if (res$depth == 1) "grey" else "rgb"
  imgdim <- c(res$height, res$width, if (res$depth == 3) res$depth else NULL)
  res <- .C("read_JPEG_file", as.character(filename),
            image=double(res$width * res$height * res$depth),
            ret=integer(1), PACKAGE="SCVA")
  img <- array(res$image, dim=imgdim)
  imagematrix(img/255, type=imgtype)
}
 
  plot(image)
  refpoints <- locator(n = 4, type = 'p', pch = 4, col = 'blue', lwd = 2)
  refpoints <- as.data.frame(refpoints)	
  datapoints <- locator(n=MT,type='p',pch=1,col='red',lwd=2,cex=2)
  datapoints <- as.data.frame(datapoints)	
  x <- refpoints$x[c(1,2)]
  y <- refpoints$y[c(3,4)]
  cx <- lm(formula=c(refX[1],refX[2])~c(x))$coeff
  cy <- lm(formula=c(refY[1],refY[2])~c(y))$coeff
  datapoints$x <- datapoints$x*cx[2]+cx[1]
  datapoints$y <- datapoints$y*cy[2]+cy[1]
  true.data <- as.data.frame(datapoints)
  plot(true.data,type='b',pch=1,col='blue',lwd=1.1,bty='l')
  rounded <- round(true.data,digits=2)
  if(save=="yes"){
    write.table(rounded,file=file.choose(new=FALSE),col.names=FALSE,row.names=FALSE,append=FALSE,sep="\t")
  }
  return(rounded)

}

