draw.barplot2D <-
function(x,y,z,width,height,scale=F,col=NULL,col.frame='black',lwd.frame=1,silent=TRUE,...){
  nx <- length(x)
  nz <- dim(z)[2]
  if (is.null(col)) 
      col <- c("#737373", "#F15A60", "#7BC36A", "#599BD3", "#F9A75B", "#9E67AB", "#CE7058", "#D77FB4")
  col <- rep(col, length.out = nz)
  if(length(y)!=nx) stop('x and y should be vectors of the same length')
  if(length(dim(z))!=2) stop('z should be a 2-dimensional array')
  if(dim(z)[1]!=nx) stop('the number of rows in of z should match as the length of x and y')
  if(length(width)!=length(height)) stop('width and height should have the same length')
  if(length(width)>1 & length(width)!=length(x)) stop('width and height should have the same length as x and y')  
  maxsumz <- max(rowSums(z,na.rm=T),na.rm=T)
  pm <- setProgressMsg(1,nx)
  for(i in 1:nx){
    xi=x[i]
    yi=y[i]
    zi=z[i,]
    if(length(width)>1) widthi <- width[i] else widthi <- width
    if(length(height)>1) heighti <- height[i] else heighti <- height
    if(scale & length(width)==1) {
      widthi <- width * sqrt(sum(zi,na.rm=T))/sqrt(maxsumz)
      heighti <- height * sqrt(sum(zi,na.rm=T))/sqrt(maxsumz)    
    }
    j=which(zi>0)
    if(sum(zi,na.rm=T)>0) barplot2D(z=zi[j],colour=col[j],x=xi,y=yi,width=widthi,height=heighti,add=T,col.frame,lwd.frame,...)
    if(!silent) pm <- progressMsg(pm,i)
    }
}

