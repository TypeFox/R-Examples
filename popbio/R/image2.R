image2<-function(x, col=c('white', rev(heat.colors(23))),  breaks, log=TRUE,
                 border=NA, box.offset=0.1, round=3, cex, text.cex=1, text.col="black",
                 mar=c(1,3,3,1),  labels=2:3, label.offset=0.1, label.cex=1, srt=90 )
{
   ## convert vector like 1:5 to matrix with 1 row (default is one column) 
   if(!is.matrix(x)){
     x <- t(as.matrix(x))  # only transpose if vector
   }  
   if(!is.numeric(x)){
     stop("A numeric matrix is required")
   }
   if(!missing(cex)){
     text.cex=cex; label.cex=cex
   }  ## cex replaces any values in text.cex or label.cex
   op <- par(mar=mar, xpd=TRUE)
   x  <- x[nrow(x):1, ,  drop=FALSE]       ## flip matrix so top row on botton   
   x1 <- ncol(x)                          ## number of columns and rows
   y1 <- nrow(x)
   # hack to get three colors needed for cut
   if(length(col)==1){col<-rep(col, 3)}
   if(length(col)==2){col<-c(col, col[-1])}
   if(missing(breaks)){  breaks<-length(col) - 1 }   ## number of breaks
   else{ if(length(breaks) != length(col)) { warning("Breaks is not the same length as colors.
Some blocks may be unfilled or some colors may not be used")}}

   ## check if any values < 0.  If so, minimum value will be in first color category.
   checkmin<-min(x, na.rm=TRUE)
   if(checkmin<0){x <- x - checkmin}
   ## first check for NAs
   missingNA<-is.na(x)
   x[x==0]<-NA                          ## set zeros to NA for log10  
   if(log){z<-cut(log10(x), breaks)}    ## cut into intervals using log10 transformation
   else{z<-cut(x, breaks)}
   z2<-matrix(z, y1, x1)                ## reshape into matrix
   x[is.na(x)]<-0     ## set NA values back to zero 
   x[missingNA]<-NA   ## reset NAs 
   
   if(checkmin<0){x<-x + checkmin}      
   x<-round(x, round)                   ## round decimal places

   ## PLOT rectangles
   offset<-box.offset/2
   plot(seq(1, x1+1, len=2), seq(1, y1+1, len=2), type='n', axes=FALSE, ann=FALSE, xaxs='i', yaxs='i')
   for(i in 1:x1)
   {
      for(j in 1:y1)
      {
        if(is.na(z2[j,i])){n1<-1}       ## if element is zero, use first color (default white)
         else{ n1<-match(z2[j,i], levels(z)) + 1}  # else match to cut levels in z2 and use that color (plus 1)   
        rect(i+offset, j+offset, i+1-offset, j+1-offset, border=border, col=col[n1])
        if(!is.na(text.col)){ text(i+.5,j+.5, x[j,i], cex=text.cex, col=text.col)}
      }
   }
   ## rownames (left and right) 
   if(2 %in% labels) text(1-label.offset,    1:y1+.5, rownames(x), pos=2, offset=0, cex=label.cex)  ## on left
   if(4 %in% labels) text(x1+1+label.offset, 1:y1+.5, rownames(x), pos=4, offset=0, cex=label.cex)
   ## colnames (bottom and top)
   if(1 %in% labels) text(1:x1+.5, 0.9-label.offset,    colnames(x), pos=2, offset=0, srt=(0+srt), cex=label.cex)
   if(3 %in% labels) text(1:x1+.5, y1+1.1+label.offset, colnames(x), pos=2, offset=0, srt=(360-srt), cex=label.cex)
   par(op)
 }
