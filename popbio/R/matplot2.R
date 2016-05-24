matplot2<-function(x, proportions=FALSE, legend="topright",  xlab=NULL, ylab=NULL, type='l', 
                   las=1, pch=c(15:18, 1:3), lwd=1, lty=1:nrow(x), col=rainbow(nrow(x)),
                   lcex=1, lbty='o', lcol=1, ltitle=NULL, lsort=TRUE, ... )
{
   n<-nrow(x)
   if (is.null(n)) {stop("x should be a matrix with at least 2 rows")}   
   if (proportions){
       x<-prop.table(x, 2)      ## Change counts to proportions using prop.table
       if(is.null(ylab)){ ylab<-"Proportion" }
   }
   if (length(col)<n) {col<-rep(col,length.out=n)}  ## line colors (repeat if necessary)
   if (length(lty)<n) {lty<-rep(lty,length.out=n)}  ## line types
   if (length(pch)<n) {pch<-rep(pch,length.out=n)}  ## point types
   if (is.null(ylab)) { ylab<-"Total"}                 ## default y label
   if (is.null(colnames(x))) { colnames(x)<-1:ncol(x)} ## No column names, use default 1:n
   ## convert colnames to number if possible
   suppressWarnings(xnames<-as.numeric( colnames(x)  )) 
   # any characters in column names?
   if (any(is.na(xnames))) {
      if (is.null(xlab)){xlab<-"Stage"}
      matplot(t(x), xlab=xlab, ylab=ylab, col=col, las=las, type=type, lty=lty, lwd=lwd, pch=pch, xaxt='n', ...)
      ## probably should add pretty breakpoints if ncol(x)>10
      axis(1, at=1:ncol(x), labels=colnames(x))
   }
   ## else numeric column names
   else {    
      if (is.null(xlab)){xlab<-"Year"}
      matplot(xnames, t(x), xlab=xlab, ylab=ylab, col=col, las=las, type=type, lty=lty, lwd=lwd, pch=pch, xaxt='n', ...)
      ## if only a few years, do not label  1992.0, 1992.5, etc
      if (ncol(x)<6){
        axis(1, at=xnames, labels=as.character(xnames))
      }
      else {
         axis(1, at=pretty(xnames))
      }
   } 
   ## order legend by decreasing order of mean number in row
   if (lsort) {
       y<-sort(apply(x,1,mean, na.rm=TRUE), index.return=TRUE, decreasing=TRUE)
   }
   else {
      y<-list(x=apply(x,1,mean, na.rm=TRUE), ix=1:n)
   }
   leg.names<-paste(names(y$x), "")  ## pad names with trailing space for extra room 
   ## if rownames are missing
   if (leg.names[1] == " ") { leg.names <- paste("row", y$ix, "")}

   # fix to pass args.legend option like barplot.default 
   ## legend box is smaller if pch is used (fits to point)
   if (type=='l'){pch<-NULL}
   legend(legend[1],legend[2], leg.names, lty=lty[y$ix],  col=col[y$ix],
               lwd=lwd, pch=pch[y$ix], cex=lcex, bty=lbty, ncol=lcol, title=ltitle)
}
