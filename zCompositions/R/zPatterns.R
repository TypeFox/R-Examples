zPatterns <- function(X,label=NULL,plot=TRUE,
                      axis.labels=c("Component","Pattern number"),
                      bar.colors=c("red3","red3"), bar.labels=FALSE,
                      cell.colors=c("dodgerblue","white"), cell.labels=c(label,paste("No",label)),
                      cex.axis=1.1, grid.color="black",
                      grid.lty="dotted", legend=TRUE, ...){
  
  plot.patterns <- function(a, ...)
  {
    
    zones <- matrix(c(2,4,1,3), ncol=2, byrow=TRUE)
    layout(zones, widths=c(4/5,1.5/5), heights=c(2/5,3.5/5))
    par(mar=c(3,3,0.5,0.5))
    a <- as.matrix(a[rev(rownames(a)),])
    image(1:ncol(a),1:nrow(a),t(a),col=rev(cell.colors),axes=F)
    mtext(side=1,text=axis.labels[1],line=1.75)
    mtext(side=2,text=axis.labels[2],line=1.75)
    par(mgp=c(3, .3, 0))
    axis(side = 1,at = seq(1,ncol(a),by=1),labels=colnames(a),tck=0,cex.axis=cex.axis)
    axis(side = 2,at = seq(1,nrow(a),by=1),labels=rownames(a),las=2,tck=0,cex.axis=cex.axis)
    box()
    grid(ncol(a),nrow(a),col=grid.color,lty=grid.lty)
    par(mar=c(0,3.25,1,0.75))
    a <- barplot(as.vector(prop_col),axes=F,col=bar.colors[1],xaxs="i",
                 ylim=c(0,max(as.vector(prop_col)+0.2*max(as.vector(prop_col)))))
    if (bar.labels==TRUE)
      text(a,as.vector(prop_col),labels=as.vector(prop_col),cex=0.85,pos=3)
    par(mar=c(3.25,0,0.75,0.75))
    a <- barplot(rev(as.vector(pat.freq.perc)),horiz=T,axes=F,col=bar.colors[2],yaxs="i",
                 xlim=c(0,max(as.vector(pat.freq.perc)+0.3*max(as.vector(pat.freq.perc)))))
    if (bar.labels==TRUE)
      text(rev(as.vector(pat.freq.perc)),a,labels=rev(as.vector(pat.freq.perc)),cex=0.85,pos=4)
    par(mar=c(0,0,3,0))
    plot.new()
    if (legend==TRUE){
    legend("topleft",cell.labels,pch=c(22,22),bty="n",
           pt.bg=cell.colors,pt.cex=2,cex=1.1)}
  } 
  
  if (is.vector(X)) stop("X must be a matrix or data.frame class object")
  if (is.null(label)) stop("A value for label must be given")
  if (!is.na(label)){
    if (!any(X==label,na.rm=T)) stop(paste("Label",label,"was not found in the data set"))
    if (label!=0 & any(X==0,na.rm=T)) stop("Zero values not labelled as censored values were found in the data set")
    if (any(is.na(X))) stop(paste("NA values not labelled as censored values were found in the data set"))
  }
  if (is.na(label)){
    if (any(X==0,na.rm=T)) stop("Zero values not labelled as censored values were found in the data set")
    if (!any(is.na(X),na.rm=T)) stop(paste("Label",label,"was not found in the data set"))
  }
  
  X <- as.data.frame(X)
  
  n <- nrow(X); p <- ncol(X)
  X[X==label] <- NA
  

  miss <- as.data.frame(is.na(X)*1)

  prop_col <- round(colSums(miss)/n*100,2)
  prop <- round(sum(miss)/(n*p)*100,2)

  miss <- cbind(miss,pat=do.call(paste,c(miss,sep="")))
  tab <- miss[!duplicated(miss),]
  tab.num <- tab
  tab[tab==1] <- "+"
  tab[tab==0] <- "-"
  pat.freq <- sort(table(miss$pat),decreasing=T)
  pat.freq.perc <- round(as.vector(pat.freq)/sum(pat.freq)*100,2)
  pat.ID <- miss$pat; levels(pat.ID) <- match(levels(miss$pat),names(pat.freq))
  pat.ID <- factor(pat.ID,levels=1:length(levels(pat.ID)))
  
  # Summary
  
  tab <- tab[match(names(pat.freq),tab$pat),1:p]
  tab.num <- tab.num[match(names(pat.freq),tab.num$pat),1:p]
  rownames(tab.num) <- 1:nrow(tab.num)
  if (plot==TRUE) plot.patterns(tab.num)
  tab <- cbind(Patt.ID=1:length(levels(pat.ID)),
               tab,
               No.Unobs=rowSums(tab[,1:p]=="+"),
               Patt.Freq=as.vector(pat.freq),
               Percent=pat.freq.perc)

  cat("Patterns ('+' means ",cell.labels[1],", '-' means ",cell.labels[2],") \n\n",sep="")
  print(tab,row.names=FALSE)
  cat("\n")
  cat("Percentage cases by component \n")
  print(prop_col)
  cat("\n")
  cat(paste("Overall percentage cases: ",prop))
  
  invisible(pat.ID)
}