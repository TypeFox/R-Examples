summMCMC <-
function(xList,EPSfileName,PDFfileName,plotInd=1,parNames,
                     columnHeadings,colourVersion=TRUE,columnCols,
                     credLevel=0.95,numerSummCex=1.3,BGRsttPos=10,
                     BGRyRange=c(0.95,1.25),BGRtickPos=1.2,
                     BGRlogTransf,BGRlogitTransf,KDExlim,KDEvertLine=TRUE,
                     KDEvertLineCol="black",addTruthToKDE=NULL)
{
   
   options(warn=-1)
   x <- lapply(xList,as.matrix)
   num.par <- ncol(x[[plotInd]])
   samp.size <- nrow(x[[plotInd]])
   num.chains <- length(x)
   
   
   # Divert figure to PDF file if filename specified.

   if (!missing(PDFfileName))
       pdf(PDFfileName,width=9,height=(num.par+1))

   if (!missing(EPSfileName))
       postscript(EPSfileName,horizontal=FALSE,width=11,height=(num.par+1))

   op <- par()
   if (missing(columnHeadings))
      columnHeadings <- c("parameter","trace","lag 1","acf","BGR",
                           "density","summary")   

   if (colourVersion)
   {

      if (missing(columnCols))
      {
         columnCols <- c("darkmagenta", "green4","darkorange","dodgerblue",
                          "darkgoldenrod1","red","navy")
         if (KDEvertLine&missing(KDEvertLineCol))
            KDEvertLineCol <- "DarkGreen"           
      }
   }

   panels.per.par <- length(columnHeadings)
   if (num.chains==1) panels.per.par <- panels.per.par - 1

   if (!colourVersion)
   {
      columnCols <- rep(NA,num.par)
      for (i in 1:length(columnHeadings))
         columnCols[i] <- "black"
   }
 
   empty.panel <- function()
   {
    plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), xaxt = "n",
        yaxt = "n", xlab = "", ylab = "", bty ="o")
    invisible()
   }

    par(mfrow=c((num.par+1),panels.per.par))

   headInds <- 1:7
   if (num.chains==1) headInds <- c(1:4,6,7)
   for (i in 1:panels.per.par)
   {
      par(ann=F,mar=c(0,0,0,0),xaxt="n",yaxt="n")
      empty.panel()
      text(0.5,0.5,columnHeadings[headInds[i]],cex=2.2,col=columnCols[headInds[i]])
   }

   for (j in 1:num.par)
   {  
      if (j %% 15 == 0) {dev.new()}
 
      # Write the variable name

      par(ann=F,mar=c(0,0,0,0),xaxt="n",yaxt="n")

      empty.panel()
      
         text(0.5,0.5,parNames[[j]],cex=2.2,col=columnCols[1])
      

      # Do the trace plot

      plot(x[[plotInd]][,j],xlab="",ylab="",type="l",col=columnCols[2])
 
      # Do the lag 1 plot

      plot(x[[plotInd]][1:(samp.size-1),j],
           x[[plotInd]][2:samp.size,j],xlab="",ylab="",type="n")

      points(x[[plotInd]][1:(samp.size-1),j],
             x[[plotInd]][2:samp.size,j],pch=1,cex=0.5,col=columnCols[3])
     
      # Do the ACF plot

      ci.col.val <- "black"
      if (colourVersion) ci.col.val <- "blue"
      acf(x[[plotInd]][,j],lag.max=20,col=columnCols[4],
          lwd=2,ci.col=ci.col.val)

      # Do the density plot

      h <- dpik(x[[plotInd]][,j])
      est <- bkde(x[[plotInd]][,j],bandwidth=h)

      if (!missing(KDExlim))
         xlim.val <-  KDExlim[[j]]   

      if (missing(KDExlim))
         xlim.val <-  range(est$x)

      lb.ht <- 0.2*(max(est$y)-min(est$y))
      ylim.val <- c(-lb.ht,1.1*max(est$y))
      plot(est,type="l",xlab="",ylab="",ylim= ylim.val,
           col=columnCols[6],xlim=xlim.val,lwd=2)
      lines(c(min(est$x),max(est$x)),rep(0,2))

      if (KDEvertLine)
         lines(rep(0,2),c(0,max(est$y)),lwd=3,err=-1,col=KDEvertLineCol)
      
      if (!is.null(addTruthToKDE))
         lines(rep(addTruthToKDE[j],2),c(0,max(est$y)),
               lwd=1,lty=2,,err=-1,col=KDEvertLineCol)
     
      
      x.labels <- pretty(x[[plotInd]][,j],n=3)
      for (i in 1:length(x.labels))
      {
         text(x.labels[i],-0.6*lb.ht,as.character(x.labels[i]),cex=0.8)
         lines(rep(x.labels[i],2),c(0,-0.1*lb.ht))
      }   

      empty.panel()
      text(0.5,0.75,
      paste("posterior mean: ",
             as.character(signif(mean(x[[plotInd]][,j]),3)),sep=""),
             col=columnCols[7],cex=numerSummCex)
      text(0.5,0.5,paste(100*credLevel,"% credible interval: ",
           sep=""),col=columnCols[7],cex=numerSummCex)

      text(0.5,0.25,
      paste("(",
      as.character(signif(quantile(x[[plotInd]][,j],(1-credLevel)/2),3)),",",
      as.character(signif(quantile(x[[plotInd]][,j],(1+credLevel)/2),3)),")",
      sep=""),col=columnCols[7],cex=numerSummCex)
   }

   par(op)
   invisible()
}
