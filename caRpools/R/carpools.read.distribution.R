carpools.read.distribution=function(dataset,namecolumn=1, fullmatchcolumn=2,breaks="", title="Title", xlab="X-Axis", ylab="Y-Axis",statistics=TRUE, col=rgb(0, 0, 0, alpha = 0.65),extractpattern=expression("^(.+?)_.+"), plotgene=NULL, type="distribution", logscale=TRUE){   
  
  # dataset = read.table frame
  # fullmatchcolumn = integer, in which column is the matchcount
  # xlim = integer, which row should be the last row to display (e.g. if random genes make up most, you can exzlude them)
  # TYPE - distribution | whisker
  
  
  
  if(!is.null(plotgene))
  {
    gene.names = sub(extractpattern,"\\1",dataset[,namecolumn],perl=TRUE)
    dataset$genes = as.character(gene.names)
    # number of genes
    numberdata = nrow(dataset[dataset$genes %in% plotgene,])
    if(is.null(numberdata) | numberdata <1)
    {
      plotgene = gene.names
      numberdata = nrow(dataset[dataset$genes %in% plotgene,])
    }
    if (breaks=="" || is.null(breaks))
    { breaks = numberdata/2 }
    
    dataset2 = dataset[dataset$genes %in% plotgene,]

    
    
    # compute statistics
    meancount=as.integer(mean(dataset2[,fullmatchcolumn]))
    mediancount= as.integer(median(dataset2[,fullmatchcolumn]))
    
    # make log2
    if(logscale==TRUE)
    {
      dataset2[,fullmatchcolumn] = log2(dataset2[,fullmatchcolumn])
      
      # remove NAs ? infinites and set them to 0
      dataset2[,fullmatchcolumn] = apply(dataset2, 1, function(x) if(is.finite(as.numeric(x[fullmatchcolumn]))) {return(as.numeric(x[fullmatchcolumn]))}
                                                                              else {return(as.numeric(0))})
    }
    
    maxcount=max(dataset2[,fullmatchcolumn])
    
    xlim=max(dataset2[,fullmatchcolumn])
    
    # make plot TEXTS changeable
    # make plot working for genes AND designs
    if(type=="distribution")
    {
    #plot without LIMIT
    histinfo = hist(dataset2[,fullmatchcolumn],breaks=breaks,xlim=c(0,xlim), col=col, xlab=xlab, ylab=ylab, main=title, border=FALSE)
    
    if(statistics)
    {
      legend("topright",c(paste("MEAN:",meancount),paste("MEDIAN:",mediancount)),cex=0.8, bty="n", text.col=c("red","orange"))
      
      if(logscale==TRUE)
      {
        meancount = log2(meancount)
        mediancount = log2(mediancount)  
      }
      
      lines(x=rep(meancount, times=max(histinfo$counts)+1), y=c(0:max(histinfo$counts)), col="red", lwd=2, lty=2)
      lines(x=rep(mediancount, times=max(histinfo$counts)+1), y=c(0:max(histinfo$counts)), col="orange", lwd=2, lty=2)
      
      # add density function as a fit
      par(new=TRUE)
      hist(dataset2[,fullmatchcolumn],breaks=breaks,xlim=c(0,xlim), col=rgb(1, 1, 1, 0),prob=TRUE, axes=F, xlab="", ylab="", main="", border=FALSE)
      lines(density(dataset2[,fullmatchcolumn]), col="navy", lwd=3)
      
      }
    } # end type

  }
  else # some special sgRNAs will be plotted
  {
    
    # number of genes
    numberdata=nrow(dataset)
    if (breaks=="")
    {breaks = numberdata/2}
    
    # compute statistics
    meancount=as.integer(mean(dataset[,fullmatchcolumn]))
    mediancount= as.integer(median(dataset[,fullmatchcolumn]))
    maxcount=max(dataset[,fullmatchcolumn])
    
    if(logscale==TRUE)
    {
      dataset[,fullmatchcolumn] = log2(dataset[,fullmatchcolumn])
      
      # remove NAs ? infinites and set them to 0
      dataset[,fullmatchcolumn] = apply(dataset, 1, function(x) if(is.finite(as.numeric(x[fullmatchcolumn]))) {return(as.numeric(x[fullmatchcolumn]))}
                                         else {return(as.numeric(0))})
    }
    
    # make plot TEXTS changeable
    # make plot working for genes AND designs
    if(type=="distribution") {
    #plot without LIMIT
    xlim=max(dataset[,fullmatchcolumn])
      
    histinfo = hist(dataset[,fullmatchcolumn],breaks=breaks,xlim=c(0,xlim), col=col, xlab=xlab, ylab=ylab, main=title, border=FALSE)
    #print(histinfo)
    if(statistics)
    {
      legend("topright",c(paste("MEAN:",meancount),paste("MEDIAN:",mediancount)),cex=0.8, bty="n", text.col=c("red","orange"))
      
      if(logscale==TRUE)
      {
        meancount = log2(meancount)
        mediancount = log2(mediancount)  
      }
      
      lines(x=rep(meancount, times=max(histinfo$counts)+1), y=c(0:max(histinfo$counts)), col="red", lwd=2, lty=2)
      lines(x=rep(mediancount, times=max(histinfo$counts)+1), y=c(0:max(histinfo$counts)), col="orange", lwd=2, lty=2)
      
      
      # add density function as a fit
      par(new=TRUE)
      hist(dataset[,fullmatchcolumn],breaks=breaks,xlim=c(0,xlim), col=rgb(1, 1, 1, 0), prob=TRUE, axes=F, xlab="", ylab="", main="", border=FALSE)
      lines(density(dataset[,fullmatchcolumn]), col="navy", lwd=3)
      
    }
    #END TYPE
    } else if(type=="whisker")
      
    {
      boxplot(dataset[,fullmatchcolumn],
                  col=col,
                  pch=16,
                  main=title,
                  horizontal=TRUE,
                  xlab=xlab,
              #range=IQR(dataset[,fullmatchcolumn]),
                  las=2)
    
      legend("topright",c(paste("MEAN:",meancount),paste("MEDIAN:",mediancount)),cex=0.8, bty="n", text.col=c("red","red","orange"))
      
    }# end type whisker
    
  }
  
}

