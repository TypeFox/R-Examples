carpools.read.depth=function(
  datasets,
  namecolumn=1 ,
  fullmatchcolumn=2,
  dataset.names=NULL,
  extractpattern=expression("^(.+?)_.+"),
  col=rgb(0, 0, 0, alpha = 0.65),
  xlab="Genes",
  ylab="Read Count per sgRNA",
  statistics=TRUE,
  labelgenes = NULL,
  controls.target = controls.target ,
  controls.nontarget=controls.nontarget,
  labelcolor="orange", waterfall=FALSE){
  
  shadowtext <- function(x, y=NULL, labels, col='white', bg='black',
                         theta= seq(pi/4, 2*pi, length.out=8), r=0.1, ... ) {
    
    xy <- xy.coords(x,y)
    xo <- r*strwidth('A')
    yo <- r*strheight('A')
    
    for (i in theta) {
      text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
    }
    text(xy$x, xy$y, labels, col="black", ... )
  }

  # Prepare for multiple datasets
  #dataset = as.data.frame(0)
  
  # Dataset names provided?
   
  if(is.null(dataset.names) || length(datasets) != length(dataset.names) ){
    stop("No dataset names provided (dataset.names =c()) or number of names and datasets do not match.")
  }    
  dataset.number = length(datasets)  
  ## Iterate
  # Backup plotting pars and set for multiple plotting
  par.old = par
  par(oma=c(0,0,0,0))  
  if(dataset.number >=2)
  {mfrow=c(2,1)} else {mfrow=c(1,1)}
  # Add Names
  names(datasets) = dataset.names  
  data.count=0  
  # go for each dataset  
  for(data.name in names(datasets))
  {
    n=datasets[[data.name]]     
    
    if(is.null(n$Gene))
    {n$Gene  = as.character(sub(extractpattern,"\\1",n[,namecolumn],perl=TRUE))}
    
    names(n)=c("sgRNAs","Readcount","gene_name")    
    names = n$sgRNAs
    # Extract gene identifiers by extractpattern and add to dataframe
    # get and Aggregate gene names
    # prepare data for plotting
    # plot genes on x -axis (aggregated) to readcount normalised for number of sgRNA on y axis
    dataset.plot = aggregate(n$Readcount,list(sub(extractpattern,"\\1",names,perl=TRUE)),function(x) {
        return(c(
          sum(x,na.rm=T)/length(x[as.numeric(x)>0])
        )
        )
      }
    )
    dataset.plot$labelgene=col
    names(dataset.plot)=c("gene.label","persgRNACount","labelgene")
    # Label genes/designs?
    if(!is.null(labelgenes))
    {
      dataset.plot$labelgene[which(dataset.plot$gene.label %in% labelgenes)]=labelcolor
    }
 
    if(!is.null(controls.target))   {
        dataset.plot$labelgene[which(dataset.plot$gene.label %in% controls.target)]="red3"
    }
      
    if(!is.null(controls.nontarget)){
        dataset.plot$labelgene[which(dataset.plot$gene.label %in% controls.nontarget)]="blue"
    }
    
        
    # In some datasets, only ONE sgRNA might be present with read Count of 0!
    # In this case, we need to remove the data from the dataset!
    todelete = which(dataset.plot$persgRNACount %in% 0)
    if(length(todelete) > 0){
      dataset.plot = dataset.plot[-todelete,]
    }
     
    dataset.plot$width = 1
    
    if( !is.null(labelgenes) || 
        !is.null(controls.target) || 
        !is.null(controls.nontarget)
        ){
        # if gene is to be labelled, make it wider, so that it can be seen
        dataset.plot$width = apply(dataset.plot, 1, function(i) { 
                                                        if(as.character(i["gene.label"]) %in% c(labelgenes,controls.target,controls.nontarget)) {
                                                            4
                                                          }else{
                                                            1
                                                          }
                                                        }
        )
    }
    
    # setup dataset information
    numberdata=nrow(dataset.plot)
    numberstats=as.integer(numberdata/4)
    
    #calculate reads per gene for MATCH, SEEDMATCH, MISMATCH -> not normalized to number of designs
    #compute statistics
    sumdata=summary(dataset.plot)
    meandata=mean(dataset.plot$persgRNACount)
    mediandata=median(dataset.plot$persgRNACount)
    mindata=min(dataset.plot$persgRNACount)
    maxdata=max(dataset.plot$persgRNACount)
    sddata=sd(dataset.plot$persgRNACount)
    
    #remove NaN or NA and set to 0!
    dataset.plot$persgRNACount = sapply(dataset.plot$persgRNACount, function(x)
      {
      if(is.na(x) || is.infinite(x))
      {return(as.numeric(0))}
      else {return(as.numeric(x))}
    })
    
    
    if(data.count==0){
      # Set minimum and maximum to have same scaling accoridng to DMSO Controls
      plot.min=0
      plot.max=max(dataset.plot$persgRNACount)
      data.count=1
    }
    
    # sort by count
    if(identical(waterfall,TRUE))
    {
      dataset.plot = dataset.plot[order(dataset.plot$persgRNACount, decreasing=TRUE),]
    }
    
    barplot(dataset.plot$persgRNACount, xlim=c(-(numberdata/17),numberdata*1.2), ylim=c(plot.min, plot.max), main=paste("Gene Read Count per sgRNA",data.name, sep="\n"),
            xlab=xlab, ylab=ylab,
            col=dataset.plot$labelgene,
            border=FALSE, legend = FALSE, width=dataset.plot$width)
    
    if(statistics)
    {  
      lines(y=rep(meandata,times=numberstats),x=c(1:numberstats), col="red",lty="solid", lwd=2)
      shadowtext(y=meandata, x=0, col="red", labels="Mean", bg="white", adj=c(1.2,0), cex=0.7)
      
      lines(y=rep(meandata+sddata,times=numberstats),x=c(1:numberstats), col="red",lty="dashed", lwd=2)
      shadowtext(y=meandata+sddata, x=0, col="red", labels="SD", bg="white",cex=0.7, adj=c(1.2,0))
      lines(y=rep(meandata-sddata,times=numberstats),x=c(1:numberstats), col="red",lty="dashed", lwd=2)
      
      lines(y=rep(mindata,times=numberstats),x=c(1:numberstats), col="blue",lty="solid", lwd=2)
      shadowtext(y=mindata, x=0, col="blue", labels="Min", bg="white", adj=c(1.2,-1),cex=0.7)
      lines(y=rep(maxdata,times=numberstats),x=c(1:numberstats), col="blue",lty="solid", lwd=2)
      shadowtext(y=maxdata, x=0, col="blue", labels="Max", bg="white", adj=c(1.2,2),cex=0.7)
    }
    
    # Plot legend
    legend("topright",c(paste("MIN:",round(mindata,digits = 2)), paste("MEAN:",floor(meandata)),paste("MAX:",round(maxdata,digits = 2))),cex=0.8, bty="n", text.col=c("red","orange"))
  }
}