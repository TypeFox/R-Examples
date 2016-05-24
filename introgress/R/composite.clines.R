composite.clines <-
function(cline.data=NULL,pdf=TRUE,out.file="comp.pdf",colors=c("#005A32","#41AB5D"),labels=c("AA","Aa")){ 
  if (is.null(cline.data)==TRUE)
    stop("error, input file was not provided")
  if (is.list(cline.data)==FALSE)
    stop("error, cline.data input file not in list format")
  if (length(colors)!=2) stop ("colors must be a vector of length two.")
  ## set up outfile
  if(pdf==TRUE) pdf(file=out.file,width=8,height=4)
  par(mfrow=c(1,2))
  par(mar=c(4.5,4.5,2.5,3))

  ## determine whether significance testing was conducted
  if(is.null(cline.data[[5]])==FALSE) sig<-1
  else sig<-0

  ## homozygote
  plot(0:1,0:1,type="n",xlab="Hybrid index",ylab=labels[1])
  n.loci<-dim(cline.data$Loci.data)[1]

  ## Get needed values
  hi<-cline.data$hybrid.index
  if (sig==1){
    ub<-cline.data$Neutral.AA[[1]][1,]
    lb<-cline.data$Neutral.AA[[2]][1,]
  }
  
  ## Plot neutral cline
  if (sig==1){
    polygon.data.matrix<-rbind(hi,lb,ub,lb,ub)[,order(hi)]
    polygon.data.matrix.rev<-rbind(hi,lb,ub,lb,ub)[,order(-hi)]
    polygon(c(polygon.data.matrix[1,],polygon.data.matrix.rev[1,]),
            c(polygon.data.matrix[2,],polygon.data.matrix.rev[3,]), col=colors[1],border=NA)
  }

  ## Plot observed Cline
  for(i in 1:n.loci){
    AA.line<-cline.data$Fitted.AA[i,]
    line.matrix<-rbind(hi,AA.line)[,order(hi)]
    lines(line.matrix[1,],line.matrix[2,],lty=1) 
  }
  
  ## heterozygote
  plot(0:1,0:1,type="n",xlab="Hybrid index",ylab=labels[2])
  n.loci<-dim(cline.data$Loci.data)[1]

  ## Get needed values
  hi<-cline.data$hybrid.index
  if (sig==1){
    ub<-cline.data$Neutral.Aa[[1]][1,]
    lb<-cline.data$Neutral.Aa[[2]][1,]
  }
  
  ## Plot neutral cline
  if (sig==1){
    polygon.data.matrix<-rbind(hi,lb,ub,lb,ub)[,order(hi)]
    polygon.data.matrix.rev<-rbind(hi,lb,ub,lb,ub)[,order(-hi)]
    polygon(c(polygon.data.matrix[1,],polygon.data.matrix.rev[1,]),
            c(polygon.data.matrix[2,],polygon.data.matrix.rev[3,]),col=colors[2],border=NA)
  }

  ## Plot observed Cline
  for(i in 1:n.loci){
    Aa.line<-cline.data$Fitted.Aa[i,]
    line.matrix<-rbind(hi,Aa.line)[,order(hi)]
    lines(line.matrix[1,],line.matrix[2,],lty=1) 
  }	
  if(pdf==TRUE){
    dev.off()
  }
}

