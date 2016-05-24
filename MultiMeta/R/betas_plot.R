#'Summary plot of the effect sizes and their correlations for a SNP of interest
#'
#'\code{betas_plot} returns a pdf file containing information about the meta-analysis combined effect sizes (beta coefficients) for 
#'a chosen SNP of interest.
#'
#'The plot produced can be a useful tool for visualizing effect sizes across all traits, and the correlation between them. The 
#'plot will have two panels:
#'\enumerate{
#'\item effect sizes with 95\% confidence interval for each trait; 
#'\item a heat-map showing correlations between betas.
#'}
#'Two different types of input are allowed: 
#'\itemize{ 
#'\item a data.frame, containing the output from the \code{multi_meta} function for at least the specified SNP. This is the best choice if results have 
#'already been loaded into the R workspace.  
#'\item a file name containing the output from \code{multi_meta} for at least the specified SNP. The separator field is set to white space. 
#'This file can be a subset of the complete results file (suggested), e.g. a file containing only significant SNPs.  
#'} 
#'
#'@param SNP A SNP of interest, e.g. a SNP that is significantly associated to the analysed traits. 
#'@param res.file File name of the \code{multi_meta} results containing effect sizes for the above \code{SNP}. This is an alternative to \code{data}.
#'@param data Object containing results from the \code{multi_meta} function and the above \code{SNP}. This is an alternative to \code{file}.
#'@param trait.names Vector of analysed trait names to appear on the plot. 
#'@param out.file File name for the output plot.
#'@param col Colours for the heat-map plot. Combinations of several colours are possible. If left blank a heat-map in scales of grey is produced.
#'@param orig.files List of files with multivariate results for each cohort.
#'@param cohort.names Names of the cohorts included in the meta-analysis.

#'@return The output is a pdf file containing the described plot. 
#'@import ggplot2 reshape2 gtable grid
#' 
#'@export



betas_plot=function(SNP, res.file=NULL, data=NULL, trait.names=NULL, out.file="Betas.pdf",col=NULL, orig.files=NULL, 
                    cohort.names=NULL)
{
  x<-NULL
  y<-NULL
  ylo<-NULL
  yhi<-NULL
  Var1<-NULL
  Var2<-NULL
  value<-NULL
  
  ins_unit=function (x, values) 
  {
    lengx <- length(x)
    if (lengx <= 0) {
      unit.c(values, x)
    }else{
      unit.c(x, values)
    }
  }
  
  bind_grobs=function (a, b) 
  {
    b$layout$l <- b$layout$l + ncol(a)
    b$layout$r <- b$layout$r + ncol(a)
    a$layout <- rbind(a$layout, b$layout)
    a$widths <- ins_unit(a$widths, b$widths)
    a$colnames <- c(a$colnames, b$colnames)
    size <- "first"
    
    a_val <- unclass(a$heights)
    b_val <- unclass(b$heights)
    a_unit <- attr(a$heights, "unit")
    a$heights <- switch(size, first = a$heights, last = b$heights, 
                        min = unit(pmin(a_val, b_val), a_unit), 
                        max = unit(pmax(a_val,b_val), a_unit))
    a$grobs <- append(a$grobs, b$grobs)
    a
  }
  
  
  #data input
  if(length(res.file) !=0){
    he=scan(res.file, nlines=1, "char", quiet=T)
    con1=pipe(paste("grep -w ",SNP," ",res.file,sep=""))
    res=scan(con1, quiet=T, "char")
    close(con1)
    names(res)=he
  }else{
    he=names(data)
    res=data[which(data$SNP==SNP),]
  }
  #compute the number of phenotypes
  betas=grep(pattern="beta", he)
  L=length(betas)
  nphen=(-3+sqrt(9+8*L))/2
  
  #first part: forest plot of the combined meta-analysis betas for each trait
  
  a=matrix(nrow=nphen, ncol=4)
  a=as.data.frame(a)
  if(length(trait.names)==0){
    trait.names=paste("Trait",c(1:nphen),sep="_")
    a[,1]=trait.names
  }else{
    a[,1]=trait.names}
  a[,2]=as.numeric(res[paste("beta_",c(1:nphen),sep="")])
  tmp=c()
  for(h in 1:nphen){
    tmp=c(tmp,as.numeric(res[paste("Vbeta_",h,"_",h,sep="")]))
  }
  a[,3]=a[,2]-1.96*sqrt(tmp)
  a[,4]=a[,2]+1.96*sqrt(tmp)
  names(a)=c("x","y","ylo","yhi")
  dimnames(a$y)=c()
  dimnames(a$ylo)=c()
  dimnames(a$yhi)=c()
  one=ggplot(data=a, aes(x=c(1:nphen), y=y, ymin=ylo, ymax=yhi), environment=environment())+
    geom_pointrange(size=1.2)+ 
    geom_hline(aes(x=0), lty=2)+
    coord_flip()+ 
    scale_x_continuous(expand=c(0,0), limits=c(0.5,(nphen+0.5)), labels=trait.names, breaks=c(1:nphen))+
    theme_grey(base_size=10)+
    ylab("Effect size with 95% confidence interval")+
    xlab(NULL)+
    theme(plot.margin=unit(c(0,0,0,0.5), "cm"), axis.ticks=element_line(colour="white"),
          axis.text=element_text(size=12), axis.title.x=element_text(size=15, vjust=0), aspect.ratio=1)
  
  if(length(orig.files)>0){
    num.coh=length(orig.files)
    orig.data=list()
    A=list()
    A.tot=c()
    for(p in 1:num.coh){
      he=scan(orig.files[p], nlines=1, "char", quiet=T)
      
      con2=pipe(paste("grep ",SNP," ",orig.files[p]))
      orig.data[[p]]=scan(con2, quiet=T, "char")
      close(con2)
      names(orig.data[[p]])=he
      
      A[[p]]=matrix(nrow=nphen, ncol=5)
      A[[p]]=as.data.frame(A[[p]])
      A[[p]][,1]=c(1:nphen-p*0.1)
      A[[p]][,2]=as.numeric(orig.data[[p]][paste("beta_",c(1:nphen),sep="")])
      A[[p]][,5]=cohort.names[p]
      tmp=c()
      for(h in 1:nphen){
        tmp=c(tmp,as.numeric(orig.data[[p]][paste("Vbeta_",h,"_",h,sep="")]))
      }
      A[[p]][,3]=A[[p]][,2]-1.96*sqrt(tmp)
      A[[p]][,4]=A[[p]][,2]+1.96*sqrt(tmp)
      names(A[[p]])=c("x","y","ylo","yhi")
      dimnames(A[[p]]$y)=c()
      dimnames(A[[p]]$ylo)=c()
      dimnames(A[[p]]$yhi)=c()
      dimnames(A[[p]])[[2]][5]="col"
      
      A.tot=rbind(A.tot,A[[p]])
    }
    one=one+geom_pointrange(data=A.tot, aes(x=x, y=y, ymin=ylo, ymax=yhi, colour=factor(A.tot$col)))+
      theme(plot.margin=unit(c(0,0,0,0), "cm"),legend.position="left",legend.margin = unit(0, "cm"), legend.text=element_text(angle=45))+
      labs(colour="Cohorts")
    
  }
  
  
  #second part
  #heatmap
  co.v=matrix(nrow=nphen,ncol=nphen) #build matrix of correlation
  for(i in 1:nphen){
    co.v[i:nphen, i]=as.numeric(unlist(res[grep(paste("Vbeta_",i,sep=""), names(res))]))
  }
  for(i in 1:nrow(co.v)) co.v[i,]=co.v[,i]
  co.r=cov2cor(co.v)
  dat=melt(co.r)
  #set the colour pattern
  if(length(col)!=0){
    color=colorRampPalette(col, space="rgb")
  }else{color=colorRampPalette(c("white","grey","black"), space="rgb")}
  
  two<-ggplot(dat, aes(x=Var1, y=Var2, fill=value))+
    geom_tile()+
    scale_fill_gradientn(colours=color(100), limits=c(-1,1), name="Correlation")+
    coord_equal()+
    theme_grey(base_size=10)+
    scale_x_continuous(expand=c(0,0), breaks=c(1:nphen), labels=trait.names, limits=c(0.5,(nphen+0.5)))+
    scale_y_discrete(expand=c(0,0), breaks=c(1:nphen), labels=trait.names, limits=c(1:nphen))+
    theme(plot.margin=unit(c(0,0.5,0,0), "cm"), axis.ticks=element_line(colour="white"),
          axis.text=element_text(size=12), axis.title.x=element_text(size=15, vjust=0), legend.text.align=1)+ 
    ylab("")+
    xlab("Correlation among effect sizes")
  
  #put the plots together 
  g1<-ggplotGrob(one)
  g2<-ggplotGrob(two)
  g<-bind_grobs(g1,g2)
  g<-gtable_add_rows(g, unit(1,"cm"), 0)
  g3<-textGrob(paste("Summary of the combined effects for: ",SNP,sep="" ), gp=gpar(fontsize=20))
  g<-gtable_add_grob(g, g3, t=1, b=1, l=7, r=7, clip='off')
  
  #build pdf file                      
  pdf(out.file,14,7)
  grid.draw(g)
  garbage<-dev.off()
  
  
  
}