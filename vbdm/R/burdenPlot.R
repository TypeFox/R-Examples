#### Inputs
#### geno is genotype matrix, with SNPs as rows and samples as columns
#### pheno is a vector of phenotypes
### annotation is a vector of annotations, must have same length as number of SNPs, you can just use annotation=rep('missense', times=dim(geno)[1]) if you want
### title is optional title, default is nothing
### order is the order of the SNPs, can be by mean, MAF, or MAF.mean
### legend is whether or not you want a legend, defaults to 'keep', any other value will not produce a legend
### type is which type of graph to produce, 'lines' showing the range of phenotype values for each SNP, or 'points' showing the actual data points
### Missing genotypes are coded as NAs


burdenPlot <- function(y, G, annotation = rep('missense',ncol(G)), title="", order='mean', legend='keep', type='lines',post=NULL,name.snp=NULL){
  geno <- G;
  pheno <- y;
  geno <- t(geno);
  geno[which(is.na(geno))]<-9;
  geno.plot <- subset(geno, apply(geno,1,sum)>0)
  MAC <- NULL
  mean <- NULL
  carriers.value <- list()
  for(i in 1:dim(geno.plot)[1]){
    index.carriers <- which(geno.plot[i,]>0 & geno.plot[i,]!=9)
    carriers.value[[i]] <- pheno[index.carriers]
    mean[i] <- mean(carriers.value[[i]])
    MAC[i] <- length(carriers.value[[i]])
  }

  gt.hom <-apply(geno.plot==0,1,sum)
  gt.het1 <-apply(geno.plot==1,1,sum)
  gt.hom1 <-apply(geno.plot==2,1,sum)
  non.missing <- apply(geno.plot!=9,1,sum)
  foo1 <- cbind(2*gt.hom+gt.het1, gt.het1+2*gt.hom1)
  MAF <- apply(foo1, 1, min)/(2*non.missing)
  
  annotation.plot <- annotation[which(apply(geno,1,sum)>0)]
  
  col <- NULL
  col[annotation.plot=='missense'] <- 'red'
  col[annotation.plot=='nonsense'] <- 'cyan'
  col[annotation.plot=='splice' | annotation.plot=='splice-3' | annotation.plot=='splice-5'] <- 'blue'
  col[annotation.plot!='missense' & annotation.plot!='nonsense' & annotation.plot!='splice-3' & annotation.plot!='splice-3' & annotation.plot!='splice'] <- 'orange'
  
  
  if(order=='mean'){
  carriers.value <- carriers.value[order(mean)]
  col <- col[order(mean)]
  if(!is.null(post)){post <- post[order(mean)]}
  }

  if(order=='MAF'){
  carriers.value <- carriers.value[order(1/MAC)]
  col <- col[order(1/MAC)]
  if(!is.null(post)){post <- post[order(1/MAC)]}
  }    

  if(order=='MAF.mean'){
  carriers.value <- carriers.value[order(1/MAC, mean)]
  col <- col[order(1/MAC, mean)]
  if(!is.null(post)){post <- post[order(1/MAC,mean)]}
  }  

  if(order=='anno'){
  carriers.value <- carriers.value[order(annotation, mean, decreasing=FALSE)]
  col <- col[order(annotation, mean, decreasing=FALSE)]}
	
#if(order=='default'){
#		carriers.value <- carriers.value
#		col <- col;
#	}
  
  n.ticks <- sum(is.na(mean)==FALSE)

  d <- density(pheno)
	mdy <- max(d$y)

  if(is.null(post)){
  plot <- plot(d$x, d$y/3, ylab=" ", xlab='Trait Value',xlim=c(min(d$x),max(d$x)),ylim=c(-mdy, mdy/3),
               type='l', lwd=2, col='blue', yaxt="n",xaxt="n")
  } else{
  plot <- plot(d$x, d$y/3, ylab=" ", xlab='Posterior Probability',xlim=c(min(d$x),max(d$x)),ylim=c(-mdy, mdy/3),
               type='l', lwd=2, col='blue', yaxt="n",xaxt="n")
  }
  if(is.null(post)){
	axis(side=1,at=seq(min(d$x),max(d$x),length.out=6),labels=signif(seq(min(d$x),max(d$x),length.out=6),2));
  }else{
	axis(side=1,at=seq(min(d$x),max(d$x),length.out=6),labels=seq(0,1,length.out=6));
  }
   #print(mdy)
  if(type=='points'){
  points(carriers.value[[1]], -(mdy)/rep(n.ticks, times=length(carriers.value[[1]]))-5e-2, pch=20)}else{
  
  lines(c(min(carriers.value[[1]]), max(carriers.value[[1]])), -(mdy)/rep(n.ticks, times=2)-5e-2, lty='solid', lwd=1)
  points(carriers.value[[1]], -(mdy)/rep(n.ticks, times=length(carriers.value[[1]]))-5e-2, pch=124,cex=.68)
  }
	if(!is.null(name.snp)){
		axis(2,at=seq(-mdy/n.ticks-5e-2,-length(carriers.value)*mdy/n.ticks,length.out=length(carriers.value)),labels=name.snp[order(mean)],cex.axis=.5,las=2);
	}
  points(mean(carriers.value[[1]]), -(mdy)/n.ticks-5e-2, pch=21, col=col[1], cex=1, bg=col[1])
 
  #posterior probability
  if(!is.null(post)){
	  points(post[[1]]*(max(d$x)-min(d$x))+min(d$x),-(mdy)/n.ticks-5e-2,pch=17,col='purple',cex=1,bg='purple');
  }
  
  if(dim(geno.plot)[1]>1){
    for(i in 2:length(carriers.value)){
      if(type=='points'){
      points(carriers.value[[i]], -(i*mdy)/rep(n.ticks, times=length(carriers.value[[i]]))-5e-2+((i-1)*5e-2)/(length(carriers.value)-1), pch=20)}else{
      lines(c(min(carriers.value[[i]]), max(carriers.value[[i]])), -(i*mdy)/rep(n.ticks, times=2)-5e-2+((i-1)*5e-2)/(length(carriers.value)-1), lty='solid', lwd=1)
      points(carriers.value[[i]], -(i*mdy)/rep(n.ticks, times=length(carriers.value[[i]]))-5e-2+((i-1)*5e-2)/(length(carriers.value)-1), pch=124,cex=.68)}
      points(mean(carriers.value[[i]]), -(i*mdy)/n.ticks-5e-2+((i-1)*5e-2)/(length(carriers.value)-1), pch=21, col=col[i], cex=1, bg=col[i])
      if(!is.null(post)){
	  points(post[[i]]*(max(d$x)-min(d$x))+min(d$x),-(i*mdy)/n.ticks-5e-2+((i-1)*5e-2)/(length(carriers.value)-1),pch=17,col='purple',cex=1,bg='purple');
      }
    }
  }
  abline(h=0, col='black')
  if(!is.null(post)){
	axis(side=1,pos=0,at=seq(min(d$x),max(d$x),length.out=6),labels=signif(seq(min(d$x),max(d$x),length.out=6),2));
        text(mean(pheno)+sd(pheno)/2,-(1*mdy)/n.ticks,'Trait Value');
  }
  #segments(x0=mean(pheno), y0=0, x1=mean(pheno), y1=max(d$y/3)+1, col='blue', lty='dashed')
  segments(x0=mean(pheno), y0= -((i+1)*mdy)/n.ticks, x1=mean(pheno), y1=max(d$y/3)+1, col='blue', lty='dashed')

  if(legend=='keep'){
   if(is.null(post)){
    legend('topleft', pch=21, legend=c('missense', 'nonsense', 'splice', 'synonymous'), cex=0.85, 
           col=c('red','cyan','blue','orange'), pt.bg=c('red','cyan','blue','orange'))
   }
   if(!is.null(post)){
    legend('topleft', pch=c(21,17), legend=c('missense',expression(p[j])), cex=0.85, 
           col=c('red','purple'), pt.bg=c('red','purple'))

   }

}
  
  title(title)
  return(plot)
}



###################
###### Example

### Genotype matrix
#MAF <- c(1/7020, 1/7020, 1/7020, 1/7020, 1/7020, 8/7020, 1/7020, 1/7020, 2/7020, 1/7020, 7/7020, 1/7020, 1/7020, 1/7020, 2/7020, 2/7020, 87/7020, 112/7020)
#geno <- matrix(0, nrow=18, ncol=2000)
#for(i in 1:length(MAF)){
#  geno[i,] <- sample(x=c(0,1,2), prob=c((1-MAF[i])^2, 2*MAF[i]*(1-MAF[i]), MAF[i]^2), replace=TRUE, size=2000)
#}
#geno.random <- replace(geno, sample(1:36000, size=1800, replace=FALSE), NA)


#### Phenotypes are 2000 random N(0,1) random variables
#y <- rnorm(2000)

### Annotation is all missense
#annotation=rep('missense', times=dim(geno)[1])

#burden.plot(geno.random, y, annotation)
  








