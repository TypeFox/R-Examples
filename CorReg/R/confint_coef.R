# ' plot and give confidence intervals on the coefficients estimated in a model or for proportions
# ' @param modele a model from lm on which to compute the confidence intervals
# ' @param n a vector of quantities associated to prop
# '  @param prop a vector of proportions (between 0 and 1)
# ' @param mean a mean value to plot
# ' @param alpha the risk (confidence 1-alpha)
# ' @param labels a vector of names to put below the bars.
# ' @param subtitle a subtitle to identify the graph
# ' @param lang if not "fr" then in english
# ' @param ylim if needed, a vector c(ymin,ymax).
# ' @export 
# ' 
# ' # '@examples
# '    \dontrun{
# '    require(CorReg)
# ' n=c(3000,1500,45)
# ' prop=c(0.12,0.15,0.11)
# ' confint_coef( n = n, prop = prop, mean = mean(prop),alpha = 0.05)
# '     }
# ' 
#confint_coef<-function(modele=NULL,n=NULL,prop=NULL,mean=NULL,alpha=0.05,labels=NULL,subtitle=NULL,ylim=NULL){  
   confint_coef<-function(n=NULL,prop=NULL,mean=NULL,alpha=0.05,labels=NULL,subtitle=NULL,ylim=NULL){  
      modele=NULL
      if(!is.null(modele)){
      confint=confint.default(modele)
      coef=cbind(modele$coefficients,confint)
      colnames(coef)=c("values","borne inf IC", "borne sup IC")
      plot(rstudent(modele),type="p",cex=0.5,ylab="standardized residuals",main="analyse des residus",sub="95pourcent des points doivent se trouver entre les lignes horizontales")
      abline(h=c(-2,2),col="red")
      barplot(modele$coefficients,ylim=c(0,max(coef[,3])),col="cyan",main="Intervalle de confiance des coefficients",sub="Interpretation : ils ne doivent pas contenir 0")
      x0=seq(from=1,length.out=nrow(coef),by=1.2)-0.3
      y0=coef[,2]
      x1=x0
      y1=coef[,3]
      arrows(x0=x0,y0=y0,x1=x1,y1=y1,angle=90,code=3,lwd=2)
      abline(h=0,col="red")
      print(paste("AIC : ",AIC(modele)))
   }else if(!is.null(prop)){
      coef=list()
      if(length(prop)>1 & !is.null(labels)){
         un=round(prop*n)
         zero=n-un
         names(zero)=labels;names(un)=labels
         X=rbind(cbind(1,rep(labels,times=un)),cbind(0,rep(labels,times=zero)))
         coef=chisq.test(X[,1],X[,2])
      }
      coef$int=c()
      for (i in 1:length(prop)){
         coef$int=rbind(coef$int,c(prop[i],icproportion(prop=prop[i],n=n[i],alpha=alpha)))
      }
      colnames(coef$int)=c("values","borne inf IC", "borne sup IC")
      if(is.null(ylim)){
         ylim=c(0,max(coef$int[,3]))
      }
      barplot(coef$int[,1],ylim=ylim,col="cyan",main="Intervalle de confiance des proportions",sub="",names.arg=labels)
      title(sub=paste(subtitle," p=",coef$p.value))
      x0=seq(from=1,length.out=nrow(coef$int),by=1.2)-0.3
      y0=coef$int[,2]
      x1=x0
      y1=coef$int[,3]
      arrows(x0=x0,y0=y0,x1=x1,y1=y1,angle=90,code=3,lwd=2)
      if(!is.null(mean)){
         abline(h=mean,col="red")
      }
   }
   
   return(coef)
}
