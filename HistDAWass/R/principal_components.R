## 1v PCA for quantiles -------
# Principal components analysis  of histogram variable based on Wasserstein distance----
#' Principal components analysis  of histogram variable based on Wasserstein distance
#' @description The function implements a Principal components analysis of histogram variable 
#' based on Wasserstein distance. It performs a centered (not standardized) PCA on a set of quantiles of a variable.
#' Being a distribution a multivalued description, the analysis performs a dimensional reduction and a visualization of distributions. 
#' It is a 1d (one dimension) becuse it is considered just one histogram variable.
#' @param data A MatH object (a matrix of distributionH).
#' @param var An integer, the variable number.
#' @param quantiles An integer, it is the number of quantiles used in the analysis.
#' @param plots a logical value. Default=TRUE plots are drawn.
#' @param listaxes A vector of integers listing the axis for the 2d factorial reperesntations.
#' @param axisequal A logical value. Default TRUE, the plot have the same scale for the x and the y axes.
#' @param qcut  a number between 0.5 and 1, it is used for the representation of densities for avoiding very peaked densities in the plot. 
#' Default=1, all the densities are considered.

#' @return a list with the results of the PCA in the MFA format of package \pkg{FactoMineR} for function MFA
#' @references Verde, R.; Irpino, A.; Balzanella, A., "Dimension Reduction Techniques for Distributional Symbolic Data," Cybernetics, IEEE Transactions on , vol.PP, no.99, pp.1,1
#' doi: 10.1109/TCYB.2015.2389653
#' keywords: {Correlation;Covariance matrices;Distribution functions;Histograms;Measurement;Principal component analysis;Shape;Distributional data;Wasserstein distance;principal components analysis;quantiles},
#' \url{http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7024099&isnumber=6352949}
#' @details In the framework of symbolic data analysis (SDA), distribution-valued data 
#' are defined as multivalued data, where each unit is described by a distribution 
#' (e.g., a histogram, a density, or a quantile function) of a quantitative variable.
#'  SDA provides different methods for analyzing multivalued data. Among them, the most relevant 
#'  techniques proposed for a dimensional reduction of multivalued quantitative variables is principal 
#'  component analysis (PCA). This paper gives a contribution in this context of analysis. 
#'  Starting from new association measures for distributional variables based on a peculiar metric 
#'  for distributions, the squared Wasserstein distance, a PCA approach is proposed for distribution-valued data,
#'   represented by quantile-variables. 
#' @examples
#' results=WH.1d.PCA(data = BLOOD,var = 1, listaxes=c(1:2))
#' @importFrom FactoMineR PCA
#' @importFrom graphics plot abline axis hist plot.new plot.window polygon segments text title
#' @importFrom grDevices dev.new rgb
#' @importFrom stats density quantile
#' @export
WH.1d.PCA=function(data,var, quantiles=10, plots=TRUE, listaxes=c(1:4),axisequal=FALSE,qcut=1){
  if (is(data)[1]!="MatH"){stop("Input data must be a MatH object (A matrix of histograms)")}
  #Build the matrix of Quantiles
  VARS=ncol(data@M)
  INDIV=nrow(data@M)
  if ((qcut<0.5)||(qcut>1)){qcut=1}
  if (missing(var)){
    var=1
    varname=colnames(data@M)[1]
    cat(paste("Var is missing, We do a PCA on variable ", varname, "\n"))
  } else{
    if (length(var)>1){
      varname=colnames(data@M)[var[1]]
      cat(paste("Var has several values, We do a PCA on variable -->", var[1], "\n"))
      var=var[1]
    }
    else{
      if (var>VARS){stop(paste("The variables are less than ",var))}
      else{
        varname=colnames(data@M)[var]
        cat(paste("We do a PCA on variable ---> ", varname, "\n"))
      }
    }
  }
  data=data[,var]
  #check indiv and remove empty rows
  toremove=numeric(0)
  for (i in 1:INDIV){
    if (length(data@M[i,1][[1]]@x)<2){
      toremove=c(toremove,i)
    }
  }
  if (length(toremove)>0){
    data@M=as.matrix(data@M[-toremove,1])
  }
  colnames(data@M)=varname
  INDIV=nrow(data@M)
  if (INDIV<2){stop("At least two individuals are necessary, please check data")}
  
  #Create the quantile matrix
  MATQ=matrix(0,nrow = INDIV,ncol = (quantiles+1))
  p=c(0, 1:quantiles)/quantiles
  namec=list("Min")
  for (j in 2:(quantiles)){
    namec=c(namec,paste0("Q(",format(p[j], digits=2, nsmall=2),")"))
  }
  namec=c(namec, "Max")
  
  for (i in 1:INDIV)
  {
    for (j in 1:(quantiles+1))
    {MATQ[i,j]=compQ(data@M[i,1][[1]],p[j])}
    
  }
  rownames(MATQ)=rownames(data@M)
  colnames(MATQ)=namec
  
  # do a PCA
  vmeans=numeric(0)
  vstds=numeric(0)
  
  for (ind in 1:INDIV){
    vmeans=c(vmeans,data@M[ind,1][[1]]@m)
    vstds=c(vstds,data@M[ind,1][[1]]@s)
  }
  minM=min(vmeans)
  maxM=max(vmeans)
  minS=min(vstds)
  maxS=max(vstds)
  
  MATF=cbind(MATQ/sqrt((quantiles+1)),vmeans,vstds)
  res.pca = PCA(MATF,quanti.sup = c((ncol(MATF)-1),ncol(MATF)), scale.unit=FALSE,  graph=F)
  #Sparse.pca.res = SPC(MATF[,1:(quantiles+1)],sumabsv=3, K=4, orth=TRUE)
  #Sparse2.pca.res = spca(MATF[,1:(quantiles+1)], K=3,para=c(0.9,0.5,0.3))
  VARW=WH.var.covar(data)
  TOTINE=sum(res.pca$eig[,1])
  
  ## plotting PCA results ----------------
  if (plots){
    #par(mfrow=c(1,1))
    # lt's define the couple of axes
    dev.new(noRStudioGD = FALSE )
    planes=matrix(0,1,2)
    for (cc in 1:floor(length(listaxes)/2)){
      planes=rbind(planes,c(listaxes[(cc*2-1)],listaxes[(cc*2)]))
    }
    planes=as.matrix(planes[-1,])
    dim(planes)=c(floor(length(listaxes)/2),2)
    if ((length(listaxes)%%2)==1){
      planes=rbind(planes, c(listaxes[(length(listaxes)-1)],listaxes[length(listaxes)]))
    }
    
    #plot Spanish-fan plot
    for (pl in 1:nrow(planes)){
      axe1=planes[pl,1];axe2=planes[pl,2]
      
      labX=paste("Axis ",axe1," (", format(res.pca$eig[axe1,2], digits=2, nsmall=2), "%)")
      labY=paste("Axis ",axe2," (", format(res.pca$eig[axe2,2], digits=2, nsmall=2), "%)")
      CVAR=res.pca$var$coord[,c(axe1,axe2)]
      plot.new()
      
      if (axisequal){
        xrange=c(min(0,min(cbind(CVAR[,1],CVAR[,2]))), max(0,max(cbind(CVAR[,1],CVAR[,2]))))
        yrange=xrange
        plot.window(xrange, yrange)
        title(xlab=labX)
        title(ylab=labY)
        
        #         plot(type="n",
        #              xlab=labX,ylab=labY)
        segments(xrange[1],0,xrange[2],0)
        segments(yrange[1],0,yrange[2],0)
      }
      else{
        plot.window(c(min(0,CVAR[,1]),max(0,CVAR[,1])),c(min(0,CVAR[,2]),max(0,CVAR[,2])))
        title(xlab=labX)
        title(ylab=labY)
        #       plot(c(min(0,CVAR[,1]),max(0,CVAR[,1])), c(min(0,CVAR[,2]),max(0,CVAR[,2])),type="n",
        #            xlab=labX,ylab=labY)Y)
        segments(min(0,CVAR[,1]),0,max(0,CVAR[,1]),0)
        segments(0,min(0,CVAR[,2]),0,max(0,CVAR[,2]))}
      title("PCA Variable plot (Spanish-fan plot)")
      if ((nrow(CVAR)%%2)==0){centr=nrow(CVAR)/2} else{centr=(nrow(CVAR)-1)/2}
      centr=nrow(CVAR)/2
      cxl=0.8
      for (tr in 2:nrow(CVAR)){
        centrality=(abs(tr-centr-1)/(centr-1))
        
        #cat(centrality,"\n")
        x=c(0, CVAR[(tr-1),1], CVAR[tr,1], 0)
        y=c(0, CVAR[(tr-1),2], CVAR[tr,2], 0)
        red=1
        green=centrality
        #cat(red,green,"\n")
        polygon(x, y, col=rgb(red, green, 0,0.7), lty = 1, lwd = 1,  border = "black")
      }
      text(x = CVAR[1,1],y= CVAR[1,2],labels = "Min",cex=cxl)
      if (nrow(CVAR)<8){
        for (tr in 2:(nrow(CVAR)-1)){
          text(x = CVAR[tr,1],y= CVAR[tr,2],labels = colnames(MATQ)[tr],cex=cxl)
        }}
      else{
        
        text(x = CVAR[ceiling(nrow(CVAR)/8),1],y= CVAR[ceiling(nrow(CVAR)/8),2],
             labels = colnames(MATQ)[ceiling(nrow(CVAR)/8)],cex=cxl)
        text(x = CVAR[ceiling(nrow(CVAR)*2/8),1],y= CVAR[ceiling(nrow(CVAR)*2/8),2],
             labels = colnames(MATQ)[ceiling(nrow(CVAR)*2/8)],cex=cxl)
        text(x = CVAR[ceiling(nrow(CVAR)*3/8),1],y= CVAR[ceiling(nrow(CVAR)*3/8),2],
             labels = colnames(MATQ)[ceiling(nrow(CVAR)*3/8)],cex=cxl)
        text(x = CVAR[centr+1,1],y= CVAR[centr+1,2],labels = colnames(MATQ)[centr+1],cex=cxl)
        text(x = CVAR[ceiling(nrow(CVAR)*5/8),1],y= CVAR[ceiling(nrow(CVAR)*5/8),2],
             labels = colnames(MATQ)[ceiling(nrow(CVAR)*5/8)],cex=cxl)
        text(x = CVAR[ceiling(nrow(CVAR)*6/8),1],y= CVAR[ceiling(nrow(CVAR)*6/8),2],
             labels = colnames(MATQ)[ceiling(nrow(CVAR)*6/8)],cex=cxl)
        text(x = CVAR[ceiling(nrow(CVAR)*7/8),1],y= CVAR[ceiling(nrow(CVAR)*7/8),2],
             labels = colnames(MATQ)[ceiling(nrow(CVAR)*7/8)],cex=cxl)
        
      }
      text(x = CVAR[nrow(CVAR),1],y= CVAR[nrow(CVAR),2],labels = "Max",cex=cxl)
      axis(1, pos = 0, cex.axis=0.7)#, at=xM,labels=labX)
      axis(2, pos = 0, cex.axis=0.7)#, at=yM,labels=labY)
      dev.new(noRStudioGD = FALSE )
      
      #plot individuals
      
      #layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
      
      xm=min(res.pca$ind$coord[,axe1])
      xM=max(res.pca$ind$coord[,axe1])
      ym=min(res.pca$ind$coord[,axe2])
      yM=max(res.pca$ind$coord[,axe2])
      Xlen=xM-xm
      Ylen=yM-ym
      
      #Compute scaling factor for domain
      matX=matrix(0,128,INDIV)
      matD=matrix(0,128,INDIV)
      for (ind in 1:INDIV){
        distr=data@M[ind,1][[1]]
        #generate 200 random points according to the QF
        rn=200
        xn=c(rep(0,rn))
        random_no=c(0:rn)/rn
        
        for (i in 1:rn){
          xn[i]=compQ(distr,random_no[i])
        }
        d <- density(xn,n=128)
        matX[,ind]=d$x
        matD[,ind]=d$y
      }
      MinX=min(matX)
      MaxX=max(matX)
      RX=MaxX-MinX
      q95D=quantile(matD,probs = qcut)
      matD[matD>q95D]=as.numeric(q95D)
      MaxD=max(matD)
      
      meanXRange=diff(range(apply(matX,2,FUN = range)))
      xfact=0.3
      yfact=0.1
      if (axisequal){xm1=min(xm,ym)
                     xM1=max(xM,yM)
                     ym1=xm1
                     yM1=xM1
                     Xlen=xM1-xm1
                     Ylen=yM1-ym1
                     xm=(xM1+xm1)/2-Xlen/2*(1+1.5*xfact)
                     xM=(xM1+xm1)/2+Xlen/2*(1+1.5*xfact)
                     ym=(yM1+ym1)/2-Ylen/2*(1+0.5*yfact) 
                     yM=(yM1+ym1)/2+Ylen/2*(1+1.5*yfact)
                     
      }else{
        xm1=xm
        xM1=xM
        ym1=ym
        yM1=yM
        xm=(xM1+xm1)/2-Xlen/2*(1+1.5*xfact)
        xM=(xM1+xm1)/2+Xlen/2*(1+1.5*xfact)
        ym=(yM1+ym1)/2-Ylen/2*(1+0.5*yfact) 
        yM=(yM1+ym1)/2+Ylen/2*(1+1.5*yfact)     
      }
      dev.new(noRStudioGD = FALSE )
      plot.new()
      plot.window(c(xm,xM),c(ym,yM))
      title(main="PCA plot of distributions", sub="Smoothed distributions")
      for (ind in 1:INDIV){
        x=matX[,ind]
        y=matD[,ind]
        x=(x-data@M[ind,1][[1]]@m)/meanXRange*xfact*Xlen+res.pca$ind$coord[ind,axe1]
        x=c(x,x[length(x)],x[1])
        y=c(y,0,0)/q95D*yfact*Ylen+res.pca$ind$coord[ind,axe2]
        polygon(x,y, col=rgb(1, 1,0,0.6))
        
      }
      for (ind in 1:INDIV){
        text(res.pca$ind$coord[ind,axe1],res.pca$ind$coord[ind,axe2],
             label=rownames(MATQ)[ind], pos=1, cex=0.7)
      }
      axis(1, pos = 0, cex.axis=0.7)#, at=xM,labels=labX)
      axis(2, pos = 0, cex.axis=0.7)#, at=yM,labels=labY)
      segments(xm,0,xM,0)
      segments(0,ym,0,yM)
      
      title(xlab=labX)
      title(ylab=labY)
      dev.new(noRStudioGD = FALSE )
      # plot.new()
      plot.new()
      plot.window(c(xm,xM),c(ym,yM))
      title(main="Coloured using the mean values")
      trasp=0.6
      
      for (ind in 1:INDIV){
        x=matX[,ind]
        y=matD[,ind]
        x=(x-data@M[ind,1][[1]]@m)/meanXRange*xfact*Xlen+res.pca$ind$coord[ind,axe1]
        x=c(x,x[length(x)],x[1])
        y=c(y,0,0)/q95D*yfact*Ylen+res.pca$ind$coord[ind,axe2]
        red=(data@M[ind,1][[1]]@m-minM)/(maxM-minM)
        if (red>0.5){
          green=1-red
          blue=0
        }
        else{green=1-red
             red=0
             blue=1-green
        }
        polygon(x,y, col=rgb(red,green,blue,trasp))
        
      }
      for (ind in 1:INDIV){
        text(res.pca$ind$coord[ind,axe1],res.pca$ind$coord[ind,axe2],
             label=rownames(MATQ)[ind], pos=1, cex=0.7)
      }
      axis(1, pos = 0, cex.axis=0.7)#, at=xM,labels=labX)
      axis(2, pos = 0, cex.axis=0.7)#, at=yM,labels=labY)
      segments(xm,0,xM,0)
      segments(0,ym,0,yM)
      
      abline(0,res.pca$quanti.sup$coord[1,2]/res.pca$quanti.sup$coord[1,1],lty=3,lwd=2, col="red")
      title(xlab=labX)
      title(ylab=labY)
      dev.new(noRStudioGD = FALSE )
      plot.new()
      plot.window(c(xm,xM),c(ym,yM))
      title(main="Coloured using std values")
      for (ind in 1:INDIV){
        
        x=matX[,ind]
        y=matD[,ind]
        x=(x-data@M[ind,1][[1]]@s)/meanXRange*xfact*Xlen+res.pca$ind$coord[ind,axe1]
        x=c(x,x[length(x)],x[1])
        y=c(y,0,0)/q95D*yfact*Ylen+res.pca$ind$coord[ind,axe2]
        red=(data@M[ind,1][[1]]@s-minS)/(maxS-minS)
        if (red>0.5){
          green=1-red
          blue=0
        }
        else{green=1-red
             red=0
             blue=1-green
        }
        polygon(x,y, col=rgb(red,green,blue,trasp))
        
      }
      for (ind in 1:INDIV){
        text(res.pca$ind$coord[ind,axe1],res.pca$ind$coord[ind,axe2],
             label=rownames(MATQ)[ind], pos=1, cex=0.7)
      }
      axis(1, pos = 0, cex.axis=0.7)#, at=xM,labels=labX)
      axis(2, pos = 0, cex.axis=0.7)#, at=yM,labels=labY)
      segments(xm,0,xM,0)
      segments(0,ym,0,yM)
      
      abline(0,res.pca$quanti.sup$coord[2,2]/res.pca$quanti.sup$coord[2,1],lty=3,lwd=2, col="red")
      
      title(xlab=labX)
      title(ylab=labY)
    }
  }
  return(list(PCAout=res.pca, WASSVARIANCE=VARW, INERTIA=TOTINE))
}
# Principal components analysis  of a set of histogram variables based on Wasserstein distance----
#' (Preliminary version) Principal components analysis  of a set of histogram variable based on Wasserstein distance
#' @description (Beta version) The function implements a Principal components analysis of a set of histogram variables 
#' based on Wasserstein distance. It performs a centered (not standardized) PCA on a set of quantiles of a variable.
#' Being a distribution a multivalued description, the analysis performs a dimensional reduction and a visualization of distributions. 
#' It is a 1d (one dimension) becuse it is considered just one histogram variable.
#' @param data A MatH object (a matrix of distributionH).
#' @param list.of.vars A list of  integers, the active variables.
#' @param quantiles An integer, it is the number of quantiles used in the analysis.
#' Default=1, all the densities are considered.

#' @return a list with the results of the PCA in the MFA format of package \pkg{FactoMineR} for function MFA
#' @details It is a first tentative of extending WH.1d.PCA to the multiple case.All plotting functions must be defined
#' in order to have interpretable results. This is an ongoing function, even if the numerical results are methodologically verified.
#' @importFrom FactoMineR MFA
#' @export
## MFA PCA for quantiles -------
WH.MultiplePCA=function(data, list.of.vars, quantiles=10){
  data=data[,list.of.vars]
  VARS=ncol(data@M)
  INDIV=nrow(data@M)
  MATQ=matrix(0,nrow = INDIV,ncol = VARS*(quantiles+1))
  p=c(0, 1:quantiles)/quantiles
  #create the names of the variables
  COLnames=list();
  for (v in 1:VARS){
    namec=list("Min")
    for (j in 2:(quantiles)){
      namec=c(namec,paste0("Q(",format(p[j], digits=2, nsmall=2),")"))
    }
    namec=c(namec, "Max")
    COLnames=c(COLnames, paste(substr(colnames(data@M)[v],1,4),namec,sep="."))
    for (i in 1:INDIV)
    {
      for (j in 1:(quantiles+1)){
        MATQ[i,((v-1)*(quantiles+1)+j)]=compQ(data@M[i,v][[1]],p[j])
      }
    }
  }
  rownames(MATQ)=rownames(data@M)
  colnames(MATQ)=COLnames
  
  # do a MFA
  MATF=cbind(MATQ/sqrt((quantiles+1)))
  MFA.res= MFA(MATF,group=c(rep(quantiles+1,length(list.of.vars))),type=c(rep("c",length(list.of.vars))),
               ncp=3,name.group=colnames(data@M),graph=TRUE)
  # do Plots!!!!
  #plot Spanish-fan plot
  #plot individuals
  return(MFA.res)
}