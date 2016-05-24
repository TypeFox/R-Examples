#' Analyse composite records by a Kruskal-Wallis ANOVA
#' 
#' The function applies a Kruskal-Wallis ANOVA on binned data issued from a
#' "pfComposite" object (of directly on "pfTransform" objects), in order to
#' test the difference in biomass burning activity between different time
#' periods.
#' 
#' 
#' @param data An object returned by \code{\link[paleofire]{pfComposite}} or
#' \code{\link[paleofire]{pfTransform}}.
#' @param p.adj Method for adjusting p values (see
#' \code{\link[stats]{p.adjust}}). Includes: "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr" and "none" (default).
#' @param alpha Numeric, confidence level.
#' @param bins Numeric, bins to use if a "pfTransform" object is provided.
#' @param verbose Logical, verbose or not...
#' @return A "pfKruskal" object containing multiple comparison results.
#' @author O. Blarquez
#' @seealso \code{\link{plot.pfKruskal}},\code{\link[agricolae]{kruskal}}
#' @examples
#' 
#' ## Composite charcoal record for Western Boreal North America:
#' ID=pfSiteSel(id_region=="WNA0" & l12==1)
#' plot(ID)
#' ## Transform data
#' res3=pfTransform(ID,method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,4000))
#' 
#' ## Composite
#' comp=pfComposite(res3,bins=seq(from=-500,to=12500,by=1000))
#' plot(comp)
#' 
#' ## Kruskal Wallis Anova
#' comparison=pfKruskal(comp)
#' 
pfKruskal=function(data,p.adj="none",
                   alpha=0.05,bins=NULL,
                   verbose=TRUE){
  
  if (!requireNamespace("agricolae", quietly = TRUE)) {
  install.packages("agricolae")
  }
  
  if(class(data)=="pfComposite"){
    data=cbind(x=rep(data$BinCentres,ncol(data$BinnedData)),y=c(data$BinnedData))
    data=data.frame(data)
  }
  
  if(class(data)=="pfTransform"){
    data=data.frame(na.omit(cbind(y=c(data$TransData), x=c(data$Age))))
    if(is.null(bins)) bins=seq(-500,12500,1000)
    xx=as.numeric(cut(data$x,breaks=bins))
    
    hbins=bins+max(diff(bins)/2)
    hbins=data.frame(cbind(1:(length(hbins)-1),hbins[1:(length(hbins)-1)]))
    id <- with(hbins, X2[match(xx,X1)])
    
    data=data.frame(na.omit(cbind(x=id,y=data$y)))
  }
  
  ## Post-hoc  
  data$x=as.factor(data$x)
  comparison=agricolae::kruskal(data$y, data$x,alpha = alpha, p.adj=p.adj, group=TRUE, main = NULL)
  
  ## Data prep
  signifL=as.data.frame(comparison$group)
  signifL$M=as.character(signifL$M)
  
  signifL$trt=as.numeric(paste(comparison$group$trt))
  signifL=signifL[order(signifL$trt) , ]

  ## Pos of letters
  
  test=boxplot(data$y~data$x)$stats[5,]
  garbage <-dev.off()
  
  signifL=data.frame(cbind(signifL,lpos=test))
  signifL$trt=as.factor(signifL$trt)
  colnames(signifL)=c("x","means","lab","lpos")
  
  # X reverse
  data$x <- factor(data$x, levels = rev(levels(data$x)))
  
  if(verbose==TRUE) print(comparison)
  
  out=list(letters=signifL,data=data,summary=comparison)
  class(out)="pfKruskal"
  return(out)
}





#' Plot a "pfKruskal" object.
#' 
#' Plot a "pfKruskal" object using boxplots and showing significant differences
#' between the periods using letters.
#' 
#' If two periods share the same letter their rank (median) is not
#' significantly different at the confidence level specified by \code{alpha}.
#' If not, equality could be rejected at the confidence level specified by
#' \code{alpha}.
#' 
#' @method plot pfKruskal
#' @export
#' @param x An object returned by \code{\link[paleofire]{pfKruskal}}.
#' @param trend Logical, show trend using linear regression?
#' @param outliers Logical, show outliers?
#' @param xlim Numeric, x axis limits.
#' @param ylim Numeric, y axis limits.
#' @param \dots ...
#' @return Return a ggplot2 "gg" object.
#' @author O. Blarquez
#' @seealso \code{\link[paleofire]{pfKruskal}}
#' @examples
#' 
#' ## Composite charcoal record for Western Boreal North America:
#' ID=pfSiteSel(id_region=="WNA0" & l12==1)
#' plot(ID)
#' ## Transform data
#' res3=pfTransform(ID,method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,4000))
#' 
#' ## Composite
#' comp=pfComposite(res3,bins=seq(from=-500,to=12500,by=1000))
#' plot(comp)
#' 
#' ## Kruskal Wallis Anova
#' comparison=pfKruskal(comp)
#' 
#' plot(comparison)
#' 
#' # p=plot(comparison)
#' # require(ggplot2)
#' # p+ggtitle("my title")
#' 
plot.pfKruskal=function(x,trend=FALSE,outliers=FALSE, xlim = NULL, ylim = NULL, ...){  
  
  y<-lpos<-lab<-NULL
  x$data=na.omit(x$data)
  ### GGplot2 boxplot
  p=ggplot(data=x$data,aes(x=x,y=y))
  if(outliers==TRUE) p=p+geom_boxplot(fill="grey90",outlier.size = 1,show_guide = FALSE)
  if(outliers==FALSE) p=p+geom_boxplot(fill="grey90",outlier.size = 0,show_guide = FALSE)
  p=p+theme_bw(base_size=18)+xlab("Age (cal BP)") + ylab("Composite")+
    geom_text(aes(x=x, y=lpos+0.1, label = lab), data = x$letters,vjust=0,size=5)+
    coord_cartesian(xlim = xlim, ylim = ylim)
  # geom_text(data=data.frame( x = 15, y= 2.7), map=aes(x=x, y=y), label = w) 
  
  if(trend==TRUE)
    # p=p+geom_text(data=data.frame( x = 1, y= 3), map=aes(x=x, y=y), label = label) 
    p=p+geom_smooth(method="lm", se=FALSE,color="black",aes(group=1))
  p
  return(p)  
}
