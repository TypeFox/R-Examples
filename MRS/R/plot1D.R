#' Plot regions of the representative tree in 1D
#'
#' This function visualizes the regions of the representative tree 
#' of the output of the \code{\link{mrs}} function.
#' For each region the posterior probability of difference (PMAP)  or the effect size is plotted.
#' 
#' @param ans An \code{mrs} object.
#' @param type What is represented at each node. 
#' The options are \code{type = c("eff", "prob")}. 
#' Default is \code{type = "prob"}.
#' @param group If \code{type = "eff"}, which group effect size is used. 
#' Default is \code{group = 1}.
#' @param dim If the data are multivariate, \code{dim} is the dimension plotted. 
#' Default is \code{dim = 1}.
#' @param regions Binary vector indicating the regions to plot. 
#' The default is to plot all regions. 
#' @param legend Color legend for type. Default is \code{legend = FALSE}.
#' @param main Overall title for the plot.   
#' @references Soriano J. and Ma L. (2014). Multi-resolution two-sample comparison 
#' through the divide-merge Markov tree. \emph{Preprint}. 
#'  \url{http://arxiv.org/abs/1404.3753}
#' @export
#' @examples
#' set.seed(1)
#' p = 1
#' n1 = 200
#' n2 = 200
#' mu1 = matrix( c(0,10), nrow = 2, byrow = TRUE)
#' mu2 = mu1; mu2[2] = mu1[2] + .01
#' sigma = c(1,.1)
#' 
#' Z1 = sample(2, n1, replace=TRUE, prob=c(0.9, 0.1))
#' Z2 = sample(2, n2, replace=TRUE, prob=c(0.9, 0.1))
#' X1 = mu1[Z1] + matrix(rnorm(n1*p), ncol=p)*sigma[Z1]
#' X2 = mu2[Z2] + matrix(rnorm(n2*p), ncol=p)*sigma[Z1]    
#' X = rbind(X1, X2)
#' G = c(rep(1, n1), rep(2,n2))
#' 
#' ans = mrs(X, G, K=10)    
#' plot1D(ans, type = "prob")    
#' plot1D(ans, type = "eff")
plot1D <- function( ans,
                    type = "prob",
                    group = 1,
                    dim = 1,
                    regions = rep(1,length(ans$RepresentativeTree$Levels)),
                    legend = FALSE,
                    main = "default")
{
  if(class(ans)!="mrs")
  {
    print("ERROR: ans should be an MRS object.")
    return(0)
  }
  .pardefault <- par(no.readonly = T)
  if(legend)
  {
    layout(matrix(c(1,2), nrow=1),widths=c(.8,.2))    
    par(mar=c(5.1, 3.1, 4.1, 0.1))
  }
  else
  {
    layout(1)    
  }
  
  if(type == "prob")
  {
    names = round(ans$RepresentativeTree$AltProbs[which(regions==1)], digits=2)    
    col_range <- colorRampPalette(c("white","darkblue"))(100)
    col = col_range[ ceiling(names*99+ 1) ]
    if (main == "default")
      main = "PMAPs"
  }
  else if(type == "eff")
  {
    names = abs(ans$RepresentativeTree$EffectSizes[which(regions==1),group])    
    #     col_range <- c(colorRampPalette(c("red","white"))(50), 
    #                    colorRampPalette(c("white","dodgerblue"))(50))
    #     col = col_range[ ceiling(names/max(abs(names)+0.01)*100/2+50) ]
    col_range <- colorRampPalette(c("white","darkred"))(100)
    col = col_range[ ceiling( names/max(names + 0.01)*99 + 1) ]
    if (main == "default")
      main = paste("Effect Size Group", group)
  }
  else 
  {
    names = rep(1, length(ans$RepresentativeTree$AltProbs[which(regions==1)]))
    col = rep(NA, length(ans$RepresentativeTree$altProbs[which(regions==1)]))
    if (main == "default")
      main = ""
  }
  
  max.level = max(ans$RepresentativeTree$Levels[which(regions==1)]) 
  plot(1, type="n", yaxt='n', xlab="", ylab="", 
       xlim = ans$Data$Omega[dim,], ylim=c(0,max.level+.75), main = main)
  mtext("level", side=2, line = 2, at= (max.level+.75)/2 )
  axis(2, at= seq(0.4, max.level+0.4),
       labels=seq(max.level, 0, by=-1) )
  
  it = 1
  for (i in which(regions == 1) )
  {
    region = ans$RepresentativeTree$Regions[i,]
    rect(xleft=(region[2*dim-1]) , 
         ybottom= (max.level - ans$RepresentativeTree$Levels[i]) , 
         xright=region[2*dim] , 
         ytop=(max.level - ans$RepresentativeTree$Levels[i] + .75), col = col[it], border = col[it])
    it = it + 1
  }
  
  if(legend)
  {
    plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n")
    rect(xleft = 1.25, 
         ybottom = head(seq(1,2,1/100),-1), 
         xright = 1.75, 
         ytop = tail(seq(1,2,1/100),-1), 
         col=col_range, border = col_range )
    rect(1.25, 1, 1.75, 2.0)
    if(type == "prob")
      mtext(formatC(seq(0,1,.1), format = "f", digits = 1),side=2,at=seq(1,2.,by=.1),las=2,cex=1, line=0)
    if(type == "eff")
      mtext(formatC(seq( 0, max(names), length.out=11), format = "f", digits = 1),
            side=2,at=seq(1,2.,by=.1),las=2,cex=1, line=0)
    
  }
  
  par(.pardefault)   
  
  
}