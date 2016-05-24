#' Plot regions of the representative tree in 2D
#'
#' This function visualizes the regions of the representative tree 
#' of the output of the \code{\link{mrs}} function.
#' 
#' @param ans An \code{mrs} object.
#' @param type Different options on how to visualize the rectangular regions. 
#' The options are \code{type = c("eff", "prob", "empty", "none")}. 
#' Default is \code{type = "prob"}.
#' @param data.points Different options on how to plot the data points. 
#'  The options are \code{data.points = c("all", "differential", "none")}. 
#'  Default is \code{data.points = "all"}. 
#' @param background Different options on the background. 
#'  The options are \code{background = c("smeared", "none") }.
#' @param group If \code{type = "eff"}, which group effect size is used. 
#'  Default is \code{group = 1}.
#' @param dim If the data are multivariate, 
#'  \code{dim} are the two dimensions plotted. Default is \code{dim = c(1,2)}.
#' @param levels Vector with the level of the regions to plot. 
#'  The default is to plot regions at all levels.
#' @param regions Binary vector indicating the regions to plot. 
#'  The default is to plot all regions. 
#' @param legend Color legend for type. Default is \code{legend = FALSE}.
#' @param main Overall title for the legend.  
#' @references Soriano J. and Ma L. (2014). Multi-resolution two-sample comparison 
#' through the divide-merge Markov tree. \emph{Preprint}. 
#'  \url{http://arxiv.org/abs/1404.3753}
#' @export
#' @examples
#' set.seed(1)
#' p = 2
#' n1 = 200
#' n2 = 200
#' mu1 = matrix( c(9,9,0,4,-2,-10,3,6,6,-10), nrow = 5, byrow=TRUE)
#' mu2 = mu1; mu2[2,] = mu1[2,] + 1
#'  
#' Z1 = sample(5, n1, replace=TRUE)
#' Z2 = sample(5, n2, replace=TRUE)
#' X1 = mu1[Z1,] + matrix(rnorm(n1*p), ncol=p)
#' X2 = mu2[Z2,] + matrix(rnorm(n2*p), ncol=p)
#' X = rbind(X1, X2)
#' colnames(X) = c(1,2)
#' G = c(rep(1, n1), rep(2,n2))
#'   
#' ans = mrs(X, G, K=8)
#' plot2D(ans, type = "prob", legend = TRUE)
#'   
#' plot2D(ans, type="empty", data.points = "differential", 
#'  background = "none") 
#'   
#' plot2D(ans, type="none", data.points = "differential", 
#'  background = "smeared", levels = 4)       
plot2D = function(  ans, 
                    type = "prob", 
                    data.points = "all", 
                    background = "none", 
                    group = 1, 
                    dim = c(1,2),
                    levels = sort(unique(ans$RepresentativeTree$Levels)),
                    regions = rep(1,length(ans$RepresentativeTree$Levels)),
                    legend = FALSE,
                    main = "default"
                ) 
{ 
  if(class(ans)!="mrs")
  {
    print("ERROR: ans should be an mrs object.")
    return(0)
  }
  
  if(ncol(ans$Data$X)<2)
  {
    print("ERROR: p>=2")
    return(0);
  }
  
  .pardefault <- par(no.readonly = T)
  
  size.tmp = ceiling(sqrt(length(levels)))
  if( size.tmp * (size.tmp -1) == length(levels) )
  {
    mat = matrix( seq(1,  size.tmp*(size.tmp-1)) , ncol=size.tmp, byrow=TRUE  )
    if(!legend)
      layout( mat = mat )
    else
    {
      mat = cbind(mat, size.tmp*(size.tmp - 1) + seq(1, size.tmp-1) )
      layout( mat = mat, widths= c(rep(0.85, ncol(mat)-1)/(ncol(mat)-1), 0.15) )      
    }    
  }
  else if( size.tmp * (size.tmp -1) > length(levels) )
  {
    mat = matrix( seq(1,  size.tmp*(size.tmp-1)) , ncol=size.tmp, byrow=TRUE  )
    layout( mat = mat )
  }
  else if( size.tmp* size.tmp == length(levels) )
  {
    mat = matrix( seq(1,  size.tmp^2), ncol=size.tmp, byrow=TRUE  )
    if(!legend)
      layout( mat )
    else
    {
      mat = cbind(mat, size.tmp^2 + seq(1, size.tmp) )
      layout( mat = mat, widths= c(rep(0.85, ncol(mat)-1)/(ncol(mat)-1), 0.15) )      
    }
  } 
  else
  {
    mat = matrix( seq(1,  size.tmp*size.tmp) , ncol=size.tmp, byrow=TRUE  ) 
    layout(mat = mat)
    #     if(!legend)
    #       layout( mat = mat )
    #     else
    #     {
    #       mat = cbind(mat, size.tmp*(size.tmp - 1) + seq(1, size.tmp-1) )
    #       layout( mat = mat )      
    #     }
  }
  
  par(mar = c(3.1, 3.1, 3.1, 1.1))
  
  
  if(type == "prob")
  {
    names = round(ans$RepresentativeTree$AltProbs, digits=2)    
    col_range <- colorRampPalette(c("white","darkblue"))(100)
    col = col_range[ ceiling( names*99 + 1) ]
    if(main == "default")
      main = "PMAPs"
    
  }
  else if(type == "eff")
  {
    names = abs(ans$RepresentativeTree$EffectSizes[,group])
    col_range <- colorRampPalette(c("white","darkred"))(100)
    col = col_range[ ceiling( names/max(names + 0.01)*99 + 1) ]
    if(main == "default")
      main = paste("Eff. Size \n Group", group)
  }
  else 
  {
    names = rep(1, length(ans$RepresentativeTree$AltProbs))
    col = rep(NA, length(ans$RepresentativeTree$AltProbs))
    if(main == "default")
      main = ""
  }
  
  
  for( i in levels )
  {
    if(background == "smeared")
    {
      Lab.palette <- colorRampPalette(c("white", "black"), space = "rgb")
      smoothScatter(ans$Data$X[,dim], nrpoints = 0, colramp = Lab.palette, 
                    xlim = ans$Data$Omega[dim[1],], ylim=ans$Data$Omega[dim[2],],
                    main = paste("level",substitute(l, list(l = i))),
                    xlab = colnames(ans$Data$X)[dim[1]], ylab = colnames(ans$Data$X)[dim[2]])
    }
    else
    {
      plot(1, type="n", xlab = colnames(ans$Data$X)[dim[1]], ylab = colnames(ans$Data$X)[dim[2]],
           xlim = ans$Data$Omega[dim[1],], ylim=ans$Data$Omega[dim[2],],   
           main = paste("level",substitute(l, list(l = i))) )
    }
    
    
    nodes = which(ans$RepresentativeTree$Levels==i & regions==1)
    
    if( type =="eff")
      nodes = nodes[sort.int( abs(ans$RepresentativeTree$EffectSizes[nodes,group]), 
                              index.return= TRUE )$ix] 
    else if(type=="prob")
      nodes = nodes[sort.int( ans$RepresentativeTree$AltProbs[nodes], index.return= TRUE )$ix]
    
    for( j in nodes )
    {
      region = ans$RepresentativeTree$Regions[j,]
      points_idx = ans$RepresentativeTree$DataPoints[[j]]
      
      if(type=="empty")
        rect(xleft=region[2*dim[1]-1] , ybottom=region[2*dim[2]-1] , 
             xright=region[2*dim[1]] , ytop=region[2*dim[2]], col=NA, border = 1 )
      else if(type!="none")
        rect(xleft=region[2*dim[1]-1] , ybottom=region[2*dim[2]-1] , 
             xright=region[2*dim[1]] , ytop=region[2*dim[2]], col=col[j], border = col[j] )
      
      if(data.points=="all")
        points(ans$Data$X[,dim], col=(ans$Data$G+1), pch = ans$Data$G )
      else if(data.points=="differential")
      {
        points(ans$Data$X[points_idx,dim], col=(ans$Data$G[points_idx,]+1), 
               pch = ans$Data$G[points_idx,] )
      }
      
    }
    
    if(legend && (i == levels[1]) && (data.points!="none") )
    {
      group_lab = rep(NA,max(ans$Data$G))
      for(gg in 1:max(ans$Data$G))
        group_lab[gg] = paste("sample", gg)
      legend("bottomright", legend=group_lab, col=seq(2,max(ans$Data$G)+1), 
             pch = seq(1,max(ans$Data$G)), bg = "white" ) 
      
    }
  }
  
  if(legend)
  {
    plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n")
    title(main, adj = 0)
    if( (type == "prob") || (type == "eff") )
    {
      # mtext(main, side = 3, at = 1, cex = 1)
      rect(xleft = 1, 
           ybottom = head(seq(1,2,1/100),-1), 
           xright = 1.2, 
           ytop = tail(seq(1,2,1/100),-1), 
           col=col_range,  border=col_range )
      rect(1, 1, 1.2, 2.0)
    }
    
    if(type == "prob")
      mtext( format( round(seq(0,1,by=0.2), digits=1), nsmall=1) ,side=2,at=seq(1,2,by=.2),las=2,cex=1, line=0)
    if(type == "eff")
      mtext( format(round(seq( 0, max(names) , length=5), digits=1), nsmall=1),
             side=2,at=seq(1,2,length=5),las=2,cex=1, line=0)
    
  }
  par(.pardefault)   
}