#' Plot nodes of the representative tree 
#'
#' This function visualizes the representative tree of the output of the \code{\link{mrs}} function.
#' For each node of the representative tree, the posterior probability of difference (PMAP) or the effect size is plotted.
#' Each node in the tree is associated to a region of the sample space. 
#' All non-terminal nodes have two children nodes obtained by partitiing the parent region with a dyadic cut along a given direction.
#' The numbers under the vertices represent the cutting direction. 
#' 
#' @param ans A \code{mrs} object.
#' @param type What is represented at each node. 
#' The options are \code{type = c("eff", "prob")}.
#' @param group If \code{type = "eff"}, which group effect size is used. 
#' @param legend Color legend for type. Default is \code{legend = FALSE}.
#' @param main Main title. Default is \code{main = ""}.
#' @param node.size Size of the nodes. Default is \code{node.size = 5}.
#' @references Soriano J. and Ma L. (2014). Multi-resolution two-sample comparison 
#' through the divide-merge Markov tree. \emph{Preprint}. 
#'  \url{http://arxiv.org/abs/1404.3753}
#' @export
#' @note The package \pkg{igraph} is required.
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
#' plotTree(ans, type = "prob", legend = TRUE)
plotTree <- function(ans, type="prob", group = 1, legend = FALSE, main = "", node.size=5)
{
  if(class(ans)!="mrs")
  {
    print("ERROR: ans should be an mrs object.")
    return(0)
  }
  .pardefault <- par(no.readonly = T)
  if(legend)
  {
    layout(matrix(c(1,2), nrow=1),widths=c(.85,.15))    
    par(mar=c(1.1, 1.1, 3.1, 0.1))
  }
  else
  {
    layout(1)    
  }
  if(type == "prob")
  {
    box.size = (ans$RepresentativeTree$AltProbs)*node.size + .2*node.size
    name = round(ans$RepresentativeTree$AltProbs, digits=2)    
    col_range <- colorRampPalette(c("white","darkblue"))(100)
    col = col_range[ ceiling( name/max(name + 0.01)*99 + 1) ]
  }
  else if(type == "eff")
  {
    box.size = abs(ans$RepresentativeTree$EffectSizes[,group])/max(abs(ans$RepresentativeTree$EffectSizes[,group]))*node.size + .2*node.size
    name = abs(ans$RepresentativeTree$EffectSizes[,group])    
    col_range <- colorRampPalette(c("white","darkred"))(100)
    col = col_range[ ceiling( name/max(name + 0.01)*99 + 1) ]
  }
  else
  {
    box.size = 5
    name = NA
    col = "white"
  }
  
  M = list()
  it = 0
  for(i in 1:length(ans$RepresentativeTree$Levels))
  {
    for(j in 1:length(ans$RepresentativeTree$Levels))
    {
      if(ans$RepresentativeTree$Levels[i] == (ans$RepresentativeTree$Levels[j] -1)
         && ( 2*ans$RepresentativeTree$Ids[i] ==  ans$RepresentativeTree$Ids[j]  
              || (2*ans$RepresentativeTree$Ids[i]+1) ==  ans$RepresentativeTree$Ids[j] ))
      {
        it = it + 1
        M[[it]] = c(i,j)
      }
    }
  }
  M = unlist(M)
  vertex.label = ans$RepresentativeTree$Directions
  vertex.label[vertex.label==0] = NA
  vertex.label[ans$RepresentativeTree$Levels==max(ans$RepresentativeTree$Levels)]= NA
  
  G = igraph::graph(M, directed=FALSE)
  co <- igraph::layout.reingold.tilford(G, params=list(root=which(ans$RepresentativeTree$Levels==0))) 
  igraph::plot.igraph(G, layout=co, 
                      vertex.size = box.size, edge.label=NA,
                      vertex.color = col, vertex.label.color = "black",
                      vertex.label = vertex.label, asp=0, main=main, vertex.label.family = "Helvetica",
                      vertex.label.dist = 1.5, vertex.label.degree = pi/2)
  
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
      mtext(formatC(seq( 0, max(name) , length=11), format = "f", digits = 1),
            side=2,at=seq(1,2.,by=.1),las=2,cex=1, line=0)    
  }
  
  par(.pardefault)   
  
}
