# ###############################################################################
# \name{plot.gnbp}
# 
# \title{Plot a Genotype-Phenotype Network}
# 
# \description{
#   Plot method for gnbp objects. Plots genotype-phenotype network in which evidence has been absorbed and propagated and maps the beliefs or Jeffrey's signed information onto the network.
# }
# 
#   \usage
#   \method{plot}{gnbp}(x, y="JSI",col.palette,col.length = 100, ncol = 1, 
#                       fontsize=10, fontcolor="black",...)
#   }
#  
# \arguments{
#   \item{x}{
#   An object of class gnbp
#   }
#   \item{y}{
#   A character string. Valid options are \code{"JSI"} (default) 
# or \code{"belief"} for Conditional Gaussian network. 
# For Discrete Bayesian networks, the argument will be ignored 
# and the phenotype states with maximum probability will be plotted.
#   }
#   \item{col.palette}{A list of character strings. 
# For Conditional Gaussian networks, a list of 6 elements 
# specifying colors for colormap.All 6 elements should be 
# character strings specifying the colour for 
#   \code{pos_high}= high end of gradient of positive values (default = "red")
#   \code{pos_low}i=low end of gradient of positive values (default = "wheat1")
#   \code{neg_high}=high end of gradient of positive values (default = "cyan")
#   \code{neg_low}=low end of gradient of positive values (default = "blue")
#   \code{dsep_col}= \emph{d}-separated nodes (default = "white")
#   \code{qtl_col}= discrete nodes (QTLs) (default = "grey")
#   \code{node_abs_col}= nodes for which evidence has 
# been absorbed (default = "palegreen2")
#   
#   For Discrete Bayesian networks, a list of 4 elements 
# specifying colors for colormap should be specified. All 4 elements 
# should be character strings specifying the colour for 
#   \code{col_nodes}- a vector of colors for phenotype states
# should be specified. The length of the vector should be
# equal to the maximum number of phenotype states possible.
#   \code{dsep_col}= \emph{d}-separated nodes (default = "white")
#   \code{qtl_col}= discrete nodes (QTLs) (default = "grey")
#   \code{node_abs_col}= nodes for which evidence has been
# absorbed (default = "palegreen2")
#   
#   }
#   \item{col.length}{
#   a positive integer (default = 100) specifying the 
# resolution of the colormap (number of colors).
#   }
#   \item{ncol}{a positive integer specifying the column number 
# of JSI / belief / FC to plot. By default, the first column will be plot}
#   \item{fontsize}{fontsize for node labels}
#   \item{fontcolor}{fontcolor for node labels}
#   \item{...}{further arguments to the function \code{\link{plot}}. 
# These will be ignored}
#   }
# 
# \value{
#   \code{x} is invisibly returned
# }
#############################################################################

plot.gnbp=function(x, y="JSI",col.palette,col.length = 100, ncol = 1, fontsize=10, fontcolor="black",...)

{
                  
  requireNamespace("RHugin") || warning("Package not loaded: RHugin");
  
  ## extract data
  Data<-RHugin::get.cases(x$gp)
  
  ## get node attributes
  class_nodes=x$gp_nodes
  
  ## get network type
  type=x$gp_flag
  
  ## get d-connected nodes
  dnodes<-RHugin::get.dconnected.nodes(x$gp,x$node)
  dnodes<-class_nodes[match(setdiff(class_nodes[match(dnodes,class_nodes),1],x$node),class_nodes),]
  dnodes<-matrix(dnodes,ncol=4)
  colnames(dnodes)<-colnames(class_nodes)
  
  ## convert to graphNEL object for use with Rgraphviz
  BNgraph<-RHugin::as.graph.RHuginDomain(x$gp)
  
  ## set node attributes
  z<-graph::nodes(BNgraph)
  names(z)<-graph::nodes(BNgraph)
  
  nAttrs <- list()
  eAttrs <- list()
  attrs<-list()
  nAttrs$label<-z
  
  attrs$node$fontsize<-fontsize
  attrs$node$fontcolor<-fontcolor
  attrs$node$fixedsize<-F
  attrs$node$height<-1.2

  
  nAttrs$shape<-rep("ellipse",length(z))
  names(nAttrs$shape)<-graph::nodes(BNgraph)
  nAttrs$shape[class_nodes[which(class_nodes[,"type"]=="geno"),1]]<-"box"
  
  if (type == "cg")
  {
  
      if(missing(col.palette))
      col.palette<-list(pos_high="red", pos_low= "wheat1", 
                        neg_high="cyan", neg_low = "blue",
                        dsep_col="white",qtl_col="grey",node_abs_col="palegreen2") 
      
      
      if (y=="JSI")
      {score=x$JSI[,ncol]} else
        {if (y=="belief")
          {score=x$belief$pheno$mean[,ncol]} else
            {warning("invalid option for the argument 'y': plotting JSI ....")
             score=x$JSI[,ncol]}}
      
      X=length(dnodes[which(dnodes[,4]=="geno"),1])
      nAttrs$fillcolor<-rep(col.palette$dsep_col,length(z))
      names(nAttrs$fillcolor)<-graph::nodes(BNgraph)
      
        
        score_no<-seq(floor(range(score)[1]),ceiling(range(score)[2]),length.out=col.length)
        score_range<-cbind(score_no[-length(score_no)],score_no[-1])
        score_int<-findInterval(score,sort(score_no))
      
        if (!is.null(x$node) || !is.null(score))
        {
          if (length(which(score_no<0))!=0)
           {colormap1<-cscale(1:(length(which(score_no<0))-1),palette=seq_gradient_pal(low=col.palette$neg_low,high=col.palette$neg_high))}   
          else
            {colormap1<-NULL}
          if (length(which(score_no>0))!=0)
            colormap2<-cscale(1:(length(which(score_no>0))-1),palette=seq_gradient_pal(low=col.palette$pos_low,high=col.palette$pos_high))   
          else
            colormap2<-NULL
          
          colormap_network<-c(colormap1,colormap2)[score_int]
          
          for (i in 1:length(x$node))
          {
            colormap_network<-c(colormap_network,col.palette$node_abs_col)
          }
          
          colormap_network<-c(colormap_network, rep(col.palette$qtl_col,X))
          
          names(colormap_network)=c(c(names(score),x$node),dnodes[which(dnodes[,4]=="geno"),1])
          nAttrs$fillcolor<-colormap_network
          
#           non_dnodes<-setdiff(class_nodes[,1],names(colormap_network))
          colormap<-c(c(colormap1,col.palette$dsep_col),colormap2)
          
          
          layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,2,0), 4, 6, byrow = F))
          
#           layout(matrix(c(1,2), 1, 2, byrow = F),c(4,1),c(1,1))

          plot(BNgraph, nodeAttrs=nAttrs,attrs=attrs) 
    
          scale = (length(colormap)-1)/(round(max(score_no))-round(min(score_no)))
          ticks=seq(round(min(score_no)),0, len=3)
          ticks<-ticks[-3]
          ticks=c(ticks,seq(0,round(max(score_no)),len=3))
          
          plot(c(0,10), c(round(min(score_no)),round(max(score_no))), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
          axis(2, pos=3,ticks, las=1,font=1)
          for (i in 1:(length(colormap)-1)) 
            {
              y = (i-1)/scale + round(min(score_no))
              rect(3.05,y,4.05,y+(1/scale), col=colormap[i], border=NA)
            }
        
       }
  }
  
  
  if (type == "db")
  {
    if(missing(col.palette))
      col.palette<-list(col_nodes=c("Red","Blue","Orange","Yellow","Magenta","Brown"),
                        dsep_col="white",qtl_col="grey",node_abs_col="palegreen2") 
   
    ## generate colors for node states
    colormap<-dscale(factor(x$FC$pheno_state[,ncol]),manual_pal(col.palette$col_nodes))
    
    ## add color for absorbed node
    for (i in 1:length(x$node))
    {
      colormap_network<-c(colormap,col.palette$node_abs_col)
    }
    

    ## add color for d-connected qtl nodes
    X=length(dnodes[which(dnodes[,4]=="geno"),1])
    colormap_network<-c(colormap_network, rep(col.palette$qtl_col,X))
    
    ## set node names to colormap_network
    names(colormap_network)=c(c(names(colormap),x$node),dnodes[which(dnodes[,4]=="geno"),1])
    
    ## set default color of all nodes
    nAttrs$fillcolor<-rep(col.palette$dsep_col,length(z))
    names(nAttrs$fillcolor)<-graph::nodes(BNgraph)
    
    ## set color for all nodes
    nAttrs$fillcolor<-colormap_network
    
    layout(matrix(c(1,2), 1, 2, byrow = F),c(4,1),c(1,1))
    
    plot(BNgraph, nodeAttrs=nAttrs,attrs=attrs)
    
    node_leg=c(paste("state",c(1:max(dnodes[,"levels"]))))
    names(node_leg)=col.palette$col_nodes[1:length(node_leg)]
    
    leg=c(node_leg,"d-separated","d-connected","absorbed node")
    leg.col=c(names(node_leg),col.palette$dsep_col,col.palette$qtl_col,col.palette$node_abs_col)  
    
    pch=c(rep(21,length(node_leg)),22,22,21)
    
    plot(c(1:length(leg)),type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='',bg=leg.col,col="black",pch=pch)
    legend("center",legend=leg,col=par("col"),pt.bg=leg.col,pch=pch,bty='n',pt.cex=1.5,cex=0.8)
    
   
    
  }
  
  invisible(x)
  
}
