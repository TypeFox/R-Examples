################################################################################
# \name{gen.evidence}
# 
# \title{Generate a sequence of evidence for a continuous
#        node in a conditional gaussian bayesian network.
# }
# \description{
#   The evidence for a node in an RHugin domain is generated 
#   as a linear sequence within the specified standard deviation 
#   from the marginal mean of the node. The evidence can be 
#   given as an input to \link{absorb.gnbp}
# }
# 
# \usage{
#   gen.evidence(gpfit, node, std = 2, length.out = 10, std.equal = TRUE)
# }
# 
# \arguments{
#   \item{gpfit}{an object of class "gpfit" obtained by using fit.gnbp
#   }
#   \item{node}{
#     a character string specifying the name of a continuous node in the domain
#   }
#   \item{std}{
#     a numeric value specifying the number of standard deviations 
#     of marginal distribution within which the evidence is generated. 
#     A numeric vector of length = number of nodes, 
#     must be specified when std.equal=FALSE. 
#   }
#   \item{length.out}{
#     a positive integer giving the desired length of the sequence.
#   }
#   \item{std.equal}{a logical value indicating whether 
#                    same number of standard deviations should 
#                    be used to generate evidence for all nodes. Default is TRUE. 
#   }
# }
# 
# \value{
#   A matrix of evidence for each specified node
# }
###############################################################################

gen.evidence=function(gpfit,node,std=2,length.out=10,std.equal=TRUE)
  
{  
  
  requireNamespace("RHugin") || warning("Package not loaded: RHugin");
  
  RHugin::retract(gpfit$gp)
  
  ## Get node attributes
  Data<-RHugin::get.cases(gpfit$gp)

  class_nodes<-gpfit$gp_nodes
  
  ## Determine type of nodes (discrete or continuous)
    if (is.element("factor",class_nodes[match(node,class_nodes),2]))
        stop("One or more nodes specified is not continuous")
  

  ## calculate marginals
  marg_node_abs<-RHugin::get.marginal(gpfit$gp,node)
  mean1<-marg_node_abs$mean
  var1<-diag(matrix(unlist(marg_node_abs$cov),nrow=length(node)))

  evidence=matrix(nrow=length(node),ncol=length.out)
  
  if (std.equal==FALSE & length(std)!=length(node))
    stop("length of std is not equal to the number of nodes") 
    

#      warning("std.equal is set to true. using the first element of the vector std")
#    
     std<-rep(std[1],length(node)) 
    
     for (i in 1:length(node))
       {
         a<-mean1[i]-std[i]*sqrt(var1[i])
         b<-mean1[i]+std[i]*sqrt(var1[i])
         evidence[i,]=seq(a,b,length.out=length.out)
       }

rownames(evidence)=node


  return(evidence)
}
