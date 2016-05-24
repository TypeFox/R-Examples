#####################################################################################
# \name{fit.gnbp}

# \title{Fit a Conditional Gaussian Bayesian Network to QTL data}

# \description{Learn the structure of a genotype-phenotype network from 
#   quantitative trait loci (QTL) data and the conditional probability table 
#   for each node in the network using the PC algorithm and the EM algorithm 
#   implemented in the RHugin package. 
# }
# 
# \usage{
#   fit.gnbp(geno,pheno,constraints,learn="TRUE",edgelist,type ="cg",
#            alpha=0.001,tol=1e-04,maxit=0)

# }
# 
# \arguments{
#   \item{geno}{a data frame of column vectors of class factor 
#   (or one that can be coerced to that class) and non-empty column names.  
#   }
#   \item{pheno}{a data frame of column vectors of class numeric if 
# \code{type = "cg"} or class factor if \code{type = "db"} and non-empty column names. 
#   }
#   \item{constraints}{an optional list of constraints on the edges for 
# specifying required and forbidden edges. See details.
#   }
#   \item{learn}{a boolean value. If TRUE (default), the network structure 
# will be learnt using the PC algorithm in RHugin package. 
# If FALSE, only conditional probabilities will be learnt.
#   }
#   \item{edgelist}{a list of edges to be provided if learn == FALSE.
#   }
#   \item{type}{specify the type of network. \code{"cg"} for 
# \code{Conditional Gaussian} (default) and \code{"db"} for \code{Discrete Bayesian}.
#   }
#   \item{alpha}{a single numeric value specifying the significance level 
# (for use with RHugin). Default is 0.001.
#   }
#   \item{tol}{a positive numeric value (optional) specifying the tolerance 
# for EM algorithm to learn conditional probability tables (for use with RHugin).
# Default value is 1e-04. See \code{learn.cpt} for details.
#              
#   }
#   \item{maxit}{a positive integer value (optional) specifying the 
# maximum number of iterations of EM algorithm to learn 
# conditional probability tables (for use with RHugin). See \code{learn.cpt} for details.
#   }
# 
#   \value{
#     Returns an object of class "gpfit". 
#     \item{gp}{a pointer to a compiled RHugin domain. 
#       There is a cpt table associated with each node in the network.}
#     \item{gp_nodes}{a data frame containing information about 
#       nodes for internal use with other functions.}
#     \item{gp_flag}{a character string specifying the type of 
#         network (\code{Conditional Gaussian} or \code{Discrete Bayesian})}
#     
#   }

###################################################################################
## Learn Bayesian Network Structure
## from QTL data
###################################################################################
fit.gnbp=function(geno,pheno,constraints,learn="TRUE",edgelist,type ="cg",
                  alpha=0.001,tol=1e-04,maxit=0)

  {
  
    requireNamespace("RHugin") || warning("Package not loaded: RHugin");
  
    ## Geno Class Check ##
    for(i in 1:dim(geno)[2])
      if (class(geno[,i])!="factor")
        {warning("column vectors of 'geno' are not of class factor. converting to factor...")
         geno[,i]<-as.factor(geno[,i])
        }
    
    ## Pheno Class Check ##
      class_pheno=lapply(pheno,class)
     
    if (type=="cg")
      { 
       X<-which(class_pheno=="numeric")
       if(length(X)!=dim(pheno)[2])
        stop("column vectors of 'pheno' should be of class numeric.")
      }
    
    if (type=="db")
    { 
      X<-which(class_pheno=="factor")
      if(length(X)!=dim(pheno)[2])
        stop("column vectors of 'pheno' should be of class factor.")
    }
    
    Data=cbind(pheno,geno)
    
    #####################################
    ## Create RHugin domain
    ## Specify nodes and constraints
    #####################################
     
    ## Create a RHugin domain
    network<-RHugin::hugin.domain()
    Nnodes<-length(colnames(Data))
    
      
    ## Determine type of nodes (discrete or continuous)
      class_nodes=matrix(nrow=dim(Data)[2],ncol=3)

    for (i in 1:dim(Data)[2])
    {
        class_nodes[i,]=c(colnames(Data)[i],class(Data[,i]),length(levels(Data[,i])))
    }
    
    class_nodes=cbind(class_nodes,t(cbind(t(rep("pheno",dim(pheno)[2])),t(rep("geno",dim(geno)[2])))))
    
    colnames(class_nodes)=c("node","class","levels","type")
    
     
    Xpheno=which(class_nodes[,"type"]=="pheno")
    Xgeno=which(class_nodes[,"type"]=="geno")

     
    ## Add nodes
    for (node in 1:length(Xpheno))
    {
      if(class_nodes[Xpheno[node],"class"]=="numeric")
        RHugin::add.node(network,class_nodes[Xpheno[node],"node"],kind="continuous") 
      if(class_nodes[Xpheno[node],"class"]=="factor")
        RHugin::add.node(network,class_nodes[Xpheno[node],"node"],states=levels(Data[,class_nodes[Xpheno[node],"node"]]))    
    }

    for (node in 1:length(Xgeno))
      RHugin::add.node(network,class_nodes[Xgeno[node],"node"],states=levels(Data[,class_nodes[Xgeno[node],"node"]]))

    
    ## Set cases
    RHugin::set.cases(network,Data)
     
    if (learn=="TRUE")
    {
 
    # Disallow interaction between discrete nodes
    qtl_constraints<-.GenConstraints(class_nodes)
    
      if(!missing(constraints))
      {
      directed.required<- constraints$directed$required
      directed.forbidden<- constraints$directed$forbidden
        
        if (is.null(constraints$undirected$required))
          {undirected.required<-qtl_constraints$undirected$required} else 
          {idx <- unique(c(names(constraints$undirected$required), names(qtl_constraints$undirected$required)))
           undirected.required<-setNames(mapply(c, constraints$undirected$required[idx], qtl_constraints$undirected$required[idx]), idx)} 
        if (is.null(constraints$undirected$forbidden))
           {undirected.forbidden<-qtl_constraints$undirected$forbidden} else 
           {idx <- unique(c(names(constraints$undirected$forbidden), names(qtl_constraints$undirected$forbidden)))
            undirected.forbidden<-setNames(mapply(c, constraints$undirected$forbidden[idx], qtl_constraints$undirected$forbidden[idx]), idx)}
              
              all_constraints=list(directed=list(required=directed.required,forbidden=directed.forbidden),
                                   undirected=list(required=undirected.required,forbidden=undirected.forbidden))
              
      } else
          all_constraints<-qtl_constraints

  
    ####################################
    ## Learn structure and cpt
    ####################################
     
    RHugin::learn.structure(network,alpha=alpha,constraints=all_constraints)
    
    }else
      
    {
      if (missing(edgelist))
        stop("if learn == TRUE, edgelist must be provided.")
      
      for (i in 1:length(edgelist))
      {
        edges<-unlist(edgelist[[i]])
        RHugin::add.edge(network,edges[,2],edges[,1])
      }

    }
  
    ##Add experience table to all nodes
    for (node in 1:Nnodes)
      RHugin::get.table(network,colnames(Data)[node],type="experience")
    
    ## Compile network and learn cpt
    RHugin::compile(network)
    RHugin::learn.cpt(network,tol=tol,maxit=maxit)
        
    gpfit<-list(gp=network,gp_nodes=class_nodes,gp_flag=type)
    
    class(gpfit)<-"gpfit"
    
    return(gpfit)   
   
  }  

.GenConstraints=function(class_nodes)
{
  X=which(class_nodes[,"type"]=="geno")
  blackL<-matrix(0,nrow=choose(length(X),2) ,ncol=2)
  k=1
  for (i in 1:(length(X)-1))
    for (j in (i+1):length(X))
    {
      blackL[k,]=cbind(class_nodes[X[i],1],class_nodes[X[j],1])
      k=k+1
    }
  

  undirected.forbidden <- vector("list", nrow(blackL))
  
  for (i in 1:nrow(blackL))
  {
    undirected.forbidden[[i]] <- blackL[i,]
  }
  
  undirected <- list(required = NULL,forbidden = undirected.forbidden)
  
  # put into constraints list
  qtl_constraints <- list(directed = NULL, undirected = undirected)
  return(qtl_constraints)
}




