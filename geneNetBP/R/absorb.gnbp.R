###############################################################################
# \name{absorb.gnbp}
# 
# \title{Absorb evidence and infer a genotype-phenotype network}
# 
# \description{
#   Absorb a single piece or a spectrum of evidence for one or more continuous nodes 
#  in a compiled RHugin domain, obtain the updated beliefs 
#  and the Jeffrey's signed information. }
# 
#   \usage{
#   absorb.gnbp(gpfit, node, evidence)
#   }
# 
#   \arguments{
# 
#   \item{gpfit}{
#   an object of class "gpfit". Output from \code{\link{fit.gnbp}}. 
#   }
#   
#   \item{node}{
#   a character vector specifying the names of the nodes for 
# which the evidence is to be absorbed.
#   }
#   \item{evidence}{
#   a matrix or a numeric vector of evidence. number of rows of the 
# matrix or the length of the vector should be equal to the length of node. 
#   }
#   }
#   
#   \value{
#   absorb.gnbp returns an object of class "gnbp". The functions 
# summary and print can be used for objects of class "gnbp". 
# An object of class "gnbp" is a list containing the following components
#   
#   \item{gp}{an RHugin domain that is triangulated,
# compiled and with the latest absorbed evidence propagated }
#   \item{gp_flag}{type of network}
#   \item{node}{a character vector specifying the nodes 
# for which evidence has been absorbed}
#   \item{marginal}{a list of marginal probabilities for 
# phenotypes (\code{pheno}) and genotypes (\code{geno})}
#   \item{belief}{a list of updated beliefs for phenotypes 
# (\code{pheno}) and genotypes (\code{geno})}
#   \item{JSI}{a matrix of Jeffrey's signed information if 
# network is \code{Conditional Gaussian}, otherwise \code{NULL} 
# if network is \code{Discrete Bayesian}}
# \item{FC}{a list of two. a matrix \code{FC} of fold changes and 
# a matrix \code{pheno_state} of phenotype node beliefs - state with 
# maxium probability.If network is \code{Conditional Gaussian},
# a \code{NULL} value is returned.}
#############################################################################

absorb.gnbp=function(gpfit,node,evidence)

{
  
  requireNamespace("RHugin") || warning("Package not loaded: RHugin");
  
  
  ## get node attributes and network
  class_nodes=gpfit$gp_nodes
  network<-gpfit$gp
  type<-gpfit$gp_flag
  
  ## get d-connected nodes
  dnodes<-RHugin::get.dconnected.nodes(network,node)
  dnodes<-class_nodes[match(setdiff(class_nodes[match(dnodes,class_nodes),1],node),class_nodes),]
  dnodes<-matrix(dnodes,ncol=4)
  
  ## get marginal distribution
  marginal<-.get.marginal.bn(network,dnodes)
  
  
  
  ## check if class of evidence is matrix
  if(class(evidence)!="matrix")
    stop("In function(absorb.gpBP),'evidence' must be of class matrix")
  
  ## create matrices to store phenotype results
  if (type == "cg")  
  {
    JSI=matrix(nrow=length(dnodes[which(dnodes[,4]=="pheno" & dnodes[,2]=="numeric"),1]),ncol=dim(evidence)[2])
    belief_mean=matrix(nrow=length(dnodes[which(dnodes[,4]=="pheno"& dnodes[,2]=="numeric"),1]),ncol=dim(evidence)[2])
    belief_var=matrix(nrow=length(dnodes[which(dnodes[,4]=="pheno"& dnodes[,2]=="numeric"),1]),ncol=dim(evidence)[2]) 
  }
  
  if (type == "db") 
  {
    FC=matrix(nrow=length(dnodes[which(dnodes[,4]=="pheno"& dnodes[,2]=="factor"),1]),ncol=dim(evidence)[2])
    pheno_state=matrix(nrow=length(dnodes[which(dnodes[,4]=="pheno"& dnodes[,2]=="factor"),1]),ncol=dim(evidence)[2])
    belief_pheno_freq_list=list()
    belief_pheno_freq=vector()
  }
  
  ## create variables to store genotype results
    belief_geno_freq_list=list()
    belief_geno_freq=vector()
  
  ## index dconnected genotypes & phenotypes
  Y<-which(dnodes[,4]=="geno" & dnodes[,2]=="factor")
  X<-which(dnodes[,4]=="pheno" & dnodes[,2]=="numeric")
  Z<-which(dnodes[,4]=="pheno" & dnodes[,2]=="factor")

  ## Absorb evidence and calculate JSI/FC
  for (i in 1:dim(evidence)[2])
  {
   
    ##absorb evidence
    for (j in 1:length(node))
    {
      RHugin::set.finding(network,node[j],evidence[j,i])
    }
    RHugin::propagate(network)
    
    ## get belief
    belief<-.get.belief.bn(network,dnodes)

    if (type == "cg")
    {
    ## calculate JSI
    JSI[,i]<-.get.jsi(marginal,belief,dnodes) 
    
    ## extract beliefs
    belief_mean[,i]<-belief[X,1]
    belief_var[,i]<-belief[X,2]
    }
    
    if (type == "db")
    {
      ## calculate FC
      FC_temp<-.get.FC(marginal,belief,dnodes)
     
      ## extract FC
      FC[,i]<- FC_temp[,2]
      
      ## extract the state with maximum probability
      pheno_state[,i]<-FC_temp[,1]
      
      ## extract beliefs
      if(length(Z)!=0)
        belief_pheno_freq = cbind(belief_pheno_freq,matrix(belief[dnodes[Z,1],3:ncol(belief)],
                                                         nrow=length(dnodes[Z,1]),
                                                         ncol=as.numeric(max(dnodes[,3])),
                                                         dimnames=(list(dnodes[Z,1],
                                                                      colnames(belief)[3:ncol(belief)]))))
    }
        
      ## extract genotype beliefs
      if(length(Y)!=0)
           belief_geno_freq = cbind(belief_geno_freq,matrix(belief[dnodes[Y,1],3:ncol(belief)],
                                           nrow=length(dnodes[Y,1]),
                                           ncol=as.numeric(max(dnodes[,3])),
                                           dimnames=(list(dnodes[Y,1],
                                                          colnames(belief)[3:ncol(belief)]))))
    
    ##retract the evidence
    RHugin::retract(network)
  }
  
 
  ## annotate rows and columns
  if (type == "cg")
  {
    rownames(JSI)=dnodes[X,1]
    rownames(belief_mean)=dnodes[X,1]
    rownames(belief_var)=dnodes[X,1]
  }

  if (type == "db")
  {
    rownames(FC)=dnodes[Z,1]
    rownames(pheno_state)=dnodes[Z,1]
    
    FC=list(FC=FC,pheno_state=pheno_state)
    
    for (j in 1:as.numeric(max(dnodes[Z,3])))
    {
      belief_pheno_freq_temp = matrix(belief_pheno_freq[,seq(j,ncol(belief_pheno_freq),by=as.numeric(max(dnodes[Z,3])))],
                                     nrow=nrow(belief_pheno_freq),
                                     ncol=dim(evidence)[2],
                                     dimnames=list(rownames(belief_pheno_freq),NULL))
      
      name<-paste("state",j,sep="")
      belief_pheno_freq_list[[name]]= belief_pheno_freq_temp
    }
    
    phenomarginal<- matrix(marginal[dnodes[Z,1],3:ncol(marginal)],
                          nrow = length(dnodes[Z,1]),
                          ncol = as.numeric(max(dnodes[Z,3])),
                          dimnames = list(dnodes[Z,1],colnames(marginal)[3:ncol(marginal)]))
  }
    
  ## create a list for belief genotype frequencies
  if(length(belief_geno_freq)!=0)  
  {
    
    
    for (j in 1:as.numeric(max(dnodes[Y,3])))
    {
      belief_geno_freq_temp = matrix(belief_geno_freq[,seq(j,ncol(belief_geno_freq),by=as.numeric(max(dnodes[Y,3])))],
                                nrow=nrow(belief_geno_freq),
                                ncol=dim(evidence)[2],
                                dimnames=list(rownames(belief_geno_freq),NULL))
      
      name<-paste("state",j,sep="")
      belief_geno_freq_list[[name]]= belief_geno_freq_temp
    }
    
   genomarginal<- matrix(marginal[dnodes[Y,1],3:ncol(marginal)],
           nrow = length(dnodes[Y,1]),
           ncol = as.numeric(max(dnodes[Y,3])),
           dimnames = list(dnodes[Y,1],colnames(marginal)[3:ncol(marginal)]))
    
    
  }else
  {
    belief_geno_freq_list<-NULL
    genomarginal<-NULL
  }

 ## diagnostic plot

 if (dim(evidence)[2] >= 3)
 {
   if (type == "cg" & length(node)==1)
   {
   colpalette=rainbow(nrow(JSI))
   
   layout(matrix(c(1,2), 1, 2, byrow = F),c(2.5,1),c(1,1))
   
   for (i in 1:nrow(JSI))
   {
     plot(evidence,JSI[i,],type="l",xlab=NA,ylab=NA,lwd=2.5,col=colpalette[i],font.axis=10,cex=0.4,xlim=range(evidence),ylim=range(JSI))
     par(new=T)
   }
   title(xlab=paste(node,"evidence"), col.lab="black")
   title(ylab="JSI", col.lab="black")
   
   par(new=F)
   
   plot(c(1:5),type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
   
   legend("center",rownames(JSI),col=colpalette,bty="n",lwd=2.5,cex=0.8)
   }
   
 }


  ## store & return results
 if (type == "cg")
 {
  results=list(gp=network,
               gp_flag="cg",
               gp_nodes=class_nodes,
               evidence=evidence,
               node=node,
               marginal=list(pheno=list(mean = as.matrix(marginal[dnodes[X,1],1]),
                                        var = as.matrix(marginal[dnodes[X,1],2])),
                             geno=list(freq = genomarginal)),
               belief=list(pheno=list(mean=belief_mean,var=belief_var),geno=belief_geno_freq_list),
               JSI=JSI,
               FC=NULL)
 }
  if (type == "db")
  {
    results=list(gp=network,
                 gp_flag="db",
                 gp_nodes=class_nodes,
                 evidence=evidence,
                 node=node,
                 marginal=list(pheno = list(freq = phenomarginal),
                               geno = list(freq = genomarginal)),
                 belief=list(pheno = belief_pheno_freq_list,geno = belief_geno_freq_list),
                 JSI=NULL,
                 FC=FC)
  }

  class(results)<-"gnbp"
  
  return(results)
  
}



.get.marginal.bn=function(network,dnodes)
  
{
  
  ## create matrices to store results
  marg_mean<-matrix(nrow=dim(dnodes)[1],ncol=1)
  marg_var<-matrix(nrow=dim(dnodes)[1],ncol=1)
  marg_freq<-matrix(nrow=dim(dnodes)[1],ncol=as.numeric(max(dnodes[,3])))
  
  
    ## get mean and var
    for (j in 1:dim(dnodes)[1])
    {  
    
      pmarginal<-RHugin::get.marginal(network,dnodes[j,1])
      if(dnodes[j,2]=="numeric")
      {
        marg_mean[j,1]<-pmarginal$mean
        marg_var[j,1]<-unlist(pmarginal$cov)
      }
    
      if(dnodes[j,2]=="factor")
      {
        marg_freq[j,]<-pmarginal$table[,2]
      }
    
    }

  ## annotate rows and columns     
  
    if(length(marg_freq)!=0)
    {
      freqnames<-matrix(nrow=as.numeric(max(dnodes[,3])),ncol=1)
      
      for (i in 1:max(dnodes[,3]))
        freqnames[i]=paste("state",i,sep="")
      
      colnames(marg_freq)=freqnames
    }
  
   if(length(marg_mean)!=0)
   {
     colnames(marg_var)=c("var")
     colnames(marg_mean)=c("mean")
   }

marginal<-cbind(cbind(marg_mean,marg_var),marg_freq)
rownames(marginal)=dnodes[,1]
  
  ## return cpt
 
    return(marginal)
  
}



.get.belief.bn=function(network,dnodes)
{
  ## create matrices for storing mean and var of marginal distribution
  belief_mean<-matrix(nrow=dim(dnodes)[1],ncol=1)
  belief_var<-matrix(nrow=dim(dnodes)[1],ncol=1)
  belief_freq<-matrix(nrow=dim(dnodes)[1],ncol=as.numeric(max(dnodes[,3])))
  
  ## predictions for different genotypes
  ## Note:each column is replicated for coding purpose
  
  for (j in 1:(dim(dnodes)[1]))
  { 
    temp<-RHugin::get.belief(network,dnodes[j,1])
    
    if(dnodes[j,2]=="numeric")
    {
      belief_mean[j,1]<-temp[1]
      belief_var[j,1]<-temp[2]
    }
    
    if(dnodes[j,2]=="factor")
      belief_freq[j,]<-temp
    
  }
  
  
  ## annotate the rows and columns
  if(length(belief_freq)!=0)
  
  {
    freqnames<-matrix(nrow=as.numeric(max(dnodes[,3])),ncol=1)
  
    for (i in 1:max(dnodes[,3]))
      freqnames[i]=paste("state",i,sep="")

    colnames(belief_freq)=freqnames
  }

  if(length(belief_mean)!=0)
  {
    colnames(belief_mean)<-c("mean")
    colnames(belief_var)<-c("var")
  }
  
  belief<-cbind(cbind(belief_mean,belief_var),belief_freq)
  rownames(belief)<-dnodes[,1]
  
  return(belief)
  
}


.get.jsi=function(marginal,belief,dnodes)

  { 
  
    X = which(dnodes[,2]=="numeric" & dnodes[,4]=="pheno")
    KLDiv1 = matrix(nrow=length(dnodes[X,1]),ncol=1)
    KLDiv2 = matrix(nrow=length(dnodes[X,1]),ncol=1)
    jsi = matrix(nrow=length(dnodes[X,1]),ncol=1)
   
    
    for (i in 1:nrow(dnodes[X,]))
    {

        mu0=marginal[rownames(marginal)==dnodes[X[i],1],1]
        mu=belief[rownames(belief)==dnodes[X[i],1],1]
        
        sigma20=marginal[rownames(marginal)==dnodes[X[i],1],2]
        sigma2=belief[rownames(belief)==dnodes[X[i],1],2]
        
        KLDiv1[i,1]=0.5*(((mu-mu0)^2)/sigma2 +(sigma20/sigma2)-log(sigma20/sigma2)-1)
        KLDiv2[i,1]=0.5*(((mu0-mu)^2)/sigma20 +(sigma2/sigma20)-log(sigma2/sigma20)-1)
        jsi[i,1]=0.5*(KLDiv1[i,1]+KLDiv2[i,1])*sign(mu-mu0)     
    
      
    }
    
    return(jsi)
  
  }

.get.FC=function(marginal,belief,dnodes)
  
{ 
  Z<- which(dnodes[,2]=="factor" & dnodes[,4]=="pheno")
  FC = matrix(nrow=length(dnodes[Z,1]),ncol=2)
  
  for (i in 1:nrow(dnodes[Z,]))
  {
      FC[i,1] = which.max(belief[dnodes[Z[i],1],])
      FC[i,2] = belief[dnodes[Z[i],1],FC[i,1]]/marginal[dnodes[Z[i],1],FC[i,1]]  
      FC[i,1] = FC[i,1]-2
  }
  
  return(FC)
  
}
