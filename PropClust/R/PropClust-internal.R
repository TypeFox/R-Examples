#

.is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


.propensityClustering.internal <-function(adjacency, initialClusters, l2bool, nClusters, initbool,
                                fastUpdates, accelerated = TRUE){
        #if(l2bool>0){
        #  print("Using L2 updates")
        #}else{
#                print("Using Poisson updates")
#        }
#        for(i in 1:length(initialClusters)){
#                if(initialClusters[i]<1){
#                        stop("Initial Cluster values must be integers greater than 0.")
#                }
#        }
#        if(!.is.wholenumber(initialClusters[1])){
#                print("Converting Cluster values to integers and storing in initialClusters.norm")
#                initialClusters.norm = as.numeric(as.factor(initialClusters))
#        }else{
#                initialClusters.norm=initialClusters
#        }
        nodes=length(initialClusters)

#        if(nodes!=length(adjacency[1,])|nodes!=length(adjacency[,1])){
#                stop("Adjacency must have same length as initial clustering")
#        }
#        if(nClusters<2){
#                stop("Number of clusters must be > 1. Use propensityDecomposition for 1 cluster.")
#        }
        Phat=rep(0.7, nodes)
        Ahat=matrix(0.5, nClusters,nClusters)
        diag(Ahat) = 1;
        fact=0.
        lnorm=0.
        if (accelerated)
        {
          compiledFnc = if(fastUpdates) "propclusttrial" else "propclustaccel";
        } else 
          compiledFnc = "propensityclustering"
   
        results = .Fortran(compiledFnc, ADJ=as.single(adjacency),
                         Clustering=as.integer(initialClusters),Propensity=as.double(Phat),
                         IntermodularAdjacency=as.double(Ahat),Factorizability=as.double(fact),
                         Criteria=as.double(lnorm),Nodes=as.integer(nodes),Clusters=as.integer(nClusters),
                         L2=as.integer(l2bool),Init=as.integer(initbool))

        resultsmod=results[2:6]
        if(l2bool>0){
          resultsmod$L2Norm=results$Criteria
        }else{
          resultsmod$Loglik=results$Criteria
        }
        dim ( resultsmod$IntermodularAdjacency ) = c(nClusters, nClusters)
#        if(!.is.wholenumber(initialClusters[1])){
#                finalClusters = .translateUsingTable(results$Clustering, 
#                                    .translationTable(initialClusters.norm, initialClusters))
#        }else{
#                finalClusters = results$Clustering
#        }
#        
#        resultsmod$Clustering=finalClusters
#        
        meanMat = matrix((results$ADJ),ncol=nodes);
        resultsmod$MeanValues=as.dist(meanMat);
        resultsmod$TailPvalues = as.dist(t(meanMat));
        return(resultsmod)
}


CPBADecomposition<-function(adjacency,
                            clustering,
                            nClusters = NULL,
                            objectiveFunction = c("Poisson", "L2norm"),
                            dropUnassigned = TRUE,
                            unassignedLabel = 0,
                            unassignedMethod = "average",
                            accelerated = TRUE,
                            parallel = FALSE)
{
  .checkAdjMat(adjacency, min = 0, max = max(adjacency, na.rm = TRUE));
  objectiveFunction = match.arg(objectiveFunction);
  nAllNodes = nNodes = ncol(adjacency);
  useNodes = rep(TRUE, nNodes);

  if (parallel) 
  {
    warning("Parallel version does not work yet... using standard accelerated calculations.");
    parallel = FALSE;
  }

  if (length(clustering)!=nNodes) 
    stop("Length of 'clustering' must be the same as the number of nodes (columns) in 'adjacency'.");

  if (is.null(nClusters))
  {
    if (any(is.na(clustering))) stop("All entries in 'clustering' must be present (non-NA)");
    if (all(clustering==unassignedLabel)) 
       stop("All entries in 'clustering' are unassigned.");
    if (dropUnassigned)
    {
      useNodes = clustering != unassignedLabel;
      adjacency = adjacency[useNodes, useNodes];
      clustering = clustering[useNodes];
      nNodes = ncol(adjacency);
    } else {
      clustering = .assignToNearestCluster(1-adjacency, labels = clustering,
                                  method = unassignedMethod, unassignedLabel = unassignedLabel);
    }
  } else {
    if (nClusters != 1) stop("If given, number of clusters must be 1.");
    clustering = rep(1, nNodes);
  }

  nClusters = length(unique(clustering));

  clustering.norm = as.numeric(as.factor(clustering))
  propensity.all = rep(NA, nAllNodes);

  Phat=rep(0.7, nNodes)
  fact=0
  lnorm=0

  if(nClusters>1)
  {
    Ahat=matrix(0.5, nClusters,nClusters)
    diag(Ahat) = 1;
    compiledFnc =  if (accelerated) "propdecompaccel" else "propensitydecomposition";
    results=.Fortran(compiledFnc, 
                     ADJ=as.single(adjacency),
                     Clustering=as.integer(clustering.norm),
                     Propensity=as.double(Phat),
                     IntermodularAdjacency=as.double(Ahat),
                     Factorizability=as.double(fact),
                     Criteria=as.double(lnorm),Nodes=as.integer(nNodes),
                     Clusters=as.integer(nClusters),L2=as.integer(objectiveFunction=="L2norm"))

    resultsmod=results[3:5]
    resultsmod$IntermodularAdjacency=matrix((resultsmod$IntermodularAdj),ncol=nClusters)
    finalClusters = .translateUsingTable(results$Clustering, .translationTable(clustering.norm, clustering))
  } else {
    results=.Fortran("singleclusterupdate",ADJ=as.single(adjacency),
                     Propensity=as.double(Phat),Factorizability=as.double(fact),
                     Criteria=as.double(lnorm),Nodes=as.integer(nNodes),
                     L2=as.integer(objectiveFunction=="L2norm"))
    resultsmod=results[2:3]
  }
  propensity.all[useNodes] = resultsmod$Propensity
  resultsmod$Propensity = propensity.all

  if (objectiveFunction=="L2norm") 
  {
     resultsmod$L2Norm=results$Criteria
  } else {
     resultsmod$Loglik=results$Criteria
  }
  meanMat = matrix(NA, nAllNodes, nAllNodes)
  meanMat[useNodes, useNodes] = matrix((results$ADJ),ncol=nNodes);
  resultsmod$ExpectedAdjacency=as.dist(meanMat);
  resultsmod$EdgePvalues = as.dist(t(meanMat));
  return(resultsmod)
}

#.clusterWordPairs<-function(tfile,clusters){
#        words=0
#results1=.Fortran("countwords",textfile=as.character(tfile),Words=as.integer(words))
#nodes=results1$Words
#Phat=mat.or.vec(nodes,1)
#Ahat=mat.or.vec(clusters,clusters)
#Phat[]=.7
#Ahat[,]=.5
#for(i in 1:length(Ahat[1,])){
#        Ahat[i,i]=1
#}
#testmodule=mat.or.vec(nodes,1)
#testmodule[]=1
#ADJ=mat.or.vec(nodes,nodes)
#ADJ[,]=0
#        loglik=0
#results=.Fortran("wordpairclusters",textfile=as.character(tfile),ADJ=as.integer(ADJ),Testmodule=as.integer(testmodule),Propensity=as.double(Phat),IntermodularAdjacency=as.double(Ahat),Loglik=as.double(loglik),Words=as.integer(nodes),Clusters=as.integer(clusters))
##results$wordlist=scan(file="propClustTemp.txt",what=character(),skip=4)
#resultsmod=results[3:7]
#resultsmod$Wordlist=scan(file="propClustTemp.txt",what=character(),skip=4)
#resultsmod$Propensity=results$Propensity[1:results$Words]
#resultsmod$Testmodule=results$Testmodule[1:results$Words]
#resultsmod$IntermodularAdjacency=matrix((results$IntermodularAdjacency),ncol=clusters)
#return(resultsmod)
#}


#.clusterDisorderPairs<-function(tfile,clusters){
#        words=0
#results1=.Fortran("countgenes",textfile=as.character(tfile),Words=as.integer(words))
#nodes=results1$Words
#nodes2=nodes
#results2=.Fortran("countdisorders",textfile=as.character(tfile),Dis=as.integer(words))
#Dis=results2$Dis
#Phat=mat.or.vec(nodes,1)
#Ahat=mat.or.vec(clusters,clusters)
#Phat[]=.7
#Ahat[,]=.5
#for(i in 1:length(Ahat[1,])){
#        Ahat[i,i]=1
#}
#testmodule=mat.or.vec(nodes,1)
#testmodule[]=1
#ADJ=mat.or.vec(nodes,nodes)
#ADJ[,]=0
#        loglik=0
#results=.Fortran("omimmorbidmap",textfile=as.character(tfile),ADJ=as.single(ADJ),Testmodule=as.integer(testmodule),Propensity=as.double(Phat),IntermodularAdjacency=as.double(Ahat),Loglik=as.double(loglik),Genes=as.integer(nodes),Disorders=as.integer(Dis),Genes2=as.integer(nodes2),Clusters=as.integer(clusters))
#resultsmod=results[2:8]
#resultsmod$IntermodularAdjacency=matrix((results$IntermodularAdjacency),ncol=clusters)
#resultsmod$ADJ=matrix((results$ADJ),ncol=nodes)
#resultsmod$ADJ=resultsmod$ADJ[1:results$Genes2,1:results$Genes2]
#res=scan(file="propClustTempOrderedDisorders.txt",what=list(character(),numeric(),integer()),skip=8+clusters,sep="|")
#resultsmod$Genelist=res[1]
#return(resultsmod)
#}


#.clusterGenePairs<-function(tfile,clusters){
#        words=0
#results1=.Fortran("countgenes",textfile=as.character(tfile),Words=as.integer(words))
#nodes=results1$Words
#nodes2=nodes
#results2=.Fortran("countdisorders",textfile=as.character(tfile),Dis=as.integer(words))
#Dis=results2$Dis
#Phat=mat.or.vec(nodes,1)
#Ahat=mat.or.vec(clusters,clusters)
#Phat[]=.7
#Ahat[,]=.5
#for(i in 1:length(Ahat[1,])){
#        Ahat[i,i]=1
#}
#testmodule=mat.or.vec(nodes,1)
#testmodule[]=1
#ADJ=mat.or.vec(nodes,nodes)
#ADJ[,]=0
#        loglik=0
#results=.Fortran("omimgenemap",textfile=as.character(tfile),ADJ=as.single(ADJ),Testmodule=as.integer(testmodule),Propensity=as.double(Phat),IntermodularAdjacency=as.double(Ahat),Loglik=as.double(loglik),Genes=as.integer(nodes),Disorders=as.integer(Dis),Genes2=as.integer(nodes2),Clusters=as.integer(clusters))
#resultsmod=results[2:8]
#resultsmod$IntermodularAdjacency=matrix((results$IntermodularAdjacency),ncol=clusters)
#resultsmod$ADJ=matrix((results$ADJ),ncol=nodes)
#resultsmod$ADJ=resultsmod$ADJ[1:results$Genes2,1:results$Genes2]
#res=scan(file="propClustTempOrderedGenes.txt",what=list(character(),numeric(),integer()),skip=8+clusters,sep="|")
#resultsmod$Genelist=res[1]
#return(resultsmod)
#}


#.clusterCompanyPairs<-function(tfile,clusters){
#        mems=0
#results2=.Fortran("countitems",textfile=as.character(tfile),Members=as.integer(mems),Position=as.integer(2))
#Members=results2$Members
#        comps=0
#results1=.Fortran("countitems",textfile=as.character(tfile),Comps=as.integer(comps),Position=as.integer(3))
#nodes=results1$Comps
#nodes2=nodes
#Phat=mat.or.vec(500,1)
#Ahat=mat.or.vec(clusters,clusters)
#Phat[]=.7
#Ahat[,]=.5
#for(i in 1:length(Ahat[1,])){
#        Ahat[i,i]=1
#}
#        testmodule=mat.or.vec(500,1)
#        testmodule[]=1
#        ADJ=mat.or.vec(500,500)
#        ADJ[,]=0
#        loglik=0
#results=.Fortran("clustercompanies",textfile=as.character(tfile),ADJ=as.single(ADJ),Testmodule=as.integer(testmodule),Propensity=as.double(Phat),IntermodularAdjacency=as.double(Ahat),Loglik=as.double(loglik),Comps=as.integer(nodes),Members=as.integer(Members),Comps2=as.integer(nodes2),Clusters=as.integer(clusters))
#resultsmod=results[2:8]
#resultsmod$IntermodularAdjacency=matrix((results$IntermodularAdjacency),ncol=clusters)
#resultsmod$ADJ=matrix((results$ADJ),ncol=500)
#resultsmod$ADJ=resultsmod$ADJ[1:results$Comps2,1:results$Comps2]
#res=scan(file="propClustTempOrderedComps.txt",what=list(character(),numeric(),integer()),skip=8+clusters,sep="|")
#resultsmod$Genelist=res[1]
#return(resultsmod)
#}



.speed1<-function(adjacency,nodes){
        results=.Fortran("speedtest1",ADJ=as.double(adjacency),Nodes=as.integer(nodes))
        resultsmod=results$Nodes
        return(resultsmod)
}

.speed2<-function(adjacency,nodes){
        results=.Fortran("speedtest2",ADJ=as.double(adjacency),Nodes=as.integer(nodes))
        resultsmod=results$Nodes
        return(resultsmod)
}
