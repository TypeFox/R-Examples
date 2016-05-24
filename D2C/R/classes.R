#' @import RBGL gRbase randomForest Rgraphviz methods foreach 


#########################################
########   class D2C.descriptor
#########################################

##' An S4 class to store the descriptor parameters
setClass("D2C.descriptor",
         slots = list(lin="logical", acc="logical",
                      struct="logical",pq="numeric",
                      bivariate="logical",ns="numeric"))

##' creation of a D2C.descriptor 
##' @param .Object : the D2C.descriptor object
##' @param lin :	TRUE OR FALSE: if TRUE it uses a linear model to assess a dependency, otherwise a local learning algorithm
##' @param acc : TRUE OR FALSE: if TRUE it uses the accuracy of the regression as a descriptor
##'  @param struct	: TRUE or FALSE to use the ranking in the markov blanket as a descriptor
##'  @param pq :a vector of quantiles used to compute the descriptors
##'  @param bivariate :TRUE OR FALSE: if TRUE it includes also the descriptors of the bivariate dependence
##'  @param ns : size of the Markov Blanket returned by the mIMR algorithm
##' @references Gianluca Bontempi, Maxime Flauder (2014) From dependency to causality: a machine learning approach. Under submission
##' @examples
##' require(RBGL)
##' require(gRbase)
##' require(foreach)
##'descr.example<-new("D2C.descriptor",bivariate=FALSE,ns=3,acc=TRUE)
##'trainDAG<-new("simulatedDAG",NDAG=2, N=50,noNodes=10,
##'              functionType = "linear", seed=0,sdn=0.5) 
##'
##' @export
setMethod("initialize",
          "D2C.descriptor",
          function(.Object, lin=TRUE, acc=TRUE,
                   struct=TRUE,pq=c(0.1, 0.25, 0.5, 0.75, 0.9),
                   bivariate=FALSE,ns=4)
          {
            
            .Object@lin <- lin
            .Object@acc <- acc
            .Object@struct <- struct
            .Object@bivariate <- bivariate
            .Object@pq <- pq
            .Object@ns <- ns
            .Object
          }
)




#########################################
########   class DAG.network
#########################################




##' An S4 class to store DAG.network
##' @param network : object of class "graph"
setClass("DAG.network",  slots = list(network = "graph"))



##' creation of a DAG.network
##' @param .Object : DAG.network object
##' @param sdn : standard deviation of aditive noise. 
##' @param sigma : function returning the additive noise
##' @param H : function describing the type of the dependency. 
##' @param network : object of class "graph"
##' @export
setMethod("initialize", signature="DAG.network",
          function(.Object, network,
                   sdn=0.5,
                   sigma=function(x) return(rnorm(n = 1,sd = sdn)),
                   H=function(x) return(H_Rn(1)) )
          {
            DAG = network
            if(!is.DAG(DAG))
            {
              stop("it is not a DAG")
            } else  {
              nodeDataDefaults(DAG,"bias") <-0
              nodeDataDefaults(DAG,"sigma") <-sigma
              edgeDataDefaults(DAG,"H") <- function(x) return(x)
              for( edge in edgeList(DAG)){
                edgeData(DAG, from=edge[1], to=edge[2], attr="weight") <- runif(1,0.5,1)*sample(c(-1,1),1)
                edgeData(DAG, from=edge[1], to=edge[2],attr="H") <- H()
                
              }     
            }   	 
            .Object@network <- DAG
            
            return(.Object)
          }
)




#' @docType methods
setGeneric("compute", function(object,...) {standardGeneric("compute")})


##' compute N samples according to the network distribution
##' @param N  numeric. the number of samples generated according to the network
##' @param object a DAG.network object
##' @return a N*nNodes matrix
##' @export
setMethod("compute", signature="DAG.network", function(object,N=50)
{
  if(!is.numeric(N))
    stop("N is not numeric")
  
  DAG = object@network
  nNodes <- numNodes(DAG)
  topologicalOrder <-tsort(DAG)
  
  D <- matrix(NA,nrow=N,ncol=nNodes)
  colnames(D) <- topologicalOrder  
  
  
  for (i in topologicalOrder){   
    bias = nodeData(DAG,n=i,attr="bias")[[1]] 
    sigma = nodeData(DAG,n=i,attr="sigma")[[1]] 
    inEdg <-  inEdges(node=i,object=DAG)[[1]]
    
    if (length(inEdg)==0){
      D[,i]<-bias + replicate(N,sigma())
    } else  {
      D[,i]<-bias
      for(j in  inEdg)
        ## it computes the linear combination of the inputs
      {
        inputWeight = edgeData(self=DAG,from=j,to=i,attr="weight")[[1]]
        H = edgeData(self=DAG,from=j,to=i,attr="H")[[1]]
        D[,i]<- D[,i] + H(D[,j]) *  inputWeight
        D[,i] <- scale(D[,i]) + replicate(N,sigma())
      }
    }    
  }
  col.numeric<-as(colnames(D),"numeric")
  D<-D[,topologicalOrder[order(col.numeric)]]
  
  return(D)
})


#########################################
########   class simulatedDAG
#########################################

##' An S4 class to store a list of DAGs and associated observations
##' @param list.DAGs : list of stored DAGs
##' @param list.observationsDAGs : list of observed datasets, each sampled from the corresponding member of list.DAGs 
##' @param NDAG  : number of DAGs. 
##' @param functionType : type of the dependency. It is of class "character" and is one of  ("linear", "quadratic","sigmoid")
##' @param seed : random seed
setClass("simulatedDAG",
         slots = list(list.DAGs="list",list.observationsDAGs="list",
                      NDAG="numeric", functionType="character",seed="numeric"))




##' creation of a "simulatedDAG" containing a list of DAGs and associated observations 
##' @param .Object : simulatedDAG object
##' @param NDAG : number of DAGs to be created and simulated
#' @param noNodes  : number of Nodes of the DAGs. If it is a two-valued vector , the value of Nodes is randomly sampled in the interval 
#' @param N  : number of sampled observations for each DAG. If it is a two-valued vector [a,b], the value of N is randomly sampled in the interval [a,b]
#' @param sdn : standard deviation of aditive noise. If it is a two-valued vector, the value of N is randomly sampled in the interval
#' @param seed : random seed
#' @param verbose : if TRUE it prints out the state of progress
#' @param functionType : type of the dependency. It is of class "character" and is one of  ("linear", "quadratic","sigmoid")
#'  @param quantize  : if TRUE it discretize the observations into two bins. If it is a two-valued vector [a,b], the value of quantize is randomly sampled in the interval [a,b]
#'  @param goParallel : if TRUE it uses parallelism
#' @references Gianluca Bontempi, Maxime Flauder (2014) From dependency to causality: a machine learning approach. Under submission
#' @examples
#' require(RBGL)
#' require(gRbase)
#' require(foreach)
#' descr=new("D2C.descriptor")
#'descr.example<-new("D2C.descriptor",bivariate=FALSE,ns=3,acc=TRUE)
#'trainDAG<-new("simulatedDAG",NDAG=10, N=c(50,100),noNodes=c(15,40),
#'              functionType = "linear", seed=0,sdn=c(0.45,0.75)) 
#' @export
#' 
#' 
setMethod("initialize",
          "simulatedDAG",
          function(.Object, NDAG=1,
                   noNodes=sample(10:20,size=1),functionType="linear",
                   quantize=FALSE,
                   verbose=TRUE,N=sample(100:500,size=1),
                   seed=1234,sdn=0.5, goParallel=FALSE)
          {
            
            ##generate a training set
            ## NDAG the number of network to use
            ##functionType example : "R1" "R2" "sigmoid1"
            
            
            
            `%op%` <- if (goParallel) `%dopar%` else `%do%`
            
            .Object@functionType=functionType
            .Object@seed=seed
            X=NULL
            Y=NULL
            list.DAGs=NULL
            list.observationsDAGs=NULL
            if (NDAG<=0)
              return(.Object)
            
           
            
            FF<-foreach (i=1:NDAG) %op%{         
        ##  for (i in 1:NDAG){  
         set.seed(seed+i)                
              N.i<-N
              if (length(N)>1)
                N.i<-sample(N[1]:N[2],1)
              
              quantize.i<-quantize
              if (length(quantize)>1)
                quantize.i<-sample(quantize,1)
              
              noNodes.i<-noNodes
              if (length(noNodes)>1)
                noNodes.i<-sample(noNodes[1]:noNodes[2],1)
              
              sdn.i<-sdn
              if (length(sdn)>1)
                sdn.i<-runif(1,sdn[1],sdn[2])
              
              functionType.i<-functionType
              if (length(functionType.i)>1)
                functionType.i<-sample(functionType,1)
              
              
              
              
              V=1:noNodes.i
              
              maxpar = sample(2:round(noNodes.i/2),size=1)
              
              
              if(functionType.i=="linear"){
                H = function() return(H_Rn(1))
                
              }else if(functionType.i=="quadratic"){
                H = function() return(H_Rn(2))
                
              }else if(functionType.i=="sigmoid"){
                H = function() return(H_sigmoid(1))
                
              }
              wgt = runif(n = 1,min = 0.65,max = 0.85)
              netwDAG<-random_dag(V,maxpar = maxpar,wgt)
              cnt<-2
              while (sum(unlist(lapply(edges(netwDAG),length)))<3 & cnt<100){
                netwDAG<-random_dag(V,maxpar = maxpar,wgt)
                wgt = runif(n = 1,min = 0.65/cnt,max = 0.85)
                cnt<-cnt+1
                
              }
              
              
              DAG = new("DAG.network",
                        network=netwDAG,H=H,sdn=sdn.i)
              
              
              observationsDAG = compute(DAG,N=N.i)
              if (quantize.i)
                observationsDAG<-apply(observationsDAG,2,quantization)
              
              if (verbose)
                cat("simulatedDAG: DAG number:",i,"generated: #nodes=", length(V), 
                    "# edges=",sum(unlist(lapply(edges(netwDAG),length))), "# samples=", N.i, "\n")
              
              
              
              
              list(observationsDAG=observationsDAG,netwDAG=netwDAG)
            } ## foreach
            
            
            .Object@list.DAGs=lapply(FF,"[[",2)
            .Object@list.observationsDAGs=lapply(FF,"[[",1)
            to.remove=which(unlist(lapply(lapply(.Object@list.DAGs,edgeList),length))==0)
            if (length(to.remove)>0){
              .Object@list.DAGs=.Object@list.DAGs[-to.remove]
              .Object@list.observationsDAGs=.Object@list.observationsDAGs[-to.remove]
            }
            .Object@NDAG=length(.Object@list.DAGs)
            
            .Object
          }
)




setGeneric("update", def=function(object,...) {standardGeneric("update")})

#' update of a "simulatedDAG" with a list of DAGs and associated observations 
#' @param object :  simulatedDAG to be updated
#' @param list.DAGs : list of stored DAGs
#' @param list.observationsDAGs : list of observed datasets, each sampled from the corresponding member of list.DAGs 
#' @export
setMethod(f="update",
          signature="simulatedDAG",
          definition=function(object,list.DAGs,list.observationsDAGs) {
            if (length(list.DAGs)!=length(list.observationsDAGs))
              stop("Lists with different lengths !")
            object@list.DAGs=c(object@list.DAGs,list.DAGs)
            object@list.observationsDAGs=c(object@list.observationsDAGs,list.observationsDAGs)
            object@NDAG=length(object@list.DAGs)
            object
            
          }
)


#########################################
########   class D2C
#########################################

setOldClass("randomForest")

#' An S4 class to store the RF model trained on the basis of the descriptors of NDAG DAGs
setClass("D2C",
         slots = list(mod="randomForest", X="matrix",Y="numeric",                    
                      descr="D2C.descriptor",features="numeric",rank="numeric",
                      allEdges="list",ratioMissingNode="numeric",ratioEdges="numeric",
                      max.features="numeric"
         ))

#' creation of a D2C object which preprocesses the list of DAGs and observations contained in sDAG and fits a  Random Forest classifier
#' @param .Object : the D2C object
#' @param sDAG : simulateDAG object
#' @param descr  : D2C.descriptor object containing the parameters of the descriptor
#' @param max.features  : maximum number of features used by the Random Forest classifier \link[randomForest]{randomForest}. The features are selected by the importance returned by the function \link[randomForest]{importance}.
#' @param ratioEdges  : percentage of existing edges which are added to the training set
#' @param ratioMissingNode  : percentage of existing nodes which are not considered. This is used to emulate latent variables.
#' @param goParallel : if TRUE it uses parallelism
#' @param verbose  : if TRUE it prints the state of progress
#' @references Gianluca Bontempi, Maxime Flauder (2014) From dependency to causality: a machine learning approach. Under submission
#' @examples
#' require(RBGL)
#' require(gRbase)
#'  require(foreach)
#' descr=new("D2C.descriptor")
#'descr.example<-new("D2C.descriptor",bivariate=FALSE,ns=3,acc=TRUE)
#'trainDAG<-new("simulatedDAG",NDAG=2, N=50,noNodes=10,
#'              functionType = "linear", seed=0,sdn=0.5) 
#' example<-new("D2C",sDAG=trainDAG, descr=descr.example)
#' @export
setMethod("initialize",
          "D2C",
          function(.Object, sDAG, 
                   descr=new("D2C.descriptor"),
                   verbose=TRUE, 
                   ratioMissingNode=0,
                   ratioEdges=1,max.features=20, goParallel=FALSE ) {
            
            #generate a training set
            # NDAG the number of network to use
            #functionType example : "R1" "R2" "sigmoid1"
            `%op%` <- if (goParallel) `%dopar%` else `%do%`
            
            .Object@descr=descr
            .Object@ratioMissingNode= ratioMissingNode
            .Object@ratioEdges= ratioEdges
            .Object@max.features= max.features
            X=NULL
            Y=NULL
            allEdges=NULL
            FF<-NULL
            
            FF<-foreach (i=1:sDAG@NDAG) %op%{
              set.seed(i)
              DAG = sDAG@list.DAGs[[i]]
              observationsDAG =sDAG@list.observationsDAGs[[i]]
              
              Nodes = nodes(DAG)
              
              sz=max(2,ceiling(length(Nodes)*(1-ratioMissingNode)))
              keepNode = sort(sample(Nodes,
                                     size = sz ,
                                     replace = F))  
              
              DAG2 =subGraph(keepNode, DAG)
              
              
              ##choose wich edge to train / predict and find the right label
              nEdge = length(edgeList(DAG))
              sz=max(1,round(nEdge*ratioEdges))
              
              edges = matrix(unlist(sample(edgeList(DAG2),
                                           size = sz,replace = F)),ncol=2,byrow = TRUE)  
              edges = rbind(edges,t(replicate(n =sz ,
                                              sample(keepNode,size=2,replace = FALSE))))
              
              nEdges =  NROW(edges)
              labelEdge = numeric(nEdges)
              for(j in 1:nEdges){
                I =edges[j,1] ; 
                J =edges[j,2] ;
                labelEdge[j] = as.numeric(I %in% inEdges(node = J,DAG2)[[1]])
              } 
              
              
              ##compute the descriptor for the edges 
              nNodes = length(labelEdge)
              
              X.out = NULL
              for(j in 1:nNodes){
                I =as(edges[j,1],"numeric") ; 
                J =as(edges[j,2],"numeric") ; 
                
                
                d<-descriptor(observationsDAG,I,J,lin=descr@lin,acc=descr@acc,
                              struct=descr@struct,bivariate=descr@bivariate,
                              pq=descr@pq,ns=descr@ns)
                
                X.out = rbind(X.out,d)
              }
              if (verbose)
                cat("D2C:  DAG", i, " processed \n")
              
              list(X=X.out,Y=labelEdge,edges=edges)
              
            }
            
            X<-do.call(rbind,lapply(FF,"[[",1))
            Y<-do.call(c,lapply(FF,"[[",2))
            allEdges<-lapply(FF,"[[",3)
            
            
            features<-1:NCOL(X)
            wna<-which(apply(X,2,sd)<0.01)
            if (length(wna)>0)
              features<-setdiff(features,wna)
            
            X<-scale(X[,features])
            .Object@features=features
            .Object@X=X
            .Object@Y=Y
            .Object@allEdges=allEdges
            RF <- randomForest(x =X ,y = factor(Y),importance=TRUE)
            IM<-importance(RF)[,"MeanDecreaseAccuracy"]
            rank<-sort(IM,decr=TRUE,ind=TRUE)$ix[1:min(max.features,NCOL(X))]
            RF <- randomForest(x =X[,rank] ,y = factor(Y))
            .Object@rank=rank
            .Object@mod=RF
            
            .Object
          }
)



#' @docType methods
setGeneric("updateD2C", def=function(object,...) {standardGeneric("updateD2C")})

#' update of a "D2C" with a list of DAGs and associated observations 
#' @param object :  D2C to be updated
#' @param sDAG : simulatedDAG object to update D2C
#' @param verbose : TRUE or FALSE
#' @param goParallel : if TRUE it uses  parallelism  
#' @export
setMethod(f="updateD2C",
          signature="D2C",
          definition=function(object,sDAG,           
                              verbose=TRUE, goParallel= FALSE){
            
            `%op%` <- if (goParallel) `%dopar%` else `%do%`
            ratioMissingNode=object@ratioMissingNode
            ratioEdges=object@ratioEdges
            descr=object@descr
            FF<-foreach (i=1:sDAG@NDAG) %op%{
              
              set.seed(i)
              DAG = sDAG@list.DAGs[[i]]
              observationsDAG =sDAG@list.observationsDAGs[[i]]
              
              Nodes = nodes(DAG)
              
              sz=max(2,ceiling(length(Nodes)*(1-ratioMissingNode)))
              keepNode = sort(sample(Nodes,
                                     size = sz ,
                                     replace = F))  
              
              DAG2 =subGraph(keepNode, DAG)
              
              
              ##choose wich edge to train / predict and find the right label
              nEdge = length(edgeList(DAG))
              sz=max(1,round(nEdge*ratioEdges))
              
              edges = matrix(unlist(sample(edgeList(DAG2),
                                           size = sz,replace = F)),ncol=2,byrow = TRUE)  
              edges = rbind(edges,t(replicate(n =sz ,
                                              sample(keepNode,size=2,replace = FALSE))))
              
              nEdges =  NROW(edges)
              labelEdge = numeric(nEdges)
              for(j in 1:nEdges){
                I =edges[j,1] ; 
                J =edges[j,2] ;
                labelEdge[j] = as.numeric(I %in% inEdges(node = J,DAG2)[[1]])
              } 
              
              
              ##compute the descriptor for the edges 
              nNodes = length(labelEdge)
              
              X.out = NULL
              for(j in 1:nNodes){
                I =as(edges[j,1],"numeric") ; 
                J =as(edges[j,2],"numeric") ; 
                
                
                d<-descriptor(observationsDAG,I,J,lin=descr@lin,acc=descr@acc,
                              struct=descr@struct,bivariate=descr@bivariate,
                              pq=descr@pq,ns=descr@ns)
                
                X.out = rbind(X.out,d)
              }
              if (verbose)
                cat("D2C:  DAG", i, " processed \n")
              
              list(X=X.out,Y=labelEdge,edges=edges)
              
            }
            
            X<-do.call(rbind,lapply(FF,"[[",1))
            Y<-do.call(c,lapply(FF,"[[",2))
            allEdges<-lapply(FF,"[[",3)
            
            
            
            X<-scale(X[,object@features],attr(object@X,"scaled:center"),attr(object@X,"scaled:scale"))
            
            object@X=rbind(object@X,X)
            object@Y=c(object@Y,Y)
            object@allEdges=c(object@allEdges,allEdges)
            RF <- randomForest(x =object@X ,y = factor(object@Y),importance=TRUE)
            IM<-importance(RF)[,"MeanDecreaseAccuracy"]
            rank<-sort(IM,decr=TRUE,ind=TRUE)$ix[1:min(object@max.features,NCOL(X))]
            RF <- randomForest(x =object@X[,rank] ,y = factor(object@Y))
            object@rank=rank
            object@mod=RF
            
            object  
            
          }
)

#' predict if there is a connection between node i and node j 
#' @param object : a D2C object 
#' @param i :  index of putative cause (\eqn{1 \le i \le n})
#' @param j  : index of putative effect (\eqn{1 \le j \le n})
#' @param data : dataset of observations from the DAG    
#' @return list with  response and prob of the prediction
#' @examples  
#' require(RBGL)
#' require(gRbase)
#' require(foreach)
#' data(example)
#'## load the D2C object
#' testDAG<-new("simulatedDAG",NDAG=1, N=50,noNodes=5,
#'            functionType = "linear", seed=1,sdn=c(0.25,0.5))
#' ## creates a simulatedDAG object for testing
#' plot(testDAG@@list.DAGs[[1]])
#' ## plot the topology of the simulatedDAG 
#' predict(example,1,2, testDAG@@list.observationsDAGs[[1]])
#' ## predict if the edge 1->2 exists
#' predict(example,4,3, testDAG@@list.observationsDAGs[[1]])
#' ## predict if the edge 4->3 exists
#' predict(example,4,1, testDAG@@list.observationsDAGs[[1]])
#' ## predict if the edge 4->1 exists
#' @references Gianluca Bontempi, Maxime Flauder (2014) From dependency to causality: a machine learning approach. Under submission
#' @export
setMethod("predict", signature="D2C",
          function(object,i,j,data)
          {
            out = list() 
            
            if (any(apply(data,2,sd)<0.01))
              stop("Error in D2C::predict: Remove constant variables from dataset. ")
            X_descriptor = descriptor(data,i,j,lin = object@descr@lin,
                                      acc = object@descr@acc,ns=object@descr@ns,
                                      struct = object@descr@struct,
                                      pq = object@descr@pq, bivariate =object@descr@bivariate)
            
            if (any(is.infinite(X_descriptor)))
              stop("Error in D2C::predict: infinite value ")
            
            X_descriptor=X_descriptor[object@features]
            
            
            X_descriptor=scale(array(X_descriptor,c(1,length(X_descriptor))),
                               attr(object@X,"scaled:center"),attr(object@X,"scaled:scale"))
            if (any(is.infinite(X_descriptor[,object@rank])))
              stop("error in D2C::predict")
            out[["response"]] = predict(object@mod, X_descriptor[,object@rank], type="response")
            out[["prob"]] = predict(object@mod, X_descriptor[,object@rank], type="prob")
            
            
            return(out)  
          })




#' Dataset example
#'@title stored D2C object
#'@description small D2C object for testing D2C functionalities
#' @name example
#' @docType data
#' @keywords data
#' @export
#' @examples
#' require(RBGL)
#' require(gRbase)
#' data(example)
#' print(example@@mod)
#' ## Random Forest
#' print(dim(example@@X))
#' ## dimension of the training set
NULL




#' Dataset alarm
#'@title Alarm dataset
#'@description contains the adjacency matrix of the Alarm DAG (\code{true.net}) and the related measured dataset (\code{dataset}). See the vignette for an utilization of the dataset
#' @name alarm
#' @docType data
#' @keywords data
#' @references  Aliferis C, Statnikov A, Tsamardinos I, Mani S, Koutsoukos X Local Causal and Markov Blanket Induction for Causal Discovery and Feature Selection for Classif ication Part II: Analysis and Extensions' by ; JMLR 2010'
#' @export
NULL

#' Benchmark alarm
#'@title Alarm benchmark
#'@description contains the adjacency matrix of the Alarm DAG (\code{true.net}) and the related measured dataset (\code{dataset}). See the vignette for an utilization of the dataset
#' @name alarm
#' @docType data
#' @keywords data
#' @references  Aliferis C, Statnikov A, Tsamardinos I, Mani S, Koutsoukos X Local Causal and Markov Blanket Induction for Causal Discovery and Feature Selection for Classif ication Part II: Analysis and Extensions' by ; JMLR 2010'
#' @export
NULL


#' Adjacency matrix of the Alarm benchmark
#'@title Adjacency matrix of the Alarm dataset
#'@description contains the adjacency matrix of the Alarm DAG. See the vignette for an utilization of the dataset
#' @name true.net
#' @docType data
#' @keywords data
#' @references  Aliferis C, Statnikov A, Tsamardinos I, Mani S, Koutsoukos X Local Causal and Markov Blanket Induction for Causal Discovery and Feature Selection for Classif ication Part II: Analysis and Extensions' by ; JMLR 2010'
#' @export
NULL

#'  Dataset of the Alarm benchmark
#'@title Dataset of the Alarm benchmark
#'@description contains the  measured dataset. See the vignette for an utilization of the dataset
#' @name dataset
#' @docType data
#' @keywords data
#' @references  Aliferis C, Statnikov A, Tsamardinos I, Mani S, Koutsoukos X Local Causal and Markov Blanket Induction for Causal Discovery and Feature Selection for Classif ication Part II: Analysis and Extensions' by ; JMLR 2010'
#' @export
NULL