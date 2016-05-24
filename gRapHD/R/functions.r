################################################################################
#   This file is part of gRapHD R package.
#
#   gRapHD R package
#   Copyright (C) 2009 Gabriel Coelho Goncalves de Abreu, Rodrigo Labouriau,
#   and David Edwards
#
#   gRapHD R package is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the Free
#   Software Foundation, either version 3 of the License, or any later version.
#
#   gRapHD R package program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

################################################################################
# Functions in this file:
#
#  SubGraph <- function(model=NULL,edges=NULL,v=NULL,p=0)
#  sortMat <- function(mat, cols)
#  MCS <- function(model=NULL,edges=NULL,v=0,p=NULL)
#  DFS <- function(model=NULL,edges=NULL, v, p=NULL)
#  minForest <- function(dataset,homog=TRUE,forbEdges=NULL,stat="BIC")
#  randTree <- function(p,seed=1)
#  chStat <- function(model,dataset,previous=NULL,forbEdges=NULL)
#  stepw <- function(model,dataset,stat="BIC",saveCH=NULL,forbEdges=NULL,
#                    exact=FALSE,initial=NULL,threshold=0,join=FALSE)
#  neighbours <- function(model=NULL,edges=NULL,v)
#  findEd <- function(edges,p,previous=NULL,numCat,from=0,exact=FALSE,join=F)
#  perfSets <- function(model=NULL,edges=NULL,p=NULL,varType=0,from=0)
#  rowProds <- function(x,na.rm=TRUE)
#  calcStat <- function(dataset,homog=TRUE,forbEdges=NULL,stat="LR")
#  plot.gRapHD <- function(x,vert=NULL,numIter=50,main="",
#                          plotVert=TRUE,labelVert=TRUE,energy=FALSE,
#                          useWeights=FALSE,vert.hl=NULL,col.hl="red",
#                          vert.radii=0.01,coord=NULL,col.ed="darkgray",lty.ed=1,
#                          lwd.ed=1,lwd.vert=1,border=0,symbol.vert=1,
#                          cex.vert.label=.40,vert.labels=NULL,asp=NA,disp=TRUE,
#                          font=par("font"),...)
#  neighbourhood <- function(model=NULL,edges=NULL,orig=NULL,rad=1)
#  adjMat <- function(model=NULL,edges=NULL,p=NULL)
#  Degree <- function(model=NULL,edges=NULL,v=NULL)
#  shortPath <- function(model=NULL, edges=NULL, v=NULL, p=NULL)
#  convData <- function(dataset)
#  fit <- function(model=NULL,edges=NULL,dataset,homog=TRUE)
#  summary.gRapHD <- function(object,...)
#  modelFormula <- function(model)
#  modelDim <- function(model)
#  matrix.gRapHD <- function(from)
#  gRapHD.graphNEL <- function(from)
#  graphNEL.gRapHD <- function(from)
#  print.gRapHD <- function(x,...)
#  show.gRapHD <- function(object,...)
#  is.gRapHD <- function(object)
#  jTree <- function(model)
#  ccoeff <- function(model=NULL,edges=NULL,p=NULL)
#  CI.test <- function(x,y,S,dataset,homog=TRUE)
#  seqLevels <- function(k,numCat)
#  normDens <- function (x, meanVec, covMat, logx = FALSE)
#  convertClass <- function(object)
################################################################################

################################################################################
# gRapHD class:
#   edges - integer, matrix with 2 columns, each row representing one edge,
#           and each column one of the vertices in the edge. Column 1 contains
#           the vertex with lower index.
#   homog - logical, TRUE if the covariance is homogeneous.
#   minForest - integer, first and last edges found with minForest.
#   numCat - integer, vector with number of levels for each variable (0 if
#            continuous).
#   numP - integer, vector with number of estimated parameters for each edge.
#   p - integer, number of variables (vertices) in the model.
#   stat.minForest - character, measure used (LR, AIC, or BIC).
#   stat.stepw - character, measure used (LR, AIC, or BIC).
#   stat.user - character, user defined
#   statSeq - numeric, vector with value of stat.minForest for each edge.
#   stepw - interger, first and last edges found with stepw.
#   userDef - integer, first and last edges defined by the user.
#   vertNames - character, vector with vertices' names.
################################################################################

################################################################################
# Description: Builds a sub-graph based on a list of vertices.
# In: model = gRapHD object
#     edges = matrix (n by 2, numeric).
#     v = vector k (vertices)
#     p = number of vertices
# Out: gRapHD object
################################################################################
SubGraph <- function(model=NULL,edges=NULL,v=NULL,p=0)
{
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
  {
    if (class(model)!="gRapHD")
      stop("Model must be of class gRapHD.")
    result <- model
    edges <- model@edges
  }
  else
  {
    p <- ifelse(NROW(edges)==0,p,max(max(edges),p))
    result <- new("gRapHD",edges=edges,p=p)
  }

  if (NROW(edges)==0)
    return(new("gRapHD",p=result@p,numCat=result@numCat,homog=result@homog))

  if (!is.null(v) & length(v)!=0)
  {
    v <- sort(unique(v))
    if (max(v)>result@p)
      stop(paste("v must be in [1,",result@p,"].",sep=""))
    ind <- (edges[,1]%in%v) & (edges[,2]%in%v)
    result@edges <- matrix(edges[ind,],ncol=2)
    result@statSeq <- result@statSeq[ind]
    result@numP <- result@numP[ind]
    result@vertNames <- result@vertNames[v]
    result@numCat <- result@numCat[v]
    aux <- rep(0,result@p)
    aux[v] <- 1:length(v)
    result@edges <- cbind(aux[result@edges[,1]],aux[result@edges[,2]])
    rm(aux)
    ind <- (1:nrow(edges))[ind]
    # if minForest was used, the resulting edges are the first in the list
    # and after come the edges from stepw
    if (sum(result@minForest)!=0)
    {
      if (length(ind)==0)
        result@minForest <- integer(2)
      else
        if (ind[1]>result@minForest[2])
          result@minForest <- integer(2)
        else
        {
          x <- which(ind<=result@minForest[2])
          x <- x[length(x)]
          result@minForest <- as.integer(c(1,x))
        }
    }
    if (sum(result@stepw)!=0)
    {
      if (length(ind)==0)
        result@stepw <- integer(2)
      else
        if (ind[length(ind)]<result@stepw[1])
          result@stepw <- integer(2)
        else
        {
          x <- which(ind<=result@stepw[2])
          x <- x[length(x)]
          result@stepw <- as.integer(c(x,length(ind)))
        }
    }
    if (sum(result@userDef)!=0)
      if (length(ind)>0)
        result@userDef <- as.integer(c(1,length(ind)))
    result@p <- length(v)
  }
  return(result)
}
################################################################################

################################################################################
# Sorts the rows of a matrix, using the cols (in sequence).
# In: mat = matrix
#     cols = vector numeric (sequence of columns to be used in the sorting)
# Out: sorted matrix
################################################################################
sortMat <- function(mat, cols)
{
  m <- do.call("order", as.data.frame(mat[, cols]))
  mat[m, ]
}
################################################################################

################################################################################
# Maximum Cardinality Search (Tarjan & Yannakakis, in Jayson Rome, 2002).
# Searches for a perfect numbering.
# In: model = gRapHD object
#     edges = matrix (n by 2, numeric).
#     v     = starting vertex;
#     p     = number of vertices.
# Out: vector = (p by 1) where each position represents a
#               variable, containing the respective number in the perfect
#               numbering; os 0 is the graph is not triangulated.
################################################################################
MCS <- function(model=NULL,edges=NULL,v=0,p=NULL)
{
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
  {
    if (class(model)!="gRapHD")
      stop("Model must be of class gRapHD.")
    edges <- model@edges
    p <- model@p
  }
  else
    v <- max(c(v,1))  

  if (NROW(edges)==0 && is.null(p))
    stop("p must not be NULL.")
  if (is.null(p))
    p <- max(edges)
  v <- ifelse(v==0,1,v)
  if (v > p)
    stop(paste("v must be in [1,",p,"].",sep=""))
  v1 <- edges[,1]
  v2 <- edges[,2]
  storage.mode(v1) <- storage.mode(v2) <- "integer"
  storage.mode(v) <- storage.mode(p) <- "integer"
  result <- .Call("mcs",v1,v2,v,p,PACKAGE="gRapHD")
  return(as.vector(result))
}
################################################################################

################################################################################
# Depth-first search, searches for all vertices reacheable from one specific
# vertex (considering that there is no cycle).
# In: model = gRapHD object
#     edges = matrix (n by 2, numeric).
#     v    = starting vertex;
#     p    = number of vertices.
# Out: K = vector containing all the vertices reacheable from v.
################################################################################
DFS <- function(model=NULL,edges=NULL, v, p=NULL)
{
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
  {
    if (class(model)!="gRapHD")
      stop("Model must be of class gRapHD.")
    edges <- model@edges
    p <- model@p
  }
  else
    if (NROW(edges)==0 || NCOL(edges)!=2)
      stop("No edges.")
  if (is.null(p))
    p <- max(edges)
  if (v>p)
    stop(paste("v must be in [1,",p,"].",sep=""))

  v1 <- edges[,1]
  v2 <- edges[,2]
  storage.mode(v1) <- storage.mode(v2) <- "integer"
  storage.mode(v) <- storage.mode(p) <- "integer"
  result <- .Call("dfs",v1,v2,v,p,PACKAGE="gRapHD")
  rm(v1,v2)
  return(as.vector(result))
}
################################################################################

################################################################################
# Searches for the tree that minimizes the LR.
# In: dataset = matrix n by p.
#     homog = boolean - if the covariance is homogeneous
#     forbEdges = matrix k by 2, forbidden edges
#     stat = "LR" (default), "BIC", or "AIC"
# As defined in Lauritzen (1996, pg 168): saturated mixed model - lets X have
# an arbritary CG distribution. The homogeneous saturated mixed model restricts
# the distribution of X to a CG distribution which is homogeneous, i.e. the
# conditional covariance given the discrete variables does not depend on the
# conditional.
# Out: gRapHD object
################################################################################
minForest <- function(dataset,homog=TRUE,forbEdges=NULL,stat="BIC",cond=NULL,...)
{
  comps <- function(cond)
  {
    for (i in 1:length(cond))
      storage.mode(cond[[i]]) <- "integer"
    result <- .Call("components",cond,PACKAGE="gRapHD")
    return(result)
  }

  if (length(stat)==0)
    stop("No valid measure. Options: LR, BIC, AIC, or a user defined function.")
  if (mode(stat)=="function")
  {
    FUN <- match.fun(stat)
    stat <- 3
  }
  else
  {
    stat <- toupper(stat)
    if (is.element(stat,c("LR","AIC","BIC")))
      stat <- switch(stat,LR=0,BIC=1,AIC=2)
    else
      stop("No valid measure. Options: LR, BIC, AIC, or a user defined function.")
  }

  gClique <- function(v)
  {
    v <- sort(unique(v))
    x <- function(v,y)
    {
      if (v<y)
        return(cbind(v,(v+1):y))
      else
        return(NULL)
    }
    n <- length(v)
    aux <- sapply(1:n,x,y=n)
    edges <- NULL
    for (i in 1:length(aux))
      edges <- rbind(edges,matrix(v[aux[[i]]],,2))
    return(edges)
  }

  p <- ncol(dataset)
  n <- nrow(dataset)
  statS <- stat
  aux <- convData(dataset)
  dataset <- aux$ds
  numCat <- aux$numCat
  result <- aux$vertNames
  rm(aux)

  if (sum(is.na(dataset))>0)
    stop("Missing values not allowed.")

  if (is.null(forbEdges))
    forbEdges <- 0

###################################
  condEdges <- NULL
  compType <- varType <- as.integer(numCat!=0)

  if (is.null(cond))
    comp <- 1:p
  else
  {
    k <- length(cond) #number of conditional sets
    condJ <- cond
    ed <- matrix(integer(0),,2)
    if (k > 1)
    {
      compJ <- comps(cond)
      condJ <- list() #the list of nodes in the original graph that are in the same connected component
      for (i in 1:max(compJ))
        condJ[[i]] <- sort(unique(unlist(cond[which(compJ==i)])))
    }
    comp <- rep(0,p)
    for (i in 1:length(condJ))
    {
      comp[condJ[[i]]] <- max(comp)+1
      x <- unique(varType[condJ[[i]]])
      compType[condJ[[i]]] <- ifelse(length(x)==1,x,2)
    }
    for (i in 1:length(cond))
      condEdges <- rbind(condEdges,gClique(cond[[i]]))
    condEdges <- unique(condEdges,MARGIN=1)
    ind <- which(comp==0)
    if (length(ind)>0)
    {
      comp[ind] <- max(comp) + (1:length(ind))
    }
  }
###################################
  if (forbEdges[1] != 0)
  {
    forbEdges <- matrix(forbEdges,,2)
    X <- function(x){sort(x)}
    forbEdges <- t(apply(forbEdges,1,X))
    x <- forbEdges[,1]
    y <- forbEdges[,2]
    forbEdges <- (x-1)*p-(x-1)*x/2 + y-x #simpler representation of the edges
  }
  values <- NULL
  if (stat == 3)
  {
    for (i in 1:(p-1))
      for (j in (i+1):p)
      {
        x <- (i-1)*p-(i-1)*i/2 + j-i
        if (!is.element(x,forbEdges))
          values <- rbind(values,FUN(c(i,j),numCat,dataset,...))
        else
          values <- rbind(values,NA)
      }
  }

  storage.mode(values) <- storage.mode(dataset) <- "double"
  storage.mode(numCat) <- storage.mode(homog) <- "integer"
  storage.mode(forbEdges) <- storage.mode(stat) <- "integer"
  storage.mode(comp) <- storage.mode(compType) <- "integer"
  aux <- .Call("minForest",dataset,numCat,homog,
               forbEdges,stat,values,comp,compType,PACKAGE="gRapHD")
  tree <- aux$tree
  storage.mode(homog) <- "logical"
  n <- NROW(tree)

  if (is.null(condEdges))
  {
    edges <- matrix(tree[,1:2],n,2)
    statSeq <- tree[,3]
    numP <- tree[,4]
  }
  else
  {
    edges <- rbind(matrix(tree[,1:2],n,2),condEdges)
    statSeq <- c(tree[,3],rep(NA,nrow(condEdges)))
    numP <- c(tree[,4],rep(NA,nrow(condEdges)))
  }
  result <- new("gRapHD",edges = edges,
                p = as.integer(p),
                stat.minForest = switch(stat+1,"LR","BIC","AIC","User's function"),
                statSeq = statSeq,
                numCat = as.integer(numCat),
                homog = homog,
                numP = as.integer(numP),
                vertNames = result,
                minForest = as.integer(c(1,nrow(edges))))
#  if (dim(aux$errors)[1]!=0)
#  {
#    result$error <- aux$errors
#    warning("Check model$errors for edges with problems.", call. = FALSE)
#  }
  return(result)
}
################################################################################

################################################################################
# Generates a random tree (no cycles) with "p" vertices.
# In: p = integer (number of vertices)
#     seed = for the random generator.
# Out: list = edges (matrix p-1 by 2, with edges)
#             seed
#             p
################################################################################
randTree <- function(p,seed=1)
{
  set.seed(seed,kind="Mersenne-Twister")
  V <- 1:p
  r <- sample(1:p,1)
  V[c(r,p)] <- V[c(p,r)]
  tree <- matrix(0,nrow=p-1,ncol=2)
  for (i in 1:(p-1))
  {
    x <- ifelse(i==p-1,1,sample(1:(p-i),1))
    y <- ifelse(i==1,p,sample((p-i+1):p,1))
    tree[i,] <- sort(c(V[x],V[y]))
    V[c(x,p-i)] <- V[c(p-i,x)]
  }
  result <- list(edges=tree,seed=seed,p=p)
  return(result)
}
################################################################################

################################################################################
# Calculates the change in the model's BIC resulting by the addition of all
# possible edges (edges that preserve the graph's triangulation).
# In: model = gRapHD object
#     dataset = matrix n by p.
#     previous = same output of chStat()
#     forbEdges = matrix k by 2, forbidden edges
# As defined in Lauritzen (1996, pg 168): the saturated mixed model lets X have
# an arbritary CG distribution. The homogeneous saturated mixed model restricts
# the distribution of X to a CG distribution which is homogeneous, i.e. the
# conditional covariance given the discrete variables does not depend on the
# conditional.
# Out: list:
#      edges.to.test = matrix (k by 4), column:
#                      1 - first vertex of the tested edge
#                      2 - second vertex (the two values are ordered)
#                      3 - index for the separator in S
#                      4 - change in the LR by adding this edge
#                      5 - number of parameters for that edge
#      S = list with the separators
################################################################################
chStat <- function(model,dataset,previous=NULL,forbEdges=NULL)
{
  if (sum(is.na(dataset))>0)
    stop("Missing values not allowed.")

  p <- ncol(dataset) # number of variables (vertices)
  n <- nrow(dataset) # number of observations
  edges.to.test <- previous$edges.to.test # (v1; v2; index of the minimal separator;
                             # a previous value for the LR change, or 0
                             # indicating that it has to be calculated;
                             # the number of parameters for that edge)
  SS <- previous$S # minimal separators
  rm(previous)

  num <- nrow(edges.to.test)

  numbPar <- function(y)
  {
    a <- sum(rowSums(y)!=0)-1
    b <- sum(colSums(y)!=0)-1
    return(a*b*(a>0))
  }
  whichSeq <- function(x,y){return(sum(abs(x-y))==0)}
  Gcliq <- function(v)
  {
    v <- sort(unique(v))
    i <- 1
    n <- length(v)
    ed <- sort(rep(i:(i+n-1),n))
    ed <- cbind(ed,rep(i:(i+n-1),n))
    ed <- matrix(ed[-unlist(lapply(1:n,function(x,n){1:x+(x-1)*n},n=n)),],,2)
    dimnames(ed) <- NULL
    ed[,1] <- v[ed[,1]]
    ed[,2] <- v[ed[,2]]
    return(ed)
  }

  if (num > 0)
    for (i in 1:num)
    {
      x <- edges.to.test[i,1]
      y <- edges.to.test[i,2]
      if (!is.na(edges.to.test[i,4]) && is.finite(edges.to.test[i,4]))
        if ((edges.to.test[i,4]==0) &
            (!is.element((x-1)*p-(x-1)*x/2+y-x,forbEdges)))
        { # the change has to be calculated
          # this is the clique that will result from the inclusion of this edge
          S <- SS[[edges.to.test[i,3]]]
          clique <- c(edges.to.test[i,1:2],S)

          # the 2 cliques thal will be merged
          clique1 <- c(edges.to.test[i,1],S)
          clique2 <- c(edges.to.test[i,2],S)

          if (sum(model@numCat[clique]!=0) == 0) # continuous
          {
            CM <- cov(dataset[,clique],use="pairwise.complete.obs")*(nrow(dataset)-1)
            if (sum(is.na(CM))) # meaning that the new edge does not have more
                                # than 1 observation (not NA)
            {
              edges.to.test[i,4] <- 0
              edges.to.test[i,5] <- 0
            }
            else
            {
              a <- match(edges.to.test[i,1:2],clique)
              # because of join in findEd
              d <- log(det(CM)) + ifelse(length(S)==0,0,log(det(matrix(CM[-a,-a],length(clique)-2))))
              ##
              d <- d - (log(det(matrix(CM[-a[1],-a[1]],length(clique)-1))) +
                        log(det(matrix(CM[-a[2],-a[2]],length(clique)-1))))
              edges.to.test[i,4] <- nrow(dataset)*d
              edges.to.test[i,5] <- 1
            }
          }
          else
            if (sum(model@numCat[clique]==0)==0) # discrete
            {
              t12 <- table(as.data.frame(dataset[,clique]))
              t1 <- margin.table(t12,match(clique1,clique))
              t2 <- margin.table(t12,match(clique2,clique))
              # because of join in findEd
              tS <- if(length(S)==0) 0 else margin.table(t12,match(S,clique))
              # to calculate the df: Whittaker, pg 223, Proposition 7.6.2 and
              # pg 224, Corolary 7.6.3
              # using MARGIN=(1:length(clique))[-(1:2)] uses the other margins but
              # not the first 2 (that are the 2 vertices being tested)

              # because of join in findEd
              if (length(S)==0)
                numP <- prod(model@numCat[clique]-1)#length(t12)-1
              else
                numP <- sum(apply(t12,MARGIN=(1:length(clique))[-(1:2)],numbPar))
              ##
              numObsUsed <- sum(t12)
              t12 <- t12[t12!=0]
              t1 <- t1[t1!=0]
              t2 <- t2[t2!=0]
              tS <- tS[tS!=0]

              if (numObsUsed>0)
              {
                # BIC(n+1) - BIC(n) = log(f(C)/(f(C1)*f(C2)/f(S)))
                v12 <- sum(t12*log(t12/numObsUsed))
                v1  <- sum(t1 *log(t1 /numObsUsed))
                v2  <- sum(t2 *log(t2 /numObsUsed))
                # because of join in findEd
                vS  <- ifelse(length(S)==0,0,sum(tS *log(tS /numObsUsed)))
                ##
                edges.to.test[i,4] <- -2*(v12 + vS - v1 - v2)
                edges.to.test[i,5] <- numP
              }
              else
              {
                edges.to.test[i,4] <- 0
                edges.to.test[i,5] <- 0
              }
            }
            else # mixed
            {
              discrete   <- intersect(clique,which(model@numCat!=0))
              continuous <- setdiff(clique,discrete)
              tab <- table(as.data.frame(dataset[,discrete]))
              ssd <- diag(0,length(continuous))
              if (sum(model@numCat[edges.to.test[i,1:2]])==0) # both are continuous
                if (model@homog)
                {
                  for (j in 1:length(tab))
                    if (tab[j]>0)
                    {
                      ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,model@numCat[discrete]))
                      if (sum(ind)>1) #otherwise there is only one observation, and the contribution to the variance is zero
                        ssd <- ssd + cov(dataset[ind,continuous],use="pairwise.complete.obs")*(tab[j]-1)
                    }
                  ind <- match(edges.to.test[i,1:2],continuous)
                  edges.to.test[i,4] <- n*(log(det(ssd)) +
                                            ifelse(length(continuous)>2,log(det(matrix(ssd[-ind,-ind],length(continuous)-2))),0) -
                                            log(det(matrix(ssd[-ind[1],-ind[1]],length(continuous)-1))) -
                                            log(det(matrix(ssd[-ind[2],-ind[2]],length(continuous)-1))))
                  edges.to.test[i,5] <- 1
                }
                else
                {
                  edges.to.test[i,4] <- 0
                  ind1 <- match(edges.to.test[i,1:2],continuous)
                  for (j in 1:length(tab))
                    if (tab[j]>0)
                    {
                      ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,model@numCat[discrete]))
                      if (sum(ind)>1) #otherwise there is only one observation, and the contribution to the variance is zero
                        ssd <- cov(dataset[ind,continuous],use="pairwise.complete.obs")*(tab[j]-1)
                      edges.to.test[i,4] <- edges.to.test[i,4] +
                                            +tab[j]*(log(det(ssd)) +
                                                    ifelse(length(continuous)>2,log(det(matrix(ssd[-ind1,-ind1],length(continuous)-2))),0) -
                                                    log(det(matrix(ssd[-ind1[1],-ind1[1]],length(continuous)-1))) -
                                                    log(det(matrix(ssd[-ind1[2],-ind1[2]],length(continuous)-1))))
                    }
                  edges.to.test[i,5] <- sum(tab>(length(continuous)+1))
                }
              else
              {
                vDiscr <- (edges.to.test[i,1:2])[model@numCat[edges.to.test[i,1:2]]!=0]
                vCont  <- setdiff(edges.to.test[i,1:2],vDiscr)
                discrMarg <- setdiff(discrete,vDiscr)
                tabMarg <- margin.table(tab,match(discrMarg,discrete))
                ssdMarg <- diag(0,length(continuous))
                if (model@homog)
                {
                  for (j in 1:length(tab))
                  {
                    if (tab[j]>0)
                    {
                      ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,model@numCat[discrete]))
                      if (sum(ind)>1) #otherwise there is only one observation, and the contribution to the variance is zero
                        ssd <- ssd + var(dataset[ind,continuous],use="pairwise.complete.obs")*(tab[j]-1)
                    }
                    if (j<=length(tabMarg))
                    {
                      if (tabMarg[j]>0)
                      {
                        if (length(discrMarg)>0)
                          ind <- apply(matrix(dataset[,discrMarg],,length(discrMarg)),1,whichSeq,y=seqLevels(j,model@numCat[discrMarg]))
                        else
                          ind <- rep(TRUE,tabMarg[j])
                        if (sum(ind)>1) #otherwise there is only one observation, and the contribution to the variance is zero
                          ssdMarg <- ssdMarg + var(dataset[ind,continuous],use="pairwise.complete.obs")*(tabMarg[j]-1)
                      }
                    }
                  }
                  ind <- match(vCont,continuous)
                  edges.to.test[i,4] <- n*(log(det(ssd)) +
                                            ifelse(length(continuous)>1,log(det(matrix(ssdMarg[-ind,-ind],length(continuous)-1))),0) -
                                            ifelse(length(continuous)>1,log(det(matrix(ssd[-ind,-ind],length(continuous)-1))),0) -
                                            log(det(ssdMarg)))
                  if (length(discrete)==1) #means that the separator is continuous
                    edges.to.test[i,5] <- model@numCat[vDiscr]-1
                  else
                  {
                    aux <- match(vDiscr,discrete) #position of the discrete variable in the new edge (in the list of discrete variables in the edge/separator)
                    # dimension 1 in aux1 is the discrete variable in the tested edge
                    aux1 <- table(as.data.frame(dataset[,c(vDiscr,discrete[-aux])])) #it's only the table with the discrete in the edge as the 1st dimension
                    aux.tab <- margin.table(aux1,2:length(discrete)) #the marginal table without the edge's discrete
                    aux.ind <- which(aux.tab>0,arr.ind=TRUE) #cells with at least one observation
                    numP <- 0
for (iii in 1:nrow(aux.ind))
  numP <- numP + (sum(aux1[cbind(1:model@numCat[vDiscr],aux.ind[iii,])]>0) - 1)
                    rm(aux,aux1,aux.tab,aux.ind)
                    edges.to.test[i,5] <- numP
                  }
                }
                else
                {
                  edges.to.test[i,4] <- 0
                  ind1 <- match(vCont,continuous)
                  for (j in 1:length(tab))
                  {
                    if (tab[j]>0)
                    {
                      ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,model@numCat[discrete]))
                      ssd <- var(dataset[ind,continuous],use="pairwise.complete.obs")*(tab[j]-1)
                      dim(ssd) <- rep(length(continuous),2)
                    }
                    if (j<=length(tabMarg))
                    {
                      if (tabMarg[j]>0)
                      {
                        if (length(discrMarg)>0)
                          ind <- apply(matrix(dataset[,discrMarg],,length(discrMarg)),1,whichSeq,y=seqLevels(j,model@numCat[discrMarg]))
                        else
                          ind <- rep(TRUE,tabMarg[j])
                        ssdMarg <- var(dataset[ind,continuous],use="pairwise.complete.obs")*(tabMarg[j]-1)
                        dim(ssdMarg) <- rep(length(continuous),2)
                      }
                    }
                    ind <- match(vCont,continuous)
                    edges.to.test[i,4] <- edges.to.test[i,4] +
                                          -tab[j]*(log(det(ssd))-ifelse(length(continuous)>1,log(det(matrix(ssd[-ind,-ind],length(continuous)-1))),0)) +
                                          +ifelse(j<=length(tabMarg),tabMarg[j],0)*(log(det(ssdMarg))-ifelse(length(continuous)>1,log(det(matrix(ssdMarg[-ind,-ind],length(continuous)-1))),0))
                  }
                  edges.to.test[i,4] <- edges.to.test[i,4] +
                                        sum(tab*log(tab),na.rm=TRUE) -
                                        sum(tabMarg*log(tabMarg),na.rm=TRUE)
                  edges.to.test[i,4] <- -edges.to.test[i,4]
                  if (length(discrete)==1) #means that the separator is continuous
                    edges.to.test[i,5] <- sum(tab>(length(continuous)+1))-1
                  else
                  {
                    aux <- match(vDiscr,discrete)
                    # dimension 1 in aux1 is the discrete variable in the tested edge
                    aux1 <- table(as.data.frame(dataset[,c(vDiscr,discrete[-aux])]))
                    aux.tab <- margin.table(aux1,2:length(discrete))
                    aux.ind <- which(aux.tab>(length(continuous)+1),arr.ind=TRUE)
                    numP <- 0
for (iii in 1:nrow(aux.ind))
  numP <- numP + (sum(aux1[cbind(1:model@numCat[vDiscr],aux.ind[iii,])]>0) - 1)
                    rm(aux,aux1,aux.tab,aux.ind)
                    edges.to.test[i,5] <- numP
                  }
                }
              }
            }
        }
    }
  return(list(edges.to.test=edges.to.test,S=SS))
}
################################################################################

################################################################################
# Does a stepwise forward, starting from model and using the change in BIC to
# update the model. Tests every possible edge and selects the one that reduces
# the BIC the most. Stops when no more redution is possible, or there is no
# more possible edge to test.
# In: model = gRapHD object
#     dataset = matrix n by p
#     stat = "BIC","AIC","LR"
#     saveCH = if not NULL save each iteration's "chStat", with saveCH as
#              part of the file name
#     forbEdges = matrix k by 2, forbidden edges
#     exact = TRUE or FALSE
#     initial = continue the algorithm from a previous point. The paramenter
#                must have the same structure as the result of chStat.
#     threshold = values greater than the threshold are not included in the
#                 model
# Out: model (the same updated object)
################################################################################
stepw <- function(model,dataset,stat="BIC",saveCH=NULL,forbEdges=NULL,
                  exact=FALSE,initial=NULL,threshold=0,join=FALSE)
{
  if (class(model)!="gRapHD")
    stop("Model is not of class gRapHD.")
  else
  {
    star <- NULL
    aux <- range(model@numCat)
    if ((aux[1]==0) & (aux[2]!=0)) #mixed
    {
      # there are at least 2 discrete variables, which means that forbidden paths
      # are possible to happen
      # star graph
      star <- cbind(which(model@numCat!=0),model@p+1)
      perfect <- perfSets(edges=rbind(model@edges,star),
                          p=model@p+1,varType=c(model@numCat,1),from=0)
      if (!is.list(perfect))
        stop("The model is not strongly decomposable.")
      rm(perfect)
    }
    else
      if (MCS(edges=model@edges,,p=model@p)[1]==0)
        stop("The model is not triangulated.")
    rm(aux)
  }
  if (length(stat)==0)
    stop("No valid measure. Options: LR, BIC, AIC, or a user defined function.")
  if (mode(stat)=="function")
  {
    FUN <- match.fun(stat)
    model@stat.stepw <- "User's function"
    stat <- "USER"
  }
  else
  {
    stat <- toupper(stat)
    if (!is.element(stat,c("LR","AIC","BIC")))
      stop("No valid measure. Options: LR, BIC, AIC, or a user defined function.")
    model@stat.stepw <- stat
  }

  edges.stepw <- nrow(model@edges)+1

  if (model@p != NCOL(dataset))
    stop("model@p doesn't agree with NCOL(dataset).")

  aux <- convData(dataset)
  dataset <- aux$ds
  numCat <- aux$numCat
  if (!all.equal(numCat,model@numCat))
    stop("model@numCat doesn't agree with data.")
  if (!identical(sort(aux$vertNames),sort(model@vertNames)))
    stop("model@vertNames doesn't agree with data.")
  dataset <- dataset[,match(model@vertNames,aux$vertNames)]
  rm(aux)
  if (sum(is.na(dataset))>0)
    stop("Missing values not allowed.")

  p <- ncol(dataset) # number of variables (vertices)
  n <- nrow(dataset) # number of observations
  stat <- toupper(stat)

  if (is.null(forbEdges))
    forbEdges <- 0
  else
  {
    forbEdges <- matrix(forbEdges,,2)
    X <- function(x){sort(x)}
    forbEdges <- t(apply(forbEdges,1,X))
    x <- forbEdges[,1]
    y <- forbEdges[,2]
    forbEdges <- (x-1)*p-(x-1)*x/2 + y-x; #simpler representation of the edges
    rm(X,x,y)
  }

  STOP <- FALSE # flag indicating when to stop
  iteration <- 0
  if (!is.null(model@minForest))
    iteration <- nrow(model@edges) - model@minForest[2]
  ch <- initial
  SS <- NULL

  while ((!STOP))
  {
    iteration <- iteration + 1
    ch <- findEd(model@edges,p,ch,model@numCat,0,exact,join)
    if (stat!="USER")
      ch <- chStat(model,dataset,ch,forbEdges)
    else
      ch <- FUN(model,dataset,ch,forbEdges)
    change <- ch$edges.to.test
    if (!is.null(saveCH))
    {
      fileN <- paste(saveCH,"_",
                     paste(rep("0",6-(as.integer(log10(iteration))+1)),sep="",
                           collapse=""),iteration,".RData",sep="")
      save(ch,model,file=fileN)
    }
    if (nrow(change) == 0) # there is no more edges possible to test
      STOP <- TRUE
    else
    {
      continue <- TRUE
      add <- FALSE
      # before adding one the edge, it has to test if it generates a cycle
      df <- change[,5]
      df[df<=0] <- NA
      statValues <- change[,4] + ((stat=="BIC")*df+(stat=="AIC")*2)*log(n)
      rm(df)
      ind <- order(statValues)
      i <- 0
      while (continue)
      {
        # the edge which reduces the most (the stat)
        i <- i+1
        aux <- ind[i]
        value <- statValues[aux]
        if (is.na(value))
          continue <- FALSE
        else
          if (is.finite(value))
          {
            if ((value < threshold) & (exact))
            {
              continue <- FALSE
              add <- TRUE
            }
            else
              if (value < threshold)
              {
                # Lauritzen (1996, pg 11): Proposition 2.6 (Leimer) - An undirected,
                # marked graph G is decomposable if and only if G* is triangulated.
                continue <- (MCS(edges=rbind(model@edges,change[aux,1:2],star),v=1,p=p+!is.null(star)))[1] == 0
                add <- !continue #statValues[aux] <- !continue*statValues[aux] + continue*threshold
              }
              else
                continue <- FALSE
          }
      }
      if (add) # it improves the model
      {
        model@edges <- rbind(model@edges,change[aux,1:2])
        model@numP <- as.integer(c(model@numP,change[aux,5]))
        model@statSeq <- c(model@statSeq,change[aux,4])
      }
      else
        STOP <- TRUE
    }
  }
  if (edges.stepw<=nrow(model@edges))
    model@stepw <- as.integer(c(edges.stepw,nrow(model@edges)))
  else
    model@stepw <- integer(2)
  return(model)
}
################################################################################

################################################################################
# Finds the first neighbours of a vertex.
# In: model = gRapHD object
#     edges = matrix (n by 2, numeric).
#     v
# Out: fn = vector (first neighbours numbers, sortered)
################################################################################
neighbours <- function(model=NULL,edges=NULL,v)
{
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
  {
    if (class(model)!="gRapHD")
      stop("Model must be of class gRapHD.")
    edges <- model@edges
    p <- model@p
  }
  else
    if (NROW(edges)==0 || NCOL(edges)!=2)
      stop("No edges.")
    else
      p <- max(edges)
  if (v > p)
    stop(paste("v must be in [1,",p,"].",sep=""))  
  fn <- integer(0)
  edges <- matrix(edges,ncol=2)
  ind <- as.matrix(which(rbind(edges,c(0,0))==v,arr.ind=TRUE)) # edges in "v"
  n <- nrow(ind) # number of edges found
  if (n > 0)
    for (i in 1:n)
      fn <- c(fn,edges[ind[i,"row"],(1:2)[-ind[i,"col"]]])
  rm(ind,n)
  return(sort(fn))
}
################################################################################

################################################################################
# Given a triangulated graph (edges), get the sets:
#               C1,...,Ck: clique ordering (from a perfect numbering)
#               H1,...,Hk: histories
#               R1,...,Rk: residues
#               S1,...,Sk: separators
# From these sets, finds the possible edges (but here, it finds some edges that
# don't preserve the triangulation, so, before including a edge, it has to be
# tested). It also gives the minimal separator, what makes easier (and much
# faster) the calculations.
# If some previous data is provided, the last column in the resulting matrix
# indicates if the value should be recalculated or not.
# In: edges = matrix k by 2
#     p = integer (number of variables in the dataset)
#     previous = list: edges.to.test = matrix w by 4 (vertex 1; vertex 2; index
#                        for the separator in S; value of the change in the LR,
#                        number of parameters for that edge)
#                S: list with the minimal separators
#     varType = vertices type
#     from = first vertex in the perfect sequence
#     exact = TRUE or FALSE (if it's 1 (or TRUE) the QR decomp is used, if 2,
#             the gauss elimination)
#     join = T or F (is disjoind components can be joined)
# Out: list: edges = matrix w1 by 4 (same structure as edges.to.test)
#            S: list with the minimal separators
################################################################################
findEd <- function(edges,p,previous=NULL,numCat,from=0,exact=FALSE,join=FALSE)
{
  if (from > p)
    stop("from > p!!")
  if (is.null(edges))
    edges <- matrix(integer(0),1,2)
  if (is.null(previous))
    prev <- 0
  else
  {
    x <- previous$edges.to.test[,1]
    y <- previous$edges.to.test[,2]
    prev <- (x-1)*p-(x-1)*x/2 + y-x #simpler representation of the edges
  }

  v1 <- edges[,1]
  v2 <- edges[,2]
  storage.mode(v1) <- storage.mode(v2) <- "integer"
  storage.mode(p) <- storage.mode(prev) <- "integer"
  storage.mode(numCat) <- "integer"
  storage.mode(exact) <- storage.mode(from) <- "integer"
  storage.mode(join) <- "integer"

  result <- .Call("findEd", v1, v2, p, prev, numCat,
                  from, exact, join,PACKAGE="gRapHD")

  i <- result$total
  if (i > 0)
  {
    result <- list(edges.to.test=result$edges[1:result$total,1:5],S=result$S)
    dim(result$edges.to.test) <- c(i,5)
    result$edges.to.test[,5] <- 0
    if (!is.null(previous))
    {
      result$edges.to.test[result$edges.to.test[,4]>0,5] <-
            previous$edges.to.test[result$edges.to.test[,4],5]
      result$edges.to.test[result$edges.to.test[,4]>0,4] <-
            previous$edges.to.test[result$edges.to.test[,4],4]
    }
  }
  else
  {
    result$edges.to.test <- integer(0)
    dim(result$edges.to.test) <- c(0,5)
  }

  return(result)
}
################################################################################

################################################################################
# Given a triangulated graph (edges), get the sets:
#               C1,...,Ck: clique ordering (from a perfect numbering)
#               H1,...,Hk: histories
#               R1,...,Rk: residues
#               S1,...,Sk: separators
# In: model = gRapHD object
#     edges = matrix (n by 2, numeric).
#     p = integer (number of variables in the dataset)
#     varType = vertices type
#     from = first vertex in the perfect sequence
# Out: list - C,S,H,R
################################################################################
perfSets <- function(model=NULL,edges=NULL,p=NULL,varType=0,from=0)
{
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
  {
    if (class(model)!="gRapHD")
      stop("Model must be of class gRapHD.")
    edges <- model@edges
    p <- model@p
    varType <- model@numCat
    varType[varType!=0] <- 1
  }
  if ((NROW(edges)==0) && is.null(p))
    stop("p must not be NULL.")
  else
    if (is.null(p))
      p <- max(edges)

  from <- max(c(0,from))
  if (is.null(edges))
    edges <- matrix(integer(0),,2)

  if (!is.element(length(varType),c(1,p)))
    stop("varType must be 0, 1, or lenght(varType)==p.")

  if (length(varType)==1)
    varType <- rep(as.integer(varType!=0),p)
  else
    varType[varType!=0] <- 1

  if (nrow(edges) > 0)
    if (max(edges)>p)
      stop("There are more vertices in edges than p.")
  v1 <- edges[,1]
  v2 <- edges[,2]
  storage.mode(v1) <- storage.mode(v2) <- "integer"
  storage.mode(p) <- storage.mode(varType) <- "integer"
  storage.mode(from) <- "integer"
  result <- .Call("perfSets",v1,v2,p,varType,from,PACKAGE="gRapHD")
  
  if (storage.mode(result)=="list")
    names(result$cliques) <- names(result$histories) <-
    names(result$separators) <- names(result$residuals) <- 1:length(result$cliques)

  return(result)
}
################################################################################

################################################################################
# Like "rowSums", evaluate the products of each row.
# In: x = matrix k by p
#     na.rm = TRUE/FALSE
# Out: vector k
################################################################################
rowProds <- function(x,na.rm=TRUE)
{
  if (is.data.frame(x))
      x <- as.matrix(x)
  if (!is.array(x) || length(dn <- dim(x)) < 2)
      stop("'x' must be an array of at least two dimensions")
      
  k <- nrow(x)
  p <- ncol(x)
  storage.mode(x) <- "double"
  storage.mode(k) <- "integer"
  storage.mode(na.rm) <- storage.mode(p) <- "integer"
  result <- .Call("rowProds",x,k,p,na.rm,PACKAGE="gRapHD")
  return(result)
}
################################################################################

################################################################################
# Calculates statistics (LR for independence, AIC, BIC) for each pair of
# variables in the dataset.
# In: dataset (n by p, numeric).
#     homog = T or F
#     forbEdges = list of forbidden edges
#     stat = LR, BIC, AIC
# Out: matrix = p*(p-1)/2 by 3
#               edge - columns 1 and 2
#               stat - column (2*log(L))
################################################################################
calcStat <- function(dataset,homog=TRUE,forbEdges=NULL,stat="LR")
{
  if (length(stat)==0)
    stop("No valid measure. Options: LR, BIC, AIC, or a user defined function.")
  if (mode(stat)=="function")
  {
    FUN <- match.fun(stat)
    stat <- 3
  }
  else
  {
    stat <- toupper(stat)
    if (is.element(stat,c("LR","AIC","BIC")))
      stat <- switch(stat,LR=0,BIC=1,AIC=2)
    else
      stop("No valid measure. Options: LR, BIC, AIC, or a user defined function.")
  }

  aux <- convData(dataset)
  dataset <- aux$ds
  numCat <- aux$numCat
  rm(aux)
  if (sum(is.na(dataset))>0)
    stop("Missing values not allowed.")

  p <- ncol(dataset)
  n <- nrow(dataset)

  if (is.null(forbEdges))
    forbEdges <- 0
  else
  {
    forbEdges <- matrix(forbEdges,,2)
    X <- function(x){sort(x)}
    forbEdges <- t(apply(forbEdges,1,X))
    x <- forbEdges[,1]
    y <- forbEdges[,2]
    forbEdges <- (x-1)*p-(x-1)*x/2 + y-x; #simpler representation of the edges
  }

  values <- NULL
  if (stat == 3)
  {
    for (i in 1:(p-1))
      for (j in (i+1):p)
      {
        x <- (i-1)*p-(i-1)*i/2 + j-i
        if (!is.element(x,forbEdges))
          values <- rbind(values,FUN(c(i,j),numCat,dataset))
        else
          values <- rbind(values,NA)
      }
  }

  storage.mode(values) <- storage.mode(dataset) <- "double"
  storage.mode(numCat) <- storage.mode(homog) <- "integer"
  storage.mode(forbEdges) <- storage.mode(stat) <- "integer"
  stat <- .Call("calcStat",dataset,numCat,homog,
                forbEdges,stat,values,PACKAGE="gRapHD")
  if (dim(stat$errors)[1]!=0)
    warning("Check $errors for edges with problems.", call. = FALSE)
  return(stat)
}
################################################################################

################################################################################
# Plot the graph using Fruchterman-Reingold algorithm.
# Fruchterman, T. M. J., & Reingold, E. M. (1991). Graph Drawing by
# Force-Directed Placement. Software: Practice and Experience, 21(11).
# In: x = gRapHD object
#     vert = list of vertices to be plotted
#     numIter = number of iterations for Fruchterman-Reingold algorithm
#     main = main title
#     plotVert = T or F (add the vertices)
#     energy = T or F (use the minimum energy as initial values)
#     useWeights = T or F (use the model@statSeq as edges lenght weights)
#     vert.hl = list of vertices to highlight
#     col.hl = colour to be used in the highlighted vertices
#     vert.radii = radii of the edges (scalar or vector with length p)
#     coord = matrix (number of vertices by 2) with initial coordenates
#     col.ed = colour of the edges
#     lty.ed = type of line to be used for the edges
#     lwd.ed = width of the edges
#     lwd.vert = width of the vertices' border
#     border = colour of the vertices borders (length 1 or number of vertices)
#     symbol.vert = symbol to be used in the vertices (length 1 or number of
#          vertices) 1 is a circle, 2 a square, 3 or higher represents the
#          number of sides
#     cex.vert.label  = numeric character expansion factor for the labels;
#            multiplied by par yields the final character size. NULL and NA are
#            equivalent to 1.0.
#     vert.labels = labels to be used in the vertices.
#     asp = the aspect ratio y/x (see plot.window for more details).
#     disp = plot (TRUE) or not (FALSE)
#
# Out: matrix (2xp) with the coordinates of the vertices
################################################################################
plot.gRapHD <- function(x,vert=NULL,numIter=50,main="",
                        plotVert=TRUE,plotEdges=TRUE,energy=FALSE,
                        useWeights=FALSE,vert.hl=NULL,col.hl="red",
                        vert.radii=0.01,coord=NULL,col.ed="darkgray",lty.ed=1,
                        lwd.ed=1,lwd.vert=1,border=0,symbol.vert=1,
                        cex.vert.label=.40,vert.labels=TRUE,asp=NA,disp=TRUE,
                        font=par("font"),col.labels=NULL,add=FALSE,...)
{
  model <- x
  rm(x)
  edges <- model@edges
  varType <- model@numCat
  p <- model@p
  originalOrder <- 1:p
  if (is.null(vert))
    vertices <- 1:p
  else
  {
    if ((max(vert) > p) || (min(vert) < 1))
      stop("vert has vertices that are not in the model.")
    if (length(vert)!=length(unique(vert)))
      stop("There are repeated vertices in vert.")
    originalOrder <- order(vert)
    vertices <- sort(vert)
#
    edges <- SubGraph(edges=edges,v=vertices,p=max(c(max(edges),vertices)))@edges
    #edges <- SubGraph(edges=edges,v=vertices,p=max(edges))@edges
#
    edges[,1] <- vertices[edges[,1]]
    edges[,2] <- vertices[edges[,2]]
  }
  p <- length(vertices)

  # vertices radii
  vert.radii <- rep(vert.radii,length.out=p)
  # test if all vertices to be highlighted are in the graph
  if (sum(is.element(vert.hl,vertices))<length(vert.hl))
     stop("Only vertices plotted can be highlighted.")
  # colour of vertices' borders
  border <- rep(border,length.out=p)
  # shape of the vertices
  symbol.vert <- rep(symbol.vert,length.out=p)

  # initialise vertices' coordenates
  if (!is.null(coord)) # if it is given by the user
  {
    if (NROW(coord)!=length(vertices))
      stop("coord must have the same number of rows as vertices.")
    w <- p^2/4
    coord <- w*(2*coord-1) # the algorithm uses the graph centered in (0,0)
  }
  else
    if (!energy) # or using a lattice
    {
      w <- floor(sqrt(p))
      x <- seq(from=-p^2/4,to=p^2/4,by=(2*p^2/4)/w)
      y <- seq(from=p^2/4,to=-p^2/4,by=-(2*p^2/4)/w)
      x <- matrix(rep(x,w+1),nrow=w+1,byrow=T)
      y <- matrix(rep(y,w+1),nrow=w+1)
      x <- t(x)
      y <- t(y)
      coord <- cbind(as.vector(x),as.vector(y))
      coord <- coord[1:p,]
    }
    else # using minimum energy
    {
      if (!useWeights)
      {
        A <- adjMat(edges=edges, p=p)
        Q <- -A
        diag(Q) <- rowSums(A)
      }
      else
      {
        A <- matrix(0,p,nrow(edges))
        A[cbind(edges[,1],1:nrow(edges))] <- -1
        A[cbind(edges[,2],1:nrow(edges))] <- 1
        Q <- A%*%diag(1/abs(model@statSeq))%*%t(A)
      }
      tt <- eigen(Q, symmetric=T)
      v1 <- tt$vectors[ ,p]
      v2 <- tt$vectors[ ,p-1]
      v3 <- tt$vectors[ ,p-2]
      u1 <- v1
      u2 <- v2 - sum(v2*u1)*u1
      u3 <- v3 - sum(v3*u1)*u1 - sum(v3*u2)*u2
      coord <- cbind(u3,u2)+cbind(1:p,1:p)/(1000*p)
    }
  aux <- rep(0,max(vertices))
  aux[vertices] <- 1:length(vertices)
  edgesL <- cbind(aux[edges[,1]],aux[edges[,2]])
  storage.mode(p) <- "integer"
  storage.mode(edgesL) <- "integer"
  storage.mode(numIter) <- "integer"
  storage.mode(coord) <- "double"
  coordV <- .Call("frucRein",edgesL,p,numIter,coord,PACKAGE="gRapHD")
  # rescale coordenates resulting from frucRein
  if (p > 1)
  {
    coordV[,1] <- ((coordV[,1]-min(coordV[,1]))/max(coordV[,1]-min(coordV[,1])))
    coordV[,2] <- ((coordV[,2]-min(coordV[,2]))/max(coordV[,2]-min(coordV[,2])))
    coordV[is.na(coordV)] <- .5
  }
  else
    coordV[1,] <- c(.5,.5*.9)

  if (disp) # if it is to be plotted
  {
    if (is.null(main)) main <- ""
    save.mar <- par("mar")
    save.oma <- par("oma")
    if (!add)
    {
      plot.new()
      par(mar=c(0,0,1.5*(main!=""),0),oma=c(0,0,0,0))
      plot.window(c(0,1),c(0,1),asp=asp)
    }
    arg <- list(...)
    if (!is.null(arg$xpd)) par(xpd=arg$xpd)
    rm(arg)
    coordE <- matrix(0,nrow(edges),4)
    lty.ed <- rep(lty.ed,length.out=length(edgesL))
    lwd.ed <- rep(lwd.ed,length.out=length(edgesL))
    col.ed <- rep(col.ed,length.out=length(edgesL))
    if ((plotEdges) & (length(edgesL)>0))
      for (i in 1:nrow(edgesL))
      {
        coordE[i,1:2] <- c(coordV[edgesL[i,1],1],coordV[edgesL[i,2],1])
        coordE[i,3:4] <- c(coordV[edgesL[i,1],2],coordV[edgesL[i,2],2])
        lines(x=coordE[i,1:2],
              y=coordE[i,3:4],
              col=col.ed[i],lty=lty.ed[i],lwd=lwd.ed[i])#col="darkgray",lty=1,lwd=1)
      }

    if (plotVert)
    {
      xy <- list(x=cos(seq(0,2*pi,2*pi/25)),y=sin(seq(0,2*pi,2*pi/25)))
      symbol.vert <- symbol.vert[originalOrder]
      if (!is.na(asp))
        vert.radii <- vert.radii*ifelse(asp<1,1,asp)
      vert.radii <- vert.radii[originalOrder]
      border <- border[originalOrder]
      Fill <- rep("lightgray",length(varType))
      Fill[varType!=0] <- "black"
      Fill[vert.hl] <- col.hl
      Fill <- Fill[vertices]
      for (j in 1:nrow(coordV))
      {
        if (symbol.vert[j]==0)
          polygon(xy$x*vert.radii[j]*2+coordV[j,1],xy$y*vert.radii[j]+coordV[j,2],border=border[j],col=Fill[j],lwd=lwd.vert)
        else
          if (symbol.vert[j]==1)
            symbols(x=coordV[j,1],y=coordV[j,2],circles=vert.radii[j],inches=FALSE,add=TRUE,bg=Fill[j],
                    fg=border[j],lwd=lwd.vert)
          else
            if (symbol.vert[j]==2)
              symbols(x=coordV[j,1],y=coordV[j,2],squares=vert.radii[j],inches=FALSE,add=TRUE,bg=Fill[j],
                      fg=border[j],lwd=lwd.vert)
            else
              if (symbol.vert[j]>=3)
              {
                z1 <- rep(vert.radii[j],symbol.vert[j])
                dim(z1) <- c(1,symbol.vert[j])
                symbols(x=coordV[j,1],y=coordV[j,2],stars=z1,inches=FALSE,add=TRUE,bg=Fill[j],
                        fg=border[j],lwd=lwd.vert)
              }
      }
    }
    pLabels <- if (vert.labels[1]==F) FALSE else TRUE
    if (pLabels)
    {
      if (is.null(col.labels))
      {
        Fill <- rep("white",length(varType))
        Fill[varType==0] <- "black"
        Fill <- Fill[vertices]
      }
      else
        Fill <- col.labels

      if (!is.logical(vert.labels[1]))
        text(x = coordV[, 1], y = coordV[, 2], vert.labels,
             cex = cex.vert.label, col = Fill, font = font)
      else
        text(x = coordV[, 1], y = coordV[, 2], model@vertNames[vertices],
             cex = cex.vert.label, col = Fill, font = font)
    }
    if (!is.null(main))
      title(main=main)
    par(mar=save.mar,oma=save.oma,xpd=FALSE)
  }

  invisible(coordV)
}
################################################################################

################################################################################
# Find the sub-graph with all vertices with less or equal distance from a
# central vertex.
# In: model = gRapHD object
#     edges = matrix (n by 2, numeric).
#     orig = central vertex
#     rad = distance
# Out: list - subEdges = matrix (ncol=2)
#             v = matrix (vertex and the distance to the orig)
################################################################################
neighbourhood <- function(model=NULL,edges=NULL,orig=NULL,rad=1)
{
  if (is.null(orig))
    stop("orig must be not NULL.")
    
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
  {
    if (class(model)!="gRapHD")
      stop("Model must be of class gRapHD.")
    edges <- model@edges
    p <- model@p
  }
  p <- 0
  if (NROW(edges)==0)
    return(list(edges=matrix(integer(0),,2),v=matrix(c(orig,0),,2)))
  else
    p <- max(p,max(edges))
    
  if (orig > p)
    stop(paste("orig must be in [1,",p,"].",sep=""))
  if (rad < 1)
    vertices <- orig
  else
  {
    vertices <- list()
    vertices[[1]] <- neighbours(edges=edges,v=orig) # first neighbours
    if (rad > 1)
      for (j in 2:rad)
      {
        vertices[[j]] <- integer(0)
        for (i in 1:length(vertices[[j-1]]))
          vertices[[j]] <- c(vertices[[j]],neighbours(edges=edges,v=vertices[[j-1]][i]))
        vertices[[j]] <- unique(setdiff(vertices[[j]],c(orig,unlist(vertices[1:(j-1)]))))
        if (length(vertices[[j]])==0)
        {
          vertices[[j]] <- NULL
          break
        }
      }
  }
  v <- c(orig,0)
  for (i in 1:length(vertices))
    if (length(vertices[[i]])>0)
      v <- rbind(v,cbind(vertices[[i]],i))
  colnames(v) <- c("vertex","rad")
  rownames(v) <- NULL
  V <- sort(unique(v[,1]))
  subEdges <- SubGraph(edges=edges,v=V,p=p)@edges
  subEdges[,1] <- V[subEdges[,1]]
  subEdges[,2] <- V[subEdges[,2]]
  return(list(edges=subEdges,v=v))
}
################################################################################

################################################################################
# Return the adjacency matrix.
# In: model = gRapHD object
#     edges = dataset (n by 2, numeric).
#     p = number of vertices
# Out: matrix (p by p)
################################################################################
adjMat <- function(model=NULL,edges=NULL,p=NULL)
{
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
  {
    if (class(model)!="gRapHD")
      stop("Model must be of class gRapHD.")
    edges <- model@edges
    p <- model@p
  }
  if ((NROW(edges)==0) && is.null(p))
    stop("p must not be NULL.")
  if (is.null(p))
    p <- max(edges)

  result <- matrix(0,p,p)
  result[edges] <- 1
  result[edges[,2:1]] <- 1
  return(result)
}
################################################################################

################################################################################
# Return the degree of a list of vertices (or all vertices).
# In: model = gRapHD object
#     edges = list of edges (n by 2, numeric).
#     v = vertices
# Out: vector (length(v))
################################################################################
Degree <- function(model=NULL,edges=NULL,v=NULL)
{
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
    if (class(model) != "gRapHD")
      stop("model must be of class gRapHD.")
    else
    {
      p <- model@p
      edges <- model@edges
    }
  else
    p <- max(c(max(edges),v))

  if (is.null(v))
    v <- 1:p
  else
    if (max(v)>p)
      stop(paste("v must be in [1,",p,"].",sep=""))
    
  aux <- sort(as.vector(edges))
  aux <- table(aux)
  result <- rep(0,length(v))
  names(result) <- v
  result[as.character(v)] <- aux[as.character(v)]
  result[is.na(result)] <- 0
  return(result)
}
################################################################################

################################################################################
# Return the shortest paths between v and all other vertices.
# Dijkstra's algorithm
# In: model = gRapHD object
#     edges = dataset (n by 2, numeric).
#     v = vertex
#     p = number of vertices
# Out: vector (max(edges))
################################################################################
shortPath <- function(model=NULL, edges=NULL, v=NULL, p=NULL)
{
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
    if (class(model) != "gRapHD")
      stop("model must be of class gRapHD.")
    else
    {
      p <- model@p
      edges <- model@edges
    }
  else
    if (is.null(p))
      p <- max(c(max(edges),v))

  v1 <- edges[,1]
  v2 <- edges[,2]
  storage.mode(v1) <- storage.mode(v2) <- storage.mode(p) <- "integer"
  if (!is.null(v))
  {
    if (max(v)>p)
      stop(paste("v must be in [1,",p,"].",sep=""))
    storage.mode(v) <- "integer"
    if (length(v1)==0)
    {
      result <- rep(p,p)
      result[v] <- 0
    }
    else  
      result <- .Call("shortPath",v1,v2,v,p,PACKAGE="gRapHD")
    names(result) <- 1:p
  }
  else
  {
    if (length(v1)==0)
    {
      result <- matrix(p,p,p)
      diag(result) <- 0
    }
    else  
      result <- .Call("spAll",v1,v2,p,PACKAGE="gRapHD")
  }
  result[result==p] <- Inf
  return(result)
}
################################################################################

################################################################################
# Converts the dataset to numeric format
# In: dataset = matrix (numeric) or data frame.
# Out: list = ds (matrix)
#             vertNames (vertices' names)
#             numCat (number of levels in each variable/column)
################################################################################
convData <- function(dataset)
{
  vertNames <- colnames(dataset)
  if (mode(dataset)=="numeric")
  {
    dim(dataset) <- c(NROW(dataset),NCOL(dataset))
    numCat <- rep(0,ncol(dataset))
  }
  else
    if (is.data.frame(dataset))
    {
      numCat <- rep(0,ncol(dataset))
      dsC <- NULL
      for (i in 1:ncol(dataset))
      {
        if (is.numeric(dataset[,i]))
          dsC <- cbind(dsC,dataset[,i])
        else
        {
          if (is.factor(dataset[,i]))
          {
            dsC <- cbind(dsC,as.numeric(dataset[,i]))
            numCat[i] <- nlevels(dataset[,i])
          }
          else
            if (is.logical(dataset[,i]))
            {
              dsC <- cbind(dsC,as.numeric(dataset[,i])+1)
              numCat[i] <- 2
            }
            else
              stop(paste("Column ",i," is not double, integer, factor or logic.",sep=""))
        }
      }
      dataset <- dsC
    }
    else
      stop("Input not numeric nor data frame.")
  if (is.null(vertNames))
    vertNames <- as.character(1:ncol(dataset))
  return(list(ds=dataset,numCat=numCat,vertNames=vertNames))
}
################################################################################

################################################################################
# Calculate a model's log-likelihood, AIC, and BIC
# In: model = gRapHD object
#     edges = dataset (n by 2, numeric).
#     dataset = matrix (nxp)
#     homog = T of F
# Out: number of parameters, -2*logL, BIC, AIC
################################################################################
fit <- function(model=NULL,edges=NULL,dataset,homog=NULL)
{
  if (is.null(edges) && is.null(model))
    stop("Edges or model must be provided.")
  if (!is.null(edges) && !is.null(model))
    stop("Only one (edges or model) must be provided.")
  if (is.null(edges) && (class(model)!="gRapHD"))
    stop("Model must be of class gRapHD.")

  if (is.null(edges))
  {
    edges <- model@edges
    homog <- ifelse(is.null(homog),model@homog,homog)
  }
  else
    homog <- ifelse(is.null(homog),TRUE,homog)

  aux <- convData(dataset)
  dataset <- aux$ds
  varType <- aux$numCat
  varType[varType!=0] <- 1
  numCat <- aux$numCat
  p <- length(numCat)
  n <- nrow(dataset)
  if (sum(is.na(dataset))>0)
    stop("Missing values not allowed.")
  if (!is.null(model))
  {
    if (!all.equal(numCat,model@numCat))
      stop("model@numCat doesn't agree with data.")
    if (!identical(sort(aux$vertNames),sort(model@vertNames)))
      stop("model@vertNames doesn't agree with data.")
    dataset <- dataset[,match(model@vertNames,aux$vertNames)]
  }

  if (length(unique(varType))==1) #all discrete or all continuous
    perfect <- perfSets(edges=edges,p=p,varType=varType,from=0)
  else
  {
    # star graph
    perfect <- perfSets(edges=rbind(edges,cbind(which(varType==1),p+1)),p=p+1,varType=c(varType,1),from=0)
    if (!is.list(perfect))
      stop("The model is not strongly decomposable.")
    perfect <- perfSets(edges=edges,p=p,varType=varType,from=0)
  }
  if (!is.list(perfect))
    stop("The model is not triangulated.")

  L <- 0
  numP <- 0

  whichSeq <- function(x,y)
  {
    yes <- TRUE
    i <- 0
    while ((yes) & (i<length(x)))
      yes <- x[i<-i+1]==y[i]
    return(yes)
  }

  for (i in 1:length(perfect$cliques))
  {
    x <- unlist(perfect$cliques[i])
    discrete <- x[which(varType[x]==1)]
    continuous <- setdiff(x,discrete)
    ds <- dataset[,continuous]
    dim(ds) <- c(NROW(ds),NCOL(ds))

    if ((length(continuous)>0) & (length(discrete)==0))
    {
      m <- colMeans(ds)
      covm <- var(ds)*(n-1)/n
      L <- L + sum(normDens(ds,m,covm,logx=TRUE))
      numP <- numP + length(continuous)*(length(continuous)+3)/2
    }
    else
      if ((length(continuous)==0) & (length(discrete)>0))
      {
        tc <- table(as.data.frame(dataset[,discrete]))
        tc <- tc[tc>0]
        L <- L + sum(tc*log(tc/sum(tc)))
        numP <- numP + length(tc)
      }
      else
        if (homog)
        {
          tc <- table(as.data.frame(dataset[,discrete]))
          covm <- matrix(0,NCOL(ds),NCOL(ds))
          numUsedObs <- 0
          for (j in 1:length(tc))
            if (tc[j]>length(continuous))
            {
              ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
              covm <- covm + var(ds[ind,])*(tc[j]-1)
              numUsedObs <- numUsedObs + tc[j]
            }
          covm <- covm/sum(tc)
          for (j in 1:length(tc))
            if (tc[j]>0)
            {
              ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
              m <- colMeans(matrix(ds[ind,],,length(continuous)))
              L <- L + sum(normDens(ds[ind,],m,covm,logx=TRUE)) + tc[j]*log(tc[j]/sum(tc))
              numP <- numP + 1 + length(continuous)
            }
          numP <- numP + length(continuous)*(length(continuous)+1)/2
        }
        else
        {
          tc <- table(as.data.frame(dataset[,discrete]))
          covm <- matrix(0,NCOL(ds),NCOL(ds))
          numUsedObs <- 0
          for (j in 1:length(tc))
            if (tc[j]>length(continuous))
            {
              ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
              covm <- var(ds[ind,])*(tc[j]-1)/tc[j]
              m <- colMeans(matrix(ds[ind,],,length(continuous)))
              L <- L + sum(normDens(ds[ind,],m,covm,logx=TRUE)) + tc[j]*log(tc[j]/sum(tc))
              numP <- numP + length(continuous)*(length(continuous)+3)/2 + 1
            }
        }

    x <- unlist(perfect$separators[i])
    if (length(x)>0)
    {
      discrete <- x[which(varType[x]==1)]
      continuous <- setdiff(x,discrete)
      ds <- dataset[,continuous]
      dim(ds) <- c(NROW(ds),NCOL(ds))

      if ((length(continuous)>0) & (length(discrete)==0))
      {
        m <- colMeans(ds)
        covm <- var(ds)*(n-1)/n
        L <- L - sum(normDens(ds,m,covm,logx=TRUE))
        numP <- numP - length(continuous)*(length(continuous)+3)/2
      }
      else
        if ((length(continuous)==0) & (length(discrete)>0))
        {
          tc <- table(as.data.frame(dataset[,discrete]))
          tc <- tc[tc>0]
          L <- L - sum(tc*log(tc/sum(tc)))
          numP <- numP - length(tc)
        }
        else
          if (homog)
          {
            tc <- table(as.data.frame(dataset[,discrete]))
            covm <- matrix(0,NCOL(ds),NCOL(ds))
            numUsedObs <- 0
            for (j in 1:length(tc))
              if (tc[j]>length(continuous))
              {
                ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
                covm <- covm + var(ds[ind,])*(tc[j]-1)
                numUsedObs <- numUsedObs + tc[j]
              }
            covm <- covm/sum(tc)
            for (j in 1:length(tc))
              if (tc[j]>0)
              {
                ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
                m <- colMeans(matrix(ds[ind,],,length(continuous)))
                L <- L - sum(normDens(ds[ind,],m,covm,logx=TRUE)) - tc[j]*log(tc[j]/sum(tc))
                numP <- numP - 1 - length(continuous)
              }
            numP <- numP - length(continuous)*(length(continuous)+1)/2
          }
          else
          {
            tc <- table(as.data.frame(dataset[,discrete]))
            covm <- matrix(0,NCOL(ds),NCOL(ds))
            numUsedObs <- 0
            for (j in 1:length(tc))
              if (tc[j]>length(continuous))
              {
                ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
                covm <- var(ds[ind,])*(tc[j]-1)/tc[j]
                m <- colMeans(matrix(ds[ind,],,length(continuous)))
                L <- L - sum(normDens(ds[ind,],m,covm,logx=TRUE)) - tc[j]*log(tc[j]/sum(tc))
                numP <- numP - length(continuous)*(length(continuous)+3)/2 - 1
              }
          }
    }
  }

  if (!is.element(sum(varType)/p,c(1,0)))
  {
    if (is.null(model))
      model <- new("gRapHD",edges=edges,p=p,numCat=numCat)
    model@homog <- homog
    numP <- modelDim(model)
  }
  result <- c(numP,-2*L,-2*L+2*numP,-2*L+numP*log(nrow(dataset)))
  names(result) <- c("Number of parameters","-2*Log-likelihood","AIC","BIC")
  return(result)
}
################################################################################

################################################################################
# To be used as useMethod by the summary function.
################################################################################
summary.gRapHD <- function(object,...)
{
  if (!is(object,"gRapHD"))
    stop("Not a gRapHD object.")

  cat(paste("Number of edges       = ",nrow(object@edges),"\n",sep=""))
  cat(paste("Number of vertices    = ",object@p,"\n",sep=""))
  p <- max(c(1,object@p))
  aux <- ifelse(sum(object@numCat)==0,"continuous",ifelse(sum(object@numCat!=0)/p==1,"discrete","mixed"))
  aux1 <- ifelse(aux=="mixed",paste(" and ",ifelse(object@homog,"homogeneous","heterogeneous"),sep=""),"")
  cat(paste("Model                 = ",aux,aux1," \n",sep=""))

  cat(paste("Statistic (minForest) = ",object@stat.minForest,"\n",sep=""))
  cat(paste("Statistic (stepw)     = ",object@stat.stepw,"\n",sep=""))
  cat(paste("Statistic (user def.) = ",object@stat.user,"\n",sep=""))
  cat(paste("Edges (minForest)     = ",object@minForest[1],"...",object@minForest[2],"\n",sep=""))
  cat(paste("Edges (stepw)         = ",object@stepw[1],"...",object@stepw[2],"\n",sep=""))
  cat(paste("Edges (user def.)     = ",object@userDef[1],"...",object@userDef[2],"\n",sep=""))
  cat("\n")
}
################################################################################

################################################################################
# see pages (Steffen's): 202-203; 215-216
# In: model - gRapHD object
# Out: list - discrete: (d,\empty)
#             linear: (d,\gamma^2)
#             quadratic: (d,\gamma) and (d,{\gamma,\mu})
#             quadratic2: : (d,c^2)
################################################################################
modelFormula <- function(model)
{
  if (class(model)!="gRapHD")
    stop("Model must be of class gRapHD.")

  psSub <- function(edges,p,v,varType)
  {
    n <- length(v)
    edgesV <- SubGraph(edges=edges,v=v,p=model@p)@edges
    pSeq <- perfSets(edges=edgesV,p=n,varType=varType,from=0)
    for (i in 1:length(pSeq$cliques))
      pSeq$cliques[[i]] <- v[pSeq$cliques[[i]]]
    return(pSeq$cliques)
  }

  edges <- model@edges
  varType <- model@numCat
  varType[varType!=0] <- 1
  p <- model@p

  discr <- which(varType==1)
  nDiscr <- length(discr)
  discrete <- list()
  # find all d terms
  if (nDiscr>0)
    discrete <- psSub(edges,p,discr,varType[discr])

  ind <- cbind(varType[edges[,1]],varType[edges[,2]])
  # using this cont, only those continuous that have at least one edge with
  # a discrete are considered
  cont <- sort(unique(c(edges[which((ind[,1]==0) & (ind[,2]==1)),1],
                        edges[which((ind[,2]==0) & (ind[,1]==1)),2])))
  otherCont <- sort(setdiff(which(varType==0),cont))
  # and like this, consider all continuous
  #cont <- which(varType==0)

  # find all terms (d,gamma^2)
  linear <- list()
  if (length(cont)>0)
    for (i in 1:length(cont))
    {
      aux <- list()
      pSeq <- psSub(edges,model@p,c(discr,cont[i]),c(rep(1,nDiscr),0))
      for (j in 1:length(pSeq))
        if (is.element(cont[i],pSeq[[j]]))
          aux <- c(aux,list(sort(pSeq[[j]])))
      linear <- c(linear,aux) #linear <- c(linear,list(aux))
    }

  # reduction of discrete terms
  if (length(linear)>0)
    for (i in 1:length(linear))
    {
      j <- 1
      while (j<=length(discrete))
      {
        if (length(setdiff(discrete[[j]],linear[[i]]))==0)
          discrete <- discrete[-j]
        else
          j <- j + 1
      }
    }

  if (model@homog)
  {
    quadratic <- linear
    linearX <- list()
    rm(linear)
    x <- sort(c(otherCont,cont))
    pSeq <- psSub(edges,model@p,x,varType[x])
    quadratic2 <- pSeq
  }
  else
  {
    # find all terms (d,{gamma,mu})
    quadratic <- list()
    if (length(cont)>0)
      for (i in 1:length(cont))
      {
        for (j in (i+1):length(cont))
          if ((i!=j) & (j<=length(cont)))
          {
            aux <- list()
            pSeq <- psSub(edges,model@p,c(discr,cont[c(i,j)]),
                                         c(rep(1,nDiscr),0,0))
            for (k in 1:length(pSeq))
              if ((is.element(cont[i],pSeq[[k]]))&(is.element(cont[j],pSeq[[k]])))
                aux <- c(aux,list(sort(pSeq[[k]])))#list(sort(pSeq[[k]])))
            quadratic <- c(quadratic,aux) #quadratic <- c(quadratic,list(aux))
          }
      }

    # find the indexes in linear and quadratic that have the same d
    # so, in position i of dLinead are all indexes of linear that have the same
    # d; idem for dQuadratic
    dLinear <- list()
    dQuadratic <- list()
    visitedL <- rep(FALSE,length(linear))
    visitedQ <- rep(FALSE,length(quadratic))
    v <- 1
    storage.mode(v) <- "integer"
    k <- 1
    while (sum(visitedL==0)>0)
    {
      d <- intersect(discr,linear[[v]])
      visitedL[v] <- TRUE
      dLinear[[k]] <- v
      for (i in 1:length(linear))
      {
        if (!visitedL[i])
          if ((length(setdiff(d,intersect(discr,linear[[i]])))==0) &&
              (length(setdiff(intersect(discr,linear[[i]]),d))==0))
          {
            dLinear[[k]] <- c(dLinear[[k]],i)
            visitedL[i] <- TRUE
          }
      }
      if (length(quadratic)>0)
      {
        k1 <- NULL
        for (i in 1:length(quadratic))
        {
          if (!visitedQ[i])
            if ((length(setdiff(d,intersect(discr,quadratic[[i]])))==0) &&
                (length(setdiff(intersect(discr,quadratic[[i]]),d))==0))
            {
              k1 <- c(k1,i)
              visitedQ[i] <- TRUE
            }
        }
        if (length(k1)>0)
          dQuadratic[[length(dQuadratic)+1]] <- k1
        else
          dQuadratic[[length(dQuadratic)+1]] <- integer(0)
      }
      k <- k + 1
      v <- which(visitedL==0)[1]
    }

    # remove redundant generators - for those which have a common element d for a
    # range of pairs \gamma, \mu \in c \subset \Gamma, the corresponding list
    # of generators are replaced for (d,c^2); but the way it's here, goes a bit
    # further than what is described in Steffen's book, getting some like what
    # is described in example 6.29 ("Here we have taken advantage of the fact...")
    quadratic2 <- list()
    linearX <- list()
    if (length(quadratic)>0)
    {
      for (i in 1:length(dLinear))
      {
        x <- sort(unique(c(unlist(linear[dLinear[[i]]]),
                           unlist(quadratic[dQuadratic[[i]]]))))
        pSeq <- psSub(edges,model@p,x,varType[x])
        for (j in 1:length(pSeq))
          if (length(intersect(cont,pSeq[[j]]))==1)
            linearX <- c(linearX,pSeq[j])
          else
          {
            quadratic2 <- c(quadratic2,pSeq[j])
            for (k in 1:length(dQuadratic[[i]]))
              quadratic[[dQuadratic[[i]][k]]] <- integer(0)
          }
      }
      for (i in 1:length(quadratic))
        if (length(quadratic[[i]])>0)
          quadratic2 <- c(quadratic2,quadratic[i])
    }
    else
      linearX <- linear
    rm(linear)
    # remove all empty elements of quadratic
    if (length(quadratic) > 0)
    {
      aux <- quadratic
      quadratic <- list()
      for (i in 1:length(aux))
        if (length(aux[[i]])>0)
        quadratic <- c(quadratic,aux[i])
    }
    # get the generators with only continuous variables
    x <- sort(c(otherCont,cont))
    pSeq <- psSub(edges,model@p,x,varType[x])
    for (i in 1:length(pSeq))
      if (length(intersect(pSeq[[i]],otherCont))>0)
        if (length(pSeq[[i]])==1)
          linearX <- c(linearX,pSeq[i])
        else
          quadratic2 <- c(quadratic2,pSeq[i])
  }
  return(list(discrete=discrete,
              linear=linearX,
              quadratic=quadratic,
              quadratic2=quadratic2))
}
################################################################################

################################################################################
# Model's dimension
# In: model - gRapHD object
# Out: number of parameters
################################################################################
modelDim <- function(model)
{
  if (class(model)!="gRapHD")
    stop("Model must be of class gRapHD.")

  mf <- modelFormula(model)
  exf <- c(rep(0,length(mf$discrete)),
           rep(2,length(mf$linear)),
           rep(1,length(mf$quadratic)),
           rep(2,length(mf$quadratic2)))
  mf <- c(mf$discrete,mf$linear,mf$quadratic,mf$quadratic2)
  for (i in 1:length(mf))
    mf[i] <- list(as.integer(mf[[i]]))

  numCat <- as.integer(model@numCat)
  storage.mode(exf) <- "integer"
  result <- .Call("modelDim",mf,exf,numCat,PACKAGE="gRapHD")
  return(result)
}
################################################################################

################################################################################
# To be used as useMethod by the print function.
################################################################################
print.gRapHD <- function(x,...)
{
  cat("gRapHD object\n")
  summary(x)
}
################################################################################

################################################################################
# To be used as useMethod by the show function.
################################################################################
show.gRapHD <- function(object)
{
  cat("gRapHD object\n")
  summary(object)
}
################################################################################

################################################################################
# Finds a junction tree.
# In: gRapHD object
# Out: separators - list with unique minimal separators.
#      juncTree - edges in the tree (each vertex is a clique in the list below).
#      sepSubSetOfSep - list in which each element gives all the separators
#                       which contain this respective separator.
#      indSepOrig - index of the original separator (in the MCS result) in the
#                   list above.
#      cliques - list with cliques.
################################################################################
jTree <- function(model)
{
  if (!is(model,"gRapHD"))
    stop("Object is not of gRapHD class.")
    
  v1 <- model@edges[,1]
  v2 <- model@edges[,2]
  p <- model@p
  varType <- model@numCat
  varType[varType!=0] <- 1
  storage.mode(v1) <- storage.mode(v2) <- "integer"
  storage.mode(p) <- "integer"
  storage.mode(varType) <- "integer"

  result <- .Call("juncTree",v1,v2,p,varType,PACKAGE="gRapHD")

  return(result)
}
################################################################################

################################################################################
# Calculates the clustering coefficient for the given graph
# In: model = gRapHD object
#     edges = dataset (n by 2, numeric)
#     p = number of vertices
# Output: C
################################################################################
ccoeff <- function(model=NULL,edges=NULL,p=NULL)
{
  if ((is.null(model) && is.null(edges)) || (!is.null(model) && !is.null(edges)))
    stop("Either model or edges must be provided.")
  if (is.null(edges))
  {
    if (class(model)!="gRapHD")
      stop("Model must be of class gRapHD.")
    edges <- model@edges
    p <- model@p
  }
  if (NROW(edges)==0)
    return(0)
  if (is.null(p))
    p <- max(edges)
  adjM <- adjMat(edges=edges,p=p)
  # select all nodes with edge to node i (all rows and columns connected to i)
  cc <- rep(0,p)
  for (i in 1:p)
  {
    j <- which(adjM[i,]==1)
    if (length(j) > 1)
    {
      aux <- adjM[j,j]
      Ei <- sum(aux)/2
      k <- nrow(aux)
      cc[i] <- 2*Ei/(k*(k-1))
    }
  }
  return(cc)
}
################################################################################

################################################################################
# Calculates the deviance for the inclusion of a new edge.
# In: x = vertex 1
#     y = vertex 2
#     S = separator
#     dataset = dataframe
#     homog = homogeneous (T) or not (F)
# Out: deviance and number of parameters
################################################################################
CI.test <- function(x,y,S,dataset,homog=TRUE)
{
  edge <- sort(c(x,y))
  p <- ncol(dataset) # number of variables (vertices)
  n <- nrow(dataset) # number of observations
  result <- numP <- NA

  numbPar <- function(y)
  {
    a <- sum(rowSums(y)!=0)-1
    b <- sum(colSums(y)!=0)-1
    return(a*b*(a>0))
  }
  whichSeq <- function(x,y){return(sum(abs(x-y))==0)}
  Gcliq <- function(v)
  {
    v <- sort(unique(v))
    i <- 1
    n <- length(v)
    ed <- sort(rep(i:(i+n-1),n))
    ed <- cbind(ed,rep(i:(i+n-1),n))
    ed <- matrix(ed[-unlist(lapply(1:n,function(x,n){1:x+(x-1)*n},n=n)),],,2)
    dimnames(ed) <- NULL
    ed[,1] <- v[ed[,1]]
    ed[,2] <- v[ed[,2]]
    return(ed)
  }

  dataset <- convData(dataset)
  numCat <- dataset$numCat
  dataset <- dataset$ds
  if (sum(is.na(dataset))>0)
    stop("Missing values not allowed.")

  # this is the clique that will result from the inclusion of this edge
  clique <- c(edge,S)

  # the 2 cliques thal will be merged
  clique1 <- c(edge[1],S)
  clique2 <- c(edge[2],S)

  if (sum(numCat[clique]) == 0) # continuous
  {
    CM <- cov(dataset[,clique],use="pairwise.complete.obs")*(n-1)
    if (sum(is.na(CM))) # meaning that the new edge does not have more than 1 observation (not NA)
      result <- NA
    else
    {
      a <- match(edge,clique)
      d <- log(det(CM)) + ifelse(length(S)==0,0,log(det(matrix(CM[-a,-a],length(clique)-2))))
      d <- d - (log(det(matrix(CM[-a[1],-a[1]],length(clique)-1))) +
                log(det(matrix(CM[-a[2],-a[2]],length(clique)-1))))
      result <- nrow(dataset)*d
      numP <- 1
    }
  }
  else
    if (sum(numCat[clique]!=0)==length(clique)) # discrete
    {
      t12 <- table(as.data.frame(dataset[,clique]))
      t1 <- margin.table(t12,match(clique1,clique))
      t2 <- margin.table(t12,match(clique2,clique))
      tS <- if(length(S)==0) 0 else margin.table(t12,match(S,clique))
      if (length(S)==0)
        numP <- prod(numCat[clique]-1)#length(t12)-1
      else
        numP <- sum(apply(t12,MARGIN=(1:length(clique))[-(1:2)],numbPar))
      numObsUsed <- sum(t12)
      t12 <- t12[t12!=0]
      t1 <- t1[t1!=0]
      t2 <- t2[t2!=0]
      tS <- tS[tS!=0]

      if (numObsUsed>0)
      {
        v12 <- sum(t12*log(t12/numObsUsed))
        v1  <- sum(t1 *log(t1 /numObsUsed))
        v2  <- sum(t2 *log(t2 /numObsUsed))
        vS  <- ifelse(length(S)==0,0,sum(tS *log(tS /numObsUsed)))
        result <- -2*(v12 + vS - v1 - v2)
      }
    }
    else # mixed
    {
      discrete   <- clique[numCat[clique]!=0] #intersect(clique,which(numCat!=0))
      continuous <- clique[numCat[clique]==0] #setdiff(clique,discrete)
      tab <- table(as.data.frame(dataset[,discrete]))
      ssd <- diag(0,length(continuous))
      if (sum(numCat[edge]==0)==2) # both are continuous
        if (homog)
        {
          for (j in 1:length(tab))
            if (tab[j]>0)
            {
              ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
              if (sum(ind)>1) #otherwise there is only one observation, and the contribution to the variance is zero
                ssd <- ssd + cov(dataset[ind,continuous],use="pairwise.complete.obs")*(tab[j]-1)
            }
          ind <- match(edge,continuous)
          result <- n*(log(det(ssd)) +
                    ifelse(length(continuous)>2,log(det(matrix(ssd[-ind,-ind],length(continuous)-2))),0) -
                    log(det(matrix(ssd[-ind[1],-ind[1]],length(continuous)-1))) -
                    log(det(matrix(ssd[-ind[2],-ind[2]],length(continuous)-1))))
          numP <- 1
        }
        else
        {
          result <- 0
          ind1 <- match(edge,continuous)
          for (j in 1:length(tab))
            if (tab[j]>0)
            {
              ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
              if (sum(ind)>1) #otherwise there is only one observation, and the contribution to the variance is zero
                ssd <- cov(dataset[ind,continuous],use="pairwise.complete.obs")*(tab[j]-1)
              result <- result + tab[j]*(log(det(ssd)) +
                        ifelse(length(continuous)>2,log(det(matrix(ssd[-ind1,-ind1],length(continuous)-2))),0) -
                        log(det(matrix(ssd[-ind1[1],-ind1[1]],length(continuous)-1))) -
                        log(det(matrix(ssd[-ind1[2],-ind1[2]],length(continuous)-1))))
            }
          numP <- sum(tab>(length(continuous)+1))
        }
      else
      {
        vDiscr <- edge[numCat[edge]!=0]
        vCont  <- setdiff(edge,vDiscr)
        discrMarg <- setdiff(discrete,vDiscr)
        tabMarg <- margin.table(tab,match(discrMarg,discrete))
        ssdMarg <- diag(0,length(continuous))
        if (homog)
        {
          for (j in 1:length(tab))
          {
            if (tab[j]>0)
            {
              ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
              if (sum(ind)>1) #otherwise there is only one observation, and the contribution to the variance is zero
                ssd <- ssd + var(dataset[ind,continuous],use="pairwise.complete.obs")*(tab[j]-1)
            }
            if (j<=length(tabMarg))
            {
              if (tabMarg[j]>0)
              {
                if (length(discrMarg)>0)
                  ind <- apply(matrix(dataset[,discrMarg],,length(discrMarg)),1,whichSeq,y=seqLevels(j,numCat[discrMarg]))
                else
                  ind <- rep(TRUE,tabMarg[j])
                if (sum(ind)>1) #otherwise there is only one observation, and the contribution to the variance is zero
                  ssdMarg <- ssdMarg + var(dataset[ind,continuous],use="pairwise.complete.obs")*(tabMarg[j]-1)
              }
            }
          }
          ind <- match(vCont,continuous)
          result <- n*(log(det(ssd)) +
                    ifelse(length(continuous)>1,log(det(matrix(ssdMarg[-ind,-ind],length(continuous)-1))),0) -
                    ifelse(length(continuous)>1,log(det(matrix(ssd[-ind,-ind],length(continuous)-1))),0) -
                    log(det(ssdMarg)))
          if (length(discrete)==1) #means that the separator is continuous
            numP <- numCat[vDiscr]-1
          else
          {
            aux <- match(vDiscr,discrete)
            # dimension 1 in aux1 is the discrete variable in the tested edge
            aux1 <- table(as.data.frame(dataset[,c(vDiscr,discrete[-aux])]))
            aux.tab <- margin.table(aux1,2:length(discrete))
            aux.ind <- which(aux.tab>0,arr.ind=TRUE)
            numP <- 0
            for (i in 1:nrow(aux.ind))
            {
              aux2 <- (1:numCat[vDiscr]) + sum((aux.ind[i,]-1)*cumprod(numCat[discrete[-aux]]))
              numP <- numP + (sum(aux1[aux2]>0) - 1)
            }
            rm(aux,aux1,aux2,aux.tab,aux.ind)
          }
        }
        else
        {
          result <- 0
          ind1 <- match(vCont,continuous)
          for (j in 1:length(tab))
          {
            if (tab[j]>0)
            {
              ind <- apply(matrix(dataset[,discrete],,length(discrete)),1,whichSeq,y=seqLevels(j,numCat[discrete]))
              ssd <- var(dataset[ind,continuous],use="pairwise.complete.obs")*(tab[j]-1)
              dim(ssd) <- rep(length(continuous),2)
            }
            if (j<=length(tabMarg))
            {
              if (tabMarg[j]>0)
              {
                if (length(discrMarg)>0)
                  ind <- apply(matrix(dataset[,discrMarg],,length(discrMarg)),1,whichSeq,y=seqLevels(j,numCat[discrMarg]))
                else
                  ind <- rep(TRUE,tabMarg[j])
                ssdMarg <- var(dataset[ind,continuous],use="pairwise.complete.obs")*(tabMarg[j]-1)
                dim(ssdMarg) <- rep(length(continuous),2)
              }
            }
            ind <- match(vCont,continuous)
            result <- result +
                      -tab[j]*(log(det(ssd))-ifelse(length(continuous)>1,log(det(matrix(ssd[-ind,-ind],length(continuous)-1))),0)) +
                      +ifelse(j<=length(tabMarg),tabMarg[j],0)*(log(det(ssdMarg))-ifelse(length(continuous)>1,log(det(matrix(ssdMarg[-ind,-ind],length(continuous)-1))),0))
          }
          result <- result +
                    sum(tab*log(tab),na.rm=TRUE) -
                    sum(tabMarg*log(tabMarg),na.rm=TRUE)
          result <- -result
          if (length(discrete)==1) #means that the separator is continuous
            numP <- sum(tab>(length(continuous)+1))-1
          else
          {
            aux <- match(vDiscr,discrete)
            # dimension 1 in aux1 is the discrete variable in the tested edge
            aux1 <- table(as.data.frame(dataset[,c(vDiscr,discrete[-aux])]))
            aux.tab <- margin.table(aux1,2:length(discrete))
            aux.ind <- which(aux.tab>(length(continuous)+1),arr.ind=TRUE)
            numP <- 0
            for (i in 1:nrow(aux.ind))
            {
              aux2 <- (1:numCat[vDiscr]) + sum((aux.ind[i,]-1)*cumprod(numCat[discrete[-aux]]))
              numP <- numP + (sum(aux1[aux2]>0) - 1)
            }
            rm(aux,aux1,aux.tab,aux.ind)
          }
        }
      }
    }

  return(list(deviance=-result,numP=numP))
}
################################################################################

################################################################################
# Given a "position" in a list of combination of variables levels, return the
# respective combination.
# E.g., if numCat=c(2,3,2)
#   comb = dim1 dim2 dim3
#             1    1    1
#             2    1    1
#             1    2    1
#             2    2    1
#             1    3    1
#             2    3    1
#             1    1    2
#             2    1    2
#             1    2    2
#             2    2    2
#             1    3    2
#             2    3    2
#  if k=7 then the result is c(1,1,2)
################################################################################
seqLevels <- function(k,numCat)
{
  n <- length(numCat)
  x <- rep(0,n)
  if (n>0)
    for (i in 1:n)
    {
      y <- k%%numCat[i]
      x[i] <- ifelse(y==0,numCat[i],y)
      k <- (k-y)/numCat[i] + ifelse(y>0,1,0)
    }
  return(x)
}
################################################################################

################################################################################
# Calculates the normal (uni or multivariate) density.
# In: x = vector or matrix of values
#     meanVec = mean vector (or scalar).
#     covMat = covariance matrix (or scalar)
#     logx = if to return the log or not
# Out: density (or log) values
################################################################################
normDens <- function (x, meanVec, covMat, logx = FALSE)
{
    if (NCOL(covMat)==1)
      result <- dnorm(x,mean=meanVec,sd=sqrt(covMat),log=logx)
    else
    {
      distM <- mahalanobis(x,center=meanVec,cov=covMat)
      logD <- sum(log(eigen(covMat,symmetric=TRUE,only.values=TRUE)$values))
      logVal <- -(ncol(x)*log(2*pi)+logD+distM)/2
      if (logx)
        result <- logVal
      else
        result <- exp(logVal)
    }
  return(result)
}
################################################################################

################################################################################
# To be used as useMethod by the as function.
################################################################################
# from matrix to gRapHD class
matrix.gRapHD <- function(from)
{
  result <- new("gRapHD",edges=from)
  return(result)
}
# gRapHD -> graphNEL
gRapHD.graphNEL <- function(from)
{
  #library(graph) #require(graph)
  edgeL <- vector("list",length=from@p)
  names(edgeL) <- from@vertNames
  I <- adjMat(edges=from@edges,p=from@p)
  for (i in 1:from@p)
    edgeL[[i]] <- list(edges=(1:from@p)[I[i,]==1],weights=rep(1,sum(I[i,])))
  g1 <- new("graphNEL", nodes=as.character(from@vertNames), edgeL=edgeL, edgemode = "undirected")
  return(g1)
}
# gRapHD <- graphNEL
graphNEL.gRapHD <- function(from)
{
  v <- from@nodes
  p <- length(from@nodes)
  edges <- NULL
  for (i in 1:p)
  {
    x <- from@edgeL[[v[i]]]$edges
    if (length(x)>0)
      edges <- rbind(edges,cbind(i,x))
  }
  edges <- t(apply(edges,1,sort))
  edges <- unique(edges,MARGIN=1)

  g1 <- new("gRapHD",edges=edges,p=p,vertNames=names(from@edgeL))
  return(g1)
}
################################################################################

################################################################################
# Function for backwards selection for decomposable models based on BIC or AIC.
# Assumes G is a gRapHD object and is strongly decomposable. D is a dataframe
# containing the data. fixed.edges is a boolean vector of length G@edges, when
# TRUE the edge is not considered for removal.
# For pure models, an edge (i,j) is eligible for deletion iff S(i,j)=adj(i) cap
# adj(j) is complete in G.
# For mixed models also need that if i and j are discrete, S(i,j) is all
# discrete
# When a delete-eligible edge (i,j) is deleted,
#  (i) for all edges (u,v) in S(i,j), S(u,v) becomes incomplete
#  (ii) for all edges (i,s) and (j,s) with s in S(i,j), S(i,s) and S(j,s) change
#       and must be recalculated
#  (iii) all other edges are unaffected.
################################################################################
stepb <- function(G, dataset, fixed.edges=NULL, stat="BIC") {
 test.del <- function(i, j, G) {
   S.complete <- FALSE; result <- NA; S <- NULL
   S <- intersect(neighbours(edges=G@edges, v=i), neighbours(edges=G@edges, v=j))
   if (length(S)==0) {S.complete<-TRUE; SD.condition<-TRUE} else {
     S.complete <- sum((G@edges[,1] %in% S) & (G@edges[,2] %in% S))== (length(S)*(length(S)-1)/2)
     SD.condition <- (min(G@numCat[c(i,j)])==0) | (min(G@numCat[S])>0) }
   valid <- S.complete & SD.condition
   if (valid) result <- unlist(nCI.test(i,j,S,dataset))
   return(list(ok=valid, result=result, S=S))
 }
 oneslice.deviance <- function(t,d1,d2) {
      dim(t) <- c(d1,d2); t1 <- addmargins(t)
      cm <- t1[d1+1,1:d2]; rm <- t1[1:d1,d2+1]; N<- t1[d1+1,d2+1]
      df <- (sum(cm>0)-1)*(sum(rm>0)-1)
      G2 <- 0
      if (df>0) {fv <- (rm %o% cm)/N; G2 <- 2*sum(t*log(t/fv), na.rm=T)}
      return(c(df=df, G2=G2))
  }
  cen.cov <- function(vset, dset=NULL, dataset) {
    N <- dim(dataset)[1]
    lv <- length(vset)
    if (is.null(dset)) S <- var(dataset[,vset])*(N-1)/N else {
    sds <- split(dataset[,vset], dataset[,dset], drop=T)
    covs <- sapply(sds, var)
    if (lv>1) {dims <- sapply(sds, dim)[1,]; covs1 <- covs[,dims>1]}
       else   {dims <- sapply(sds, length); covs1 <- covs[dims>1]}
    dims1 <- dims[dims>1]
    if (length(dims1)>1) S <- (covs1 %*% (dims1-1))/N else S <- covs1*(dims1-1)/N }
    dim(S) <- c(lv,lv)
    return(S)
  }
  ld <- function(S, x) {
    lds <- 0
    if (length(x)==1) lds <- log(S[x,x])
    if (length(x)>1) lds <- log(ndet(S[x,x]))
    return(lds)
  }
  # determinant function, trapping almost singular matrices. Using condition number would be better.
  ndet <- function(S) {x <- det(S); if (x<.Machine$double.eps) x<-0; x}
nCI.test <- function(c1,c2,S=integer(0),dataset) {
   V <- union(union(c1,c2),S)
   N <- nrow(dataset)
   isf <- sapply(dataset[,V], is.factor)
   isd <- c(is.factor(dataset[,c1]), is.factor(dataset[,c2]))
   type <- 'mixed';
   if (all(isd)) type<- 'discrete' else if (all(!isd)) type <- 'cont'
   if (type=='discrete') {
     if (any (S %in% V[!isf])) stop("if c1 and c2 are discrete S must also be discrete")
     d1 <- nlevels(dataset[,c1])
     d2 <- nlevels(dataset[,c2])
     ds <- dataset[,c(c1,c2,S)]
     if (length(S)==0) {
       m <- ftable(ds, col.vars=2:1)
       ans <- oneslice.deviance(m, d1, d2)
       df <- as.integer(ans[1])
       dev <- ans[2]
     } else {
       m <- ftable(ds, col.vars=2:1, row.vars=3:(length(S)+2))
       m1 <- m[rowSums(m)>1,]
       if (class(m1)=='integer') dim(m1) <- c(length(m1)/(d1*d2), d1*d2)
       ans <- apply(m1, 1, oneslice.deviance, d1, d2)
       df <- as.integer(sum(ans[1,]))
       dev <- sum(ans[2,])}
   } else if (type=='mixed') {
       Delta <- V[isf]
       Gamma <- setdiff(V, Delta)
       q <- length(Gamma)
       X <- setdiff(Gamma, union(c1,c2))
       K <- setdiff(Delta, union(c1,c2))
       if (length(K)==0) K <- NULL
       S.ik <- cen.cov(Gamma, Delta, dataset)
       S.k  <- cen.cov(Gamma, K, dataset)
       mX <- match(X, Gamma)
       ldS.ik <- log(ndet(S.ik))
       ldS.k <- log(ndet(S.k))
       ldSik.mX <- ld(S.ik, mX)
       ldSk.mX <- ld(S.k, mX)
       dev <- - N*(ldS.ik - ldSik.mX - ldS.k + ldSk.mX)
       if (is.factor(dataset[,c2])) {c.tmp <- c1; c1 <- c2; c2 <- c.tmp}
       p1 <- length(K); p <- length(X)
       if (p1>0) {
          m <- ftable(dataset[,c(c1,K)], col.vars=1, row.vars=2:(p1+1))
         if (G@homog) {m1 <- m[rowSums(m)>0,] ; df <- sum(rowSums(m1>0)-1)} else
         {m1 <- m[rowSums(m)>p+1,];  df <- sum(rowSums(m1>p+1)-1 +(p+1)*(rowSums(m1>p)-1)) }
       } else {
         m <- ftable(dataset[,c1], col.vars=1)
         if (G@homog) {df <- rowSums(m>0)-1} else
         {df <- rowSums(m>p+1)-1 +(p+1)*(rowSums(m>p)-1) }
       }
   } else {  #type==cont
       K <- V[isf]
       Q0 <- setdiff(S, K)
       Q <- setdiff(V, K)
       if (length(K)==0) K <- NULL
       S.full  <- cen.cov(Q, K, dataset)
       red1 <- match(union(Q0,c1), Q)
       red2 <- match(union(Q0,c2), Q)
       int <- match(Q0, Q)
       ldS.full <- log(ndet(S.full))
       ldS.red1 <- ld(S.full, red1)
       ldS.red2 <- ld(S.full, red2)
       ldS.int <- ld(S.full, int)
       dev <- - N*(ldS.full - ldS.red1 - ldS.red2 + ldS.int)
       df <- 1
       if (!G@homog) {
         Sd <- c(); p <- 0;
         if (!is.null(S)) {Sd <- S[sapply(dataset[,S], is.factor)]; p <- length(setdiff(S, Sd))}
         p1 <- length(Sd)
         if (p1>0) {t <- table(dataset[,Sd]); df <- sum(t>p+1)}}
   }
   return(list(deviance=dev, df=df))
}
 # Now we start stepb as such

 if (class(G) != "gRapHD") stop("G is not of class gRapHD.")
 if (class(dataset) != "data.frame") stop("dataset is not a data.frame.")
 if (G@p != ncol(dataset)) stop("G and dataset are inconsistent.")
 if (any(G@numCat != sapply(dataset,nlevels))) stop("G and dataset are inconsistent.")
 if (stat=="AIC") penalty <- 2 else penalty <- log(nrow(dataset))
 if (!is.null(fixed.edges)) {
   if (!((class(fixed.edges) == "logical") && (length(fixed.edges) == nrow(G@edges)))) stop("fixed.edges incorrectly specified.") }
 # et is a table with available edges, whether currently delete-eligible, and if so BIC or AIC statistics
 et <- data.frame(G@edges)
 names(et) <- c('i','j')
 et$eligible <- NA
 et$BIC <- NA
 if (!is.null(fixed.edges)) et <- et[!fixed.edges,]
 for (i in 1:nrow(et)) {
   tt <- test.del(et[i,1], et[i,2], G)
   et[i,3] <- tt$ok
   if (et[i,3]) et[i,4] <- tt$result[1] - penalty*tt$result[2]
 }
 k <- which.min(et[,4])
 while (et[k,4] < 0) {
   i <- et[k,1]; j <- et[k,2]
   tt <- test.del(i, j, G)
   S <- tt$S
   # remove edge
   kk <- which((G@edges[,1]==i) & (G@edges[,2]==j))
   G@edges <- G@edges[-kk,]
   # update et as necessary
   et[k,4] <- NA
   SS <- c(i,j,S)
   update.indices <- setdiff(which((et[,1] %in% SS) & (et[,2] %in% SS)),k)
   for (kk in update.indices) {
     if ((et[kk,1] %in% S) & (et[kk,2] %in% S)) {et[kk,3] <- FALSE; et[kk,4] <- NA} else {
        tt <- test.del(et[kk,1], et[kk,2], G)
        et[kk,3] <- tt$ok
        if (et[kk,3]) et[kk,4] <- tt$result[1] - penalty*tt$result[2] else
            et[kk,4] <- NA
     }
   }
   k <- which.min(et[,4])
 }
 return(G)
}
################################################################################
