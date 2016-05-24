# tests for temporal degree functions
library(tsna)
library(testthat)
library(sna)

data(windsurfers)
tDegree(windsurfers)
data(McFarland_cls33_10_16_96)
tDegree(cls33_10_16_96)

library(networkDynamicData)
data(concurrencyComparisonNets)

# alternate version using sna's degree function for comparison testing
require(sna)
tSnaDegree<-function(nd,start, end, time.interval=1,cmode=c('freeman','indegree','outdegree')){
  if(missing(start) | missing(end)){
    times <- get.change.times(nd)
    if (length(times) == 0) {
      warning("network does not appear to have any dynamic information. Using start=0 end=1")
      start = 0
      end = 0
    }
    times[times == Inf] <- NA
    times[times == -Inf] <- NA
    start = min(times, na.rm = T)
    end = max(times, na.rm = T)
  }
  
  # figure out the times where we will do evaluations
  times<-seq(from = start, to=end,by = time.interval)
  
  # figure out gmode
  if(is.directed(nd)){
    gmode<-'digraph'
  } else {
    gmode<-'graph'
  }
  
  cmode<-match.arg(cmode)
  
  # allocate a matrix to hold the degrees of each vertex at each evaluation time
  degrees<-matrix(NA,nrow = network.size(nd),ncol=length(times),dimnames = list(paste('v',seq_len(network.size(nd)),sep='')),times)
  for (i in seq_len(length(times))){
    at<-times[i]
    slice<-network.extract(nd,at=at)
    vids<-which(is.active(nd,at=at,v = seq_len(network.size(nd))))
    degreesAt<-degree(as.matrix(slice),gmode=gmode,diag=has.loops(nd),cmode=cmode)
    # deal with case of no active vertices where sapply returns list()
    if (!is.list(degreesAt)){
      degrees[vids,i]<-degreesAt
    }
    
  }
  # transpose to match new orientation and change class to ts
  degrees<- ts(t(degrees), start = start, end = times[length(times)], 
     deltat = time.interval)
  return(degrees)
}

# check class of returned object
expect_true(is.ts(tDegree(windsurfers)))

# comparison with vertex dynamics
expect_equal(tDegree(windsurfers),tSnaDegree(windsurfers),check.attributes=FALSE)

# comparison with directed network
expect_equal(tDegree(cls33_10_16_96),tSnaDegree(cls33_10_16_96),check.attributes=FALSE)

expect_equal(tDegree(cls33_10_16_96,cmode='indegree'),tSnaDegree(cls33_10_16_96,cmode='indegree'),check.attributes=FALSE)
expect_equal(tDegree(cls33_10_16_96,cmode='outdegree'),tSnaDegree(cls33_10_16_96,cmode='outdegree'),check.attributes=FALSE)

# comparison with subset of undirected network
expect_equal(tDegree(base,start=0,end=103,time.interval = 10),tSnaDegree(base,start=0,end=103,time.interval = 10),check.attributes=FALSE)

# comparison with loops
loopy<-network.initialize(5,loops = TRUE)
add.edges.active(loopy,1,1,onset=0,terminus=1)
add.edges.active(loopy,2,3,onset=0,terminus=1)
expect_equal(tDegree(loopy),tSnaDegree(loopy),check.attributes=FALSE)

# comparison with undeclared loops
loopy2<-network.initialize(5,loops = FALSE)
add.edges.active(loopy2,1,1,onset=0,terminus=1)
add.edges.active(loopy2,2,3,onset=0,terminus=1)
# this gives error because sna version doesn't know to calculate the diagonal
#expect_equal(tDegree(loopy2),tSnaDegree(loopy2),info = 'network with undeclared loops',check.attributes=FALSE)


# evaluation on a resonably sized network using start end and time interval
degs<-tDegree(base,start = 0,end=102,time.interval = 10)
expect_equal(nrow(degs),11)
expect_equal(ncol(degs),network.size(base))

# check labeling
expect_equal(colnames(degs),as.character(1:1000))

# check results close to expected value for the network
expect_equal(mean(colMeans(degs)),0.7609091)


