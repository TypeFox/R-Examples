source("helperFunctions.R")
graphExamples = list()

## Empty graph
L = t(matrix(
  c(0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O=O+t(O)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                       globalId = 1, genId = 1, htcId = 1)))


## Ex. 3a HTC id
L = t(matrix(
  c(0, 1, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 0, 1, 1, 0,
    0, 0, 0, 1, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O=O+t(O)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 1, htcId = 1)))


## Ex 3b HTC id
L = t(matrix(
  c(0, 1, 0, 0, 1,
    0, 0, 1, 1, 0,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 0, 1, 1, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O=O+t(O)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 1, htcId = 1)))


## Ex 3c HTC inc
L = t(matrix(
  c(0, 1, 1, 1, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 1, 1, 1, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O=O+t(O)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 1, htcId = -1)))


## Ex 3d HTC id
L = t(matrix(
  c(0, 1, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 1, 0, 1, 0,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 1, 0, 1, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O=O+t(O)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 1, htcId = 1)))

## Ex 3e HTC inc
L = t(matrix(
  c(0, 1, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 0, 0, 1, 1,
    0, 0, 0, 0, 1,
    1, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 0, 1, 0, 0,
    0, 0, 1, 1, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 0), 5, 5))
O=O+t(O)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 1, htcId = -1)))

## Ex 4a HTC nonid
L = t(matrix(
  c(0, 1, 1, 1, 1,
    0, 0, 1, 1, 0,
    0, 0, 0, 1, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 0, 0, 0, 0,
    0, 0, 1, 0, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O=O+t(O)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 0, htcId = 0)))


## Ex 4b HTC inc
L = t(matrix(
  c(0, 1, 1, 1, 1,
    0, 0, 1, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 1, 1, 1, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O=O+t(O)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 0, htcId = -1)))


## Ex 4c HTC nonid
L = t(matrix(
  c(0, 1, 1, 0, 0,
    1, 0, 0, 0, 0,
    1, 0, 0, 0, 0,
    1, 0, 0, 0, 0,
    1, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 1, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O=O+t(O)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 0, htcId = 0)))


## Ex 4d HTC inc
L = t(matrix(
  c(0, 0, 0, 1, 1,
    0, 0, 1, 0, 1,
    1, 1, 0, 1, 0,
    0, 0, 0, 0, 1,
    0, 1, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 0), 5, 5))
O=O+t(O)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 0, htcId = -1)))


## Ex 5a HTC inc
L = t(matrix(
  c(0, 1, 1, 1, 1,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 1, 1, 1, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O=O+t(O)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 0, htcId = -1)))


## Ex 5b HTC inc
L = t(matrix(
  c(0, 1, 0, 1, 1,
    0, 0, 1, 0, 1,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 1, 1, 0, 1,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O=O+t(O)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 0, htcId = -1)))


## Ex 5c HTC inc
L = t(matrix(
  c(0, 1, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 1,
    1, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O=O+t(O)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 0, htcId = -1)))

## Ex 5d HTC inc

L = t(matrix(
  c(0, 1, 0, 1, 1,
    0, 0, 1, 0, 1,
    1, 0, 0, 1, 0,
    0, 1, 0, 0, 0,
    1, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 1, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O=O+t(O)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 0, htcId = -1)))


## Ex 7a (Fig 9a)
L = t(matrix(
  c(0, 1, 1, 0, 0,
    0, 0, 1, 1, 1,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 1,
    0, 0, 0, 0, 0), 5, 5))
O = t(matrix(
  c(0, 0, 0, 1, 0,
    0, 0, 1, 0, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0), 5, 5))
O=O+t(O)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 1, htcId = -1,
                            tianId = 1)))

## Sink node marginalization example from ancestral decomposition paper
dG = graph.edgelist(matrix(c(1,2, 1,3, 1,6, 2,3, 2,4, 2,5, 2,6, 3,4, 4,5),
                           ncol=2, byrow=T))
bG = graph.edgelist(matrix(c(1,6, 1,4, 2,3, 2,5, 2,6),
                           ncol=2, byrow=T), directed=F)
L = as.matrix(get.adjacency(dG))
O = as.matrix(get.adjacency(bG))
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 1, htcId = -1,
                            tianId = -1, ancId = 1)))

## Simple 2 node unidentifiable
L = t(matrix(
  c(0, 1,
    0, 0), 2, 2))
O = t(matrix(
  c(0, 1,
    1, 0), 2, 2))
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 0, htcId = 0)))

## Instrumental variable model
L = t(matrix(
  c(0, 1, 0,
    0, 0, 1,
    0, 0, 0), 3, 3))
O = t(matrix(
  c(0, 0, 0,
    0, 0, 1,
    0, 1, 0), 3, 3))
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 1, htcId = 1)))

## Several examples where ancestor decomposition was seen to be helpful
n = 6
p = .1
p1 = .2
set.seed(176796)
O = rConnectedAdjMatrix(n, p)
L = rDirectedAdjMatrix(n, p1)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 1, htcId = -1, ancId = 1)))
set.seed(335911)
O = rConnectedAdjMatrix(n, p)
L = rDirectedAdjMatrix(n, p1)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 1, htcId = -1, ancId = 1)))

set.seed(762097)
O = rConnectedAdjMatrix(n, p)
L = rDirectedAdjMatrix(n, p1)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 1, htcId = -1, ancId = 1)))

n = 8
p = .3
p1 = .4
set.seed(501)
O = rConnectedAdjMatrix(n, p)
L = rDirectedAdjMatrix(n, p1)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 1, htcId = -1, ancId = 1)))

n = 10
p = .2
p1 = .5
set.seed(3178)
O = rConnectedAdjMatrix(n, p)
L = rDirectedAdjMatrix(n, p1)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = 1, htcId = -1, ancId = 1)))

set.seed(5536)
O = rConnectedAdjMatrix(n, p)
L = rDirectedAdjMatrix(n, p1)
graphExamples = c(graphExamples,
                  list(list(L = L, O = O,
                            globalId = 0, genId = -1, htcId = -1, ancId = -1)))
