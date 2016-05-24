library(catnet)

numnodes <- 24
numcats <- 3
maxpars <- 3
cn <- cnRandomCatnet(numnodes, maxpars, numcats)

ps <- cnSamples(cn, 500)

## finds a node with descendants and reduce its number of categories
mpars <- cnMatParents(cn)
for(j in 1:numnodes)
  if(sum(mpars[,j]) > 0)
  break
if(j < numnodes)
  cn@categories[[j]] <- cn@categories[[j]][1:(numcats-1)]
cn <- cnSetProb(cn, ps)

ps <- cnSamples(cn, 500)

res1 <- cnSearchOrder(data=ps, perturbations = NULL,
                      maxParentSet = maxpars, parentSizes = NULL, 
                      maxComplexity = 0,
                      nodeOrder = cnOrder(cn),  
                      parentsPool = NULL, fixedParents = NULL,
                      edgeProb = NULL, 
                      echo = FALSE)
res1
anet1 <- cnFind(res1, cnComplexity(cn))
cnCompare(cn, anet1)

matliks <- matrix(rep(0.5, numnodes*numnodes), nrow=numnodes)
##matliks <- matrix(runif(numnodes*numnodes, 0.49, 0.51), nrow=numnodes)
matliks <- 0.5 + cnMatParents(cn) / 4

res2 <- cnSearchOrder(data=ps, perturbations = NULL,
                      maxParentSet = maxpars, parentSizes = NULL, 
                      maxComplexity = 0,
                      nodeOrder = cnOrder(cn),  
                      parentsPool = NULL, fixedParents = NULL,
                      edgeProb = matliks, 
                      echo = TRUE)
res2
anet2 <- cnFind(res2, cnComplexity(cn))
cnCompare(cn, anet2)
cnCompare(anet1, anet2)

##cnDot(list(cn, anet2), "cn")

sares <- cnSearchSA(data=ps, perturbations=NULL,
                    maxParentSet=maxpars, parentSizes = NULL, 
                    maxComplexity = 0,
                    parentsPool = NULL, fixedParents = NULL, edgeProb = matliks, 
                    selectMode = "BIC", 
                    tempStart = 1, tempCoolFact = 0.9, tempCheckOrders = 20, 
                    maxIter = 100, orderShuffles = -1, stopDiff = 1,
                    numThreads = 2, 
                    priorSearch = NULL,
                    echo=TRUE)
anet3 <- cnFind(sares, cnComplexity(cn))
cnCompare(cn, anet3)

sares2 <- cnSearchSA(data=ps, perturbations=NULL,
                    maxParentSet=maxpars, parentSizes = NULL, 
                    maxComplexity = 0,
                    parentsPool = NULL, fixedParents = NULL, edgeProb = matliks, 
                    selectMode = "BIC", 
                    tempStart = 1, tempCoolFact = 0.9, tempCheckOrders = 20, 
                    maxIter = 100, orderShuffles = -1, stopDiff = 1,
                    numThreads = 2, 
                    priorSearch = sares,
                    echo=TRUE)
anet4 <- cnFind(sares2, cnComplexity(cn))
cnCompare(cn, anet4)

mm <- cnSearchHist(data=ps, perturbations=NULL,  
                         maxParentSet=maxpars, parentSizes = NULL,
                         maxComplexity=0,
                         parentsPool = NULL, fixedParents = NULL,
                         selectMode = "BIC", 
                         maxIter = 60, numThreads = 2, echo=TRUE)
matliks <- 0.5 + mm/ (4*max(mm))
mm <- mm > quantile(mm, 0.9)
cnDot(mm, "mm")

mmt <- matliks > quantile(matliks, 0.8)
cnDot(mmt, "mmt")

sares3 <- cnSearchSA(data=ps, perturbations=NULL,
                    maxParentSet=maxpars, parentSizes = NULL, 
                    maxComplexity = 0,
                    parentsPool = NULL, fixedParents = NULL, edgeProb = matliks, 
                    selectMode = "BIC", 
                    tempStart = 1, tempCoolFact = 0.9, tempCheckOrders = 20, 
                    maxIter = 100, orderShuffles = -1, stopDiff = 1,
                    numThreads = 4, 
                    priorSearch = NULL,
                    echo=TRUE)
anet5 <- cnFind(sares3, cnComplexity(cn))
cnCompare(cn, anet5)

####################################################################################
## effectivelly nil network

library(catnet)
cnet <- cnNew(nodes = c("a", "b", "c"), cats = list(c("1", "2"),
     c("1", "2"), c("1", "2")),
              parents = list(NULL, c(1), c(1,
     2)),
              probs = list(c(0.5, 0.5), list(c(0.5, 0.5), c(0.5, 0.5)),
     list(list(c(0.5, 0.5), c(0.5, 0.5)), list(c(0.5, 0.5), c(0.5,
         0.5)))))
cnet
cnComplexity(cnet)
cnKLComplexity(cnet)
