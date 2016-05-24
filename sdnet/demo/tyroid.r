library(sdnet)

data(tyr1)
data(tyr2)

n1 <- ncol(tyr1$cdata)
n2 <- ncol(tyr2$cdata)
N <- nrow(tyr1$cdata)
clslevs <- range(tyr1$cls)
n11 <- sum(tyr1$cls==clslevs[1])
n12 <- sum(tyr1$cls==clslevs[2])
ncats <- 3
genenames <- rownames(tyr1$cdata)
nodeCats <- lapply(1:N, function(i) 1:ncats)
names(nodeCats) <- genenames

## soft discretization
q <- sdnet::cnDiscretize(cbind(tyr1$cdata,tyr2$cdata), numcats=ncats, mode="soft", marginal="quantile", learnset=1:n1, cover=0.95)
ptest <- q$pdata[,(n1+1):(n1+n2)]
dtest <- q$ddata[,(n1+1):(n1+n2)]
plearn <- q$pdata[,1:n1]
rm(q)

bnet <- sdnet::cnNew(genenames, cats=nodeCats, pars=vector("list",N), probs=NULL, dagonly=TRUE)

## sets P(X=k) \propto \sum_s q_k(y^s) = P(y^s|X=k)/(\sum_m P(y^s|X=m))
bnet <- sdnet::cnSetProb(bnet, plearn, nodeCats=nodeCats[bnet@nodes], softmode=TRUE)
net1 <- sdnet::cnSetProb(bnet, plearn[, tyr1$cls==clslevs[1]], nodeCats=nodeCats[bnet@nodes], softmode=TRUE)
net2 <- sdnet::cnSetProb(bnet, plearn[, tyr1$cls==clslevs[2]], nodeCats=nodeCats[bnet@nodes], softmode=TRUE)
rm(plearn)

catnet.res <- cnSearchOrder(ptest, NULL, 2, NULL, 0, 
                              nodeOrder=NULL, nodeCats=bnet@cats, 
                              echo=FALSE, dagsOnly=TRUE, softmode=TRUE, classes=NULL, clsdist=2)
  ## test H0: G == G0, df(G0) ~ nn  
catnet.res
catnet.res@loglik[1]
abs(catnet.res@loglik[1] + 50.6301) < 1e-6

  ss <- catnet.res@loglik/catnet.res@complx
  id <- which(ss == max(ss))
  bnet0 <- cnFind(catnet.res, catnet.res@complx[id])

abs(bnet0@loglik + 2.320366e-14) < 1e-6
bnet0@complx == 1682

bnet0 <- cnFindAIC(catnet.res)
bnet0@complx == 418

mres <- cnSearchOrder(dtest, NULL, 2, NULL, 0, 
                      nodeOrder=NULL, nodeCats=bnet@cats, 
                      echo=FALSE, dagsOnly=TRUE, softmode=FALSE)
mres

bnet1 <- cnFindAIC(mres)
bnet1@complx == 366

##mres <- cnSearchHist(dtest, NULL, 2, NULL, maxComplexity=680, numThreads = 4, score = "AIC", echo=TRUE)
##sum(mres) > 0

pl1 <- -sapply(1:net1@numnodes, function(j) sum(bnet@probs[[j]]*log(net1@probs[[j]]/bnet@probs[[j]])))
pl1 <- 2*(n1+n2)*(n1/n2)*pl1
ind <- order(pl1, decreasing = TRUE)
hc <- sapply(1:length(pl1), function(k) sum(pl1[ind[1:k]]) - 0.5*log(n1+n2)*2*k)
kmax <- which(hc==max(hc[!is.nan(hc) & abs(hc)<Inf]))[1]
hc <- hc[1:kmax]
if(is.na(kmax) || kmax<2) kmax <- 2
ind <- ind[1:kmax]
rm(bnet)

predict <- NULL
for(k in 1:ncol(ptest)) {
  dd <- t(matrix(ptest[,k], nrow=ncats))
  pl1 <- sapply(ind, function(j) sum(net1@probs[[j]]*log(dd[j,]/net1@probs[[j]])))
  pl2 <- sapply(ind, function(j) sum(net2@probs[[j]]*log(dd[j,]/net2@probs[[j]])))
  pdl <- pl2- pl1
  pdl <- sum(pdl[!is.nan(pdl)&abs(pdl)<Inf])
  sel <- clslevs[1]
  if(pdl >= 0)
    sel <- clslevs[2]
  predict <- c(predict, sel)
}
rm(ptest)

## compare the prediction to the true classes 
r <- cbind(tyr2$cls, predict)
colnames(r) <- c("tyr2$cls", "prediction")

acc <- 0.5*(sum(predict==clslevs[1]&tyr2$cls==clslevs[1])/sum(tyr2$cls==clslevs[1]) + sum(predict==clslevs[2]&tyr2$cls==clslevs[2])/sum(tyr2$cls==clslevs[2]))
cat("balanced accuracy of predicting (tyr1 -> tyr2) = ", acc, "\n")

