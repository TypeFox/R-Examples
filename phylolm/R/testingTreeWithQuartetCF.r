# TICR: Tree Incongruence Checking with R
# test of ILS on a population tree. 
# Cecile Ane, September 2014 - December 2015
# input: concordance factors on many (or all) 4-taxon set
#        guide tree

# See how these functions are used in file example.r

#--------------------------------------------------------#
#       standalone functions                             #
#      to analyze just one tree                          #
#--------------------------------------------------------#

test.tree.preparation <- function(cf, guidetree){
 # change all external edges of the guide tree to 0
 Ntip <- length(guidetree$tip.label)
 external.edge = which(guidetree$edge[,2] <= Ntip)
 # if the root leads to one leaf on one side and one clade on the other, then
 # the edge from the root to this clade is also 'external': no quartets map to it.
 root.edges <- which(guidetree$edge[,1] == Ntip+1) # the root has index Ntip+1 by convention in ape.
 if (length(root.edges)==2){
  if (guidetree$edge[root.edges[1],2] <= Ntip) # the root leads to a leaf
   external.edge <- c(external.edge, root.edges[2]) # the other edge should be considered external
  if (guidetree$edge[root.edges[2],2] <= Ntip) # the root leads to a leaf
   external.edge <- c(external.edge, root.edges[1]) # the other edge should be considered external
 }
 internal.edge = (1:nrow(guidetree$edge))[-external.edge]

 guidetree$edge.length[external.edge]<-0
 # this is modifying guidetree in the current environment only, not in gloval env.
 tax2id <- 1:Ntip; names(tax2id)=guidetree$tip.label
 M = nrow(cf) # number of quartets with data
 dat.names <- as.matrix(cf[,1:4]) # -> taxon names as characters, not factors
 dat.leaves <- matrix(tax2id[dat.names],nrow=M,ncol=4)
 #cat("dat.names: \n"); print(head(dat.names))
 #cat("dat.leaves:\n"); print(head(dat.leaves))

 tab0 = c(.01,.04,.05,.90)                # expected proportions of p-values
 names(tab0)=c(".01",".05",".10","large")

 cat("determining tree traversal node post-order... ")
 nodes.postorder = rep(NA,Ntip+guidetree$Nnode)
 guidetree$node.order = rep(NA,Ntip+guidetree$Nnode)
 nextnode <- 1
 set.node.postordertraversal = function(nodeid){
  # external: nextnode, nodes.postorder, guidetree$node.order
  if (nodeid > Ntip) { # not a leaf
   ces = which(guidetree$edge[,1]==nodeid) # child edges
   for (ce in ces){
    set.node.postordertraversal(guidetree$edge[ce,2])
   }
  }
  guidetree$node.order[nodeid] <<- nextnode
  nodes.postorder[nextnode] <<- nodeid
  nextnode <<- nextnode+1
 }
 set.node.postordertraversal(Ntip+1) # ID for the root
 cat("done.\n")
 #plot(guidetree);tiplabels(adj=0);nodelabels(adj=1);edgelabels()
 #cat("nodes.postorder:\n"); print(nodes.postorder)
 #cat("guidetree$node.order:\n"); print(guidetree$node.order)

 cat("calculating matrix of descendant relationships... ")
 node2descendant = matrix(FALSE, Ntip+guidetree$Nnode, Ntip+guidetree$Nnode)
 # entry [i,j] is TRUE if j is a descendant to (or equal to) i.
 for (n in nodes.postorder){
  node2descendant[n,n] <- TRUE
  if (n > Ntip) { # not a leaf
   ces <- which(guidetree$edge[,1]==n) # child edges
   for (cn in guidetree$edge[ces,2]){ # child nodes
    node2descendant[n,] <- node2descendant[cn,] | node2descendant[n,]
   }
  }
 }
 cat("done.\n")

 get.mrca <- function(nodes){
   i1 <- max(guidetree$node.order[nodes]) # first candidate node (postorder index)
   for (n in nodes.postorder[i1:(Ntip+guidetree$Nnode)])
     if (all(node2descendant[n,nodes])){ return(n) }
 }
 cat("calculating matrix of edges spanned by each quartet... ")
 dominant.table  = c(1:3,1:3) # ie tree displays 12|34 if dominant = 1
                              #                  13|24               2
                              #                  14|23               3
 dominant = rep(NA,M)
 quartet2edge = matrix(FALSE,M,length(guidetree$edge.length)) # big matrix!

 for (i in 1:M){
  x <- dat.leaves[i,]
  ancestor <- sapply(list(x[1:2],x[c(1,3)],x[c(1,4)],x[3:4],x[c(2,4)],x[2:3]),get.mrca)
  ind = which.min(guidetree$node.order[ancestor])
  a1 = ancestor[ind]
  dominant[i] = dominant.table[ind]
  mrca = ancestor[which.max(guidetree$node.order[ancestor])]
  ancestor[ind] = Ntip+1 # replacing by the root to find 2nd best, next
  a2 = ancestor[which.min(guidetree$node.order[ancestor])]
  mynode = a1
  while (mynode!=a2 & mynode!=mrca){
   parentE = which(guidetree$edge[,2]==mynode)
   quartet2edge[i,parentE]=TRUE
   mynode = guidetree$edge[parentE,1]
  }
  if (mynode == mrca){ # the mrca of all 4 taxa is on path between both ancestors
   mynode = a2
   while (mynode!=mrca){
    parentE = which(guidetree$edge[,2]==mynode)
    quartet2edge[i,parentE]=TRUE
    mynode = guidetree$edge[parentE,1]
   }
  }
 } # i=4288; dat.names[i,]; quartet2edge[i,]
 colnames(quartet2edge)=as.character(1:length(guidetree$edge.length))
 cat("done.\n")

 total.ew = rep(1,M) %*% quartet2edge[,internal.edge] # number quartet per edge
 qw =  1/quartet2edge %*% rep(1,dim(quartet2edge)[2])
 # quartet weight: 1/Ne for each edge, where Ne= # edges for that quartet
 qw.ew = t(qw) %*% quartet2edge[,internal.edge]
 score.edge = function(quartetID){
  ues = rep(1,length(quartetID)) %*% quartet2edge[quartetID,internal.edge] / total.ew
  wes = qw[quartetID] %*% quartet2edge[quartetID,internal.edge] / qw.ew
  return(rbind(ues,wes))
 }
 return(list(quartet2edge = quartet2edge,
             dominant     = dominant))
}

.plot.species.tree <- function(guidetree,edge.keep){
  sp.tree = guidetree # best if external edge lengths are 0 already
  if (length(edge.keep)==0){
   warning("Complete panmixia: the tree looks like a single tip. Won't plot this.")
   return(NULL)
  }
  sp.tree$edge.length[-edge.keep] <- 0
  plot(sp.tree)
  longedge <- which(sp.tree$edge.length > 2)
  if (length(longedge)>0)
    edgelabels(round(sp.tree$edge.length[longedge],2), longedge, frame="n", adj=c(.5,0))
  #add.scale.bar(col="red",lcol="red")
}

test.one.species.tree <- function(cf,guidetree,prep,edge.keep,plot=TRUE,shape.correction = TRUE){
  # prep should be the result of this, which takes a little while
  # so it's best to do it once and re-use it multiple times later:
  # prep = test.tree.preparation(cf,guidetree)

  if (shape.correction)
    add2shape.outlier <- add2shape.alpha <-"adaptive" else
    add2shape.outlier <- add2shape.alpha <- 0
  # add2shape: to make sure all shapes are >=1 for the Dirichlet,
  #            either to get p-values (.outlier) or to estimate alpha (.alpha).
  #            avoids Dirichlet density going to infinity at 0 or 1.
  # with add2shape=0: problems from 4-taxon sets with expected or observed CFs near 1,0,0.
  # add2shape = "adaptive" or: 0 or 1 to add 0/1 irrespective of the expected CFs.
  if (add2shape.outlier != "adaptive"){
    if (add2shape.outlier<0 | add2shape.outlier>1)
      stop("add2shape.outlier must be 'adaptive' or a numeric value between 0 and 1.")
  }
  if (add2shape.alpha != "adaptive"){
    if (add2shape.alpha<0 | add2shape.alpha>1)
      stop("add2shape.alpha must be 'adaptive' or a numeric value between 0 and 1.")
  }

  # Branches with 0 discordance cause infinite numbers and NaN likelihoods,
  # so we need to avoid that.
  maxBranchLength <- 30 # arbitrary, to replace infinite coalescent units

  # change external edge lengths to 0. Useful for plot, and to check long edges.
  guidetree$edge.length[guidetree$edge[,2]<=length(guidetree$tip.label)] <- 0
  # if the root leads to 1 leaf and 1 clade, the internal edge that should be considered
  # external is not changed: length kept as is. But not a problem for the test.

  M = nrow(cf) # number of quartets with data
  cf.obs = as.matrix(cf[,5:7])

  # if any observed CF == 0.0, change it to some very small value, to take log later
  maxBranchLength <- max(maxBranchLength, max(guidetree$edge.length))
  minObsCF <- exp(-maxBranchLength)/3 # minimum expected minor CF
  minObsCF <- min(minObsCF, min(cf.obs[cf.obs>0]))
  cf.obs[cf.obs==0] <- minObsCF
  # not re-normalizing each row to sum up to 1, but minObsCF should be < 3.2e-14
  # so the sum of each row should be almost equal to 1 up to rounding error.

  logCFsum = sum(log(cf.obs))
  tab0 = c(.01,.04,.05,.90)                # expected proportions of p-values
  names(tab0)=c(".01",".05",".10","large") # tab0*M:  .01    .05     .10   large
                                           #       274.05 1096.2 1370.25 24664.5

  #----------- get expected concordance factors -----------------------#
  if (length(edge.keep)==0) { tk = numeric(M) } else {
    tk = prep$quartet2edge[,edge.keep] %*% cbind(guidetree$edge.length[edge.keep])
    maxtk <- max(tk)
    #if (maxtk>3)
    #  warning(paste0("Some quartets display a long internal branch on the species tree",
    #    "\n  (largest: ",round(maxtk,3),") and the current implementation of TICR",
    #    "\n  tends to falsely call these quartets outliers"))
  }
  cf.exp = matrix(exp(-tk)/3, M,3) # minor CFs only, so far
  cf.exp[1:M + (prep$dominant-1)*M] = 1-exp(-tk)*2/3

  #--------------------- estimate alpha -------------------------#
  logCFtilde = sum(log(cf.obs)*cf.exp) # true mean would have this /M
  tmp = rle(sort(tk))
  tk.unique = tmp$values
  nk = tmp$lengths # number of quartets with internal edge of length tk.
  pk.minor = exp(-tk.unique)/3
  pk.domin = 1-2*pk.minor
  logCFsum.bytk = as.vector(tapply(log(cf.obs[,1])+log(cf.obs[,2])+log(cf.obs[,3]), tk, sum))
  # we should have sum(logCFsum.bytk) = logCFsum .
  # f1 = -logPL with term logCFbar ommitted ---but only with add2shape.alpha = 0.
  # f1 = function(alpha){-lgamma(alpha)*M - alpha * logCFtilde +
  #                        sum(nk*lgamma(alpha*pk.domin)) +
  #                      2*sum(nk*lgamma(alpha*pk.minor))   }
  # f2 = derivative of f1 = d(f1)/d(alpha)
  # f2 = function(alpha){-digamma(alpha)*M - logCFtilde +
  #                        sum(nk*pk.domin*digamma(alpha*pk.domin)) +
  #                      2*sum(nk*pk.minor*digamma(alpha*pk.minor)) }
  # now need to solve f2=0, or minimize f1
  # sa = (2/9)*(3*M)/sum((cf.obs-cf.exp)^2) # starting alpha value
  # optim(sa,f1,gr=f2, method="L-BFGS-B",lower=0,upper=+Inf)
  # alpha=res$par;  val=res$value
  # More general definition of f1 = -logPL below, allowing for option 'add2shape'.
  f1 = function(alpha){
    if (add2shape.alpha == "adaptive"){ # shape_j = expCF_j*alpha + b >=1 for all j=1,2,3.
      # adaptive: when b depends on the 4-taxon sets.
      b <- pmax(1 - alpha*pk.minor, rep(0,length(pk.minor))) # b of length < M.
    } else {b <- rep(add2shape.alpha,length(pk.minor))}
    -sum(nk*lgamma(alpha + 3*b)) - alpha * logCFtilde +
     sum(nk*lgamma(alpha*pk.domin+b)) +
   2*sum(nk*lgamma(alpha*pk.minor+b)) + sum((1-b)*logCFsum.bytk)
  }
  alpha.lower.bound <- 0
  #if (length(edge.keep)==0 & add2shape.alpha == "adaptive"){ alpha.lower.bound <- 3 }
  res=optimize(f1,interval=c(alpha.lower.bound,10^5))
  alpha=res$minimum;
  # minus.pll = res$objective + logCFsum # for previous def of f1 when add2shape.alpha=0
  minus.pll = res$objective       # - pseudo-log-likelihood
  if (plot){
    if (length(edge.keep)>0){
      layout(matrix(c(1:5,5),2,3,byrow=T))
      plot(cf.exp,cf.obs,ylab="Observed CFs",xlab="Expected CFs");
      abline(a=0,b=1,col="red")
    } else {
      layout(matrix(c(1:4),2,2,byrow=T))
      hist(cf.obs,xlab="Observed CFs", breaks=50, col="tan",xlim=0:1, main="",
           ylab="Number of quartets")
      mtext("Expected CFs=1/3",side=3,at=1/3,adj=0, line=0.1)
      abline(v=1/3)
    }
    aa = alpha.lower.bound + seq(.01,3*(alpha-alpha.lower.bound),length.out=500)
    plot(aa,sapply(aa,f1),type="l",xlab="alpha",ylab="neg. pseudo log-likelihood");
    points(alpha,minus.pll,col="red",pch=16)
  }
  #--------------------- get a p-value for each quartet -------------------#
  p.exp = pmax(cf.exp[,1], cf.exp[,2],cf.exp[,3])
  ind.tree = which(tk>0)
  ind.panmixia = which(tk==0)
  p.obs = pmax(cf.obs[,1], cf.obs[,2],cf.obs[,3]) # max for panmixia test
  # cf.obs[i+(j-1)*M] is same as cf.obs[i,j]
  p.obs[ind.tree] = cf.obs[ind.tree+(prep$dominant[ind.tree]-1)*M]

  if (add2shape.outlier == "adaptive"){ # to get expCF*alpha + shapeAdd >=1 for all quartets.
    shapeAdd <- pmax(1 - (1-p.exp)*alpha/2, rep(0,M)) # minor CF = (1-p.exp)/2
  } else {shapeAdd <- rep(add2shape.outlier,M)}

  pval = rep(NA,M) # p-values
  pval[ind.panmixia] = 3*pbeta(p.obs[ind.panmixia],lower.tail=F,
                               shape1=alpha/3+shapeAdd[ind.panmixia],
                               shape2=(alpha/3+shapeAdd[ind.panmixia])*2)
  dev = abs(p.exp[ind.tree]-p.obs[ind.tree]) # deviations
  pe = p.exp[ind.tree] # expected proportions
  pval[ind.tree] =
    pbeta(pe+dev,shape1=alpha*pe+shapeAdd[ind.tree],shape2=alpha*(1-pe)+2*shapeAdd[ind.tree],lower.tail=F) +
    pbeta(pe-dev,shape1=alpha*pe+shapeAdd[ind.tree],shape2=alpha*(1-pe)+2*shapeAdd[ind.tree],lower.tail=T)
  pval[pval>1] = 1
  if (plot) {
    hist(pval, breaks=100, xlim=c(0,.5),col="tan", main="",ylab="Number of 4-taxon sets",
         xlab="p-value (showing those < 0.5 only)")
    if (length(edge.keep)>0) {
      cf2plot <- cf.exp[1:M + (prep$dominant-1)*M]
    } else {
      cf2plot <- cf.obs[1:M + (prep$dominant-1)*M] }
    plot(pval~jitter(cf2plot,amount=0.005),cex=0.5,ylim=c(0,0.01),
         xlab=paste(ifelse(length(edge.keep)>0,"Expected","Observed"),
                    "CF of dominant quartet\n(jittered)"),
         ylab="p-value of 4-taxon set")
    mtext("(zooming to p-values<0.01 only)",side=3,line=0.01,cex=0.7)
    if (length(edge.keep) > 0)
      .plot.species.tree(guidetree,edge.keep)
  }

  pcat = cut(pval, breaks=c(0,.01,.05,.10,1),
                   labels=c(".01",".05",".10","large"),include.lowest=T)
  tab = table(pcat)
  res = chisq.test(tab, p=tab0) # p-value from X-squared and df=3
  if (res$p.value<0.05) {
    chisq.message <- "There is an excess of outlier quartets (with outlier p-value<=0.01).\nThis pattern could indicate an inadequate population tree or reticulate evolution.\n"
    if (tab[".01"] < M*tab0[".01"])
      chisq.message <- "The chi-square test is significant at level .05,\nbut there is a deficit of outlier quartets (with outlier p-value<=0.01).\nThis pattern does not have a simple evolutionary explanation.\n"
  } else {
    chisq.message <- "The chi-square test is not significant:\nthe population tree fits the quartet concordance factors adequately\n"}
  cat(chisq.message)
  # return tk as well? can be obtained from cf.exp
  return(list(alpha=alpha,minus.pll=minus.pll,X2=as.numeric(res$statistic),chisq.pval=res$p.value,
              chisq.conclusion=chisq.message,
              outlier.table=rbind(observed=tab,expected=tab0*M),
              outlier.pvalues=pval,cf.exp=cf.exp))
}



#--------------------------------------------------------#
#      Search through space of species trees             #
#      Many functions redefined to avoid passing/copying #
#      and garbage-collecting the very large CF data     #
#--------------------------------------------------------#

stepwise.test.tree = function(cf, guidetree, search="both", method="PLL", kbest=5,
                              maxiter=100, startT="panmixia",shape.correction=TRUE){
# search: "heuristic" (method etc. ignored) or
#         "both" directions (method etc. are used)
# kbest: lower value for faster, less thorough search.
# startT: starting partial tree for stepwise search. One of:
#         "panmixia", "fulltree", or numeric vector of edge numbers.

 if (shape.correction)
   add2shape.outlier <- add2shape.alpha <-"adaptive" else
   add2shape.outlier <- add2shape.alpha <- 0
 # add2shape: to make sure all shapes are >=1 for the Dirichlet,
 #            either to get p-values (.outlier) or to estimate alpha (.alpha).
 #            avoids Dirichlet density going to infinity at 0 or 1.
 # with add2shape=0: problems from 4-taxon sets with expected or observed CFs near 1,0,0.
 # add2shape = "adaptive" or: 0 or 1 to add 0/1 irrespective of the expected CFs.
 if (add2shape.outlier != "adaptive"){
 	if (add2shape.outlier<0 | add2shape.outlier>1)
 	  stop("add2shape.outlier must be 'adaptive' or a numeric value between 0 and 1.")
 }
 if (add2shape.alpha != "adaptive"){
  if (add2shape.alpha<0 | add2shape.alpha>1)
    stop("add2shape.alpha must be 'adaptive' or a numeric value between 0 and 1.")
 }

 # Branches with 0 discordance cause infinite numbers and NaN likelihoods,
 # so we need to avoid that.
 maxBranchLength <- 30 # arbitrary, to replace infinite coalescent units

 # change all external edges of the guide tree to 0
 Ntip <- length(guidetree$tip.label)
 external.edge = which(guidetree$edge[,2] <= Ntip)
 guidetree$edge.length[external.edge] <- 0
 # modifying guidetree in the current environment only, not in gloval env.
 root.edges <- which(guidetree$edge[,1] == Ntip+1)
 if (length(root.edges)==2){
  if (guidetree$edge[root.edges[1],2] <= Ntip)
   external.edge <- c(external.edge, root.edges[2])
  if (guidetree$edge[root.edges[2],2] <= Ntip)
   external.edge <- c(external.edge, root.edges[1])
 }
 internal.edge = (1:nrow(guidetree$edge))[-external.edge]
 tax2id = 1:Ntip; names(tax2id)=guidetree$tip.label
 M = nrow(cf) # 27405: number of quartets with data
 cf.obs = as.matrix(cf[,5:7])

 # if any observed CF == 0.0, change it to some very small value.
 maxBranchLength <- max(maxBranchLength, max(guidetree$edge.length))
 minObsCF <- exp(-maxBranchLength)/3 # minimum expected minor CF
 minObsCF <- min(minObsCF, min(cf.obs[cf.obs>0]))
 cf.obs[cf.obs==0] <- minObsCF
 # not re-normalizing each row to sum up to 1, but minObsCF should be < 3.2e-14
 # so the sum of each row should be almost equal to 1 up to rounding error. 

 dat.names = as.matrix(cf[,1:4]) # -> taxon names as characters, not factors
 dat.leaves = matrix(tax2id[dat.names],nrow=M,ncol=4)

 tab0 = c(.01,.04,.05,.90)                # expected proportions of p-values
 names(tab0)=c(".01",".05",".10","large") # tab0*M:  .01    .05     .10   large
                                          #       274.05 1096.2 1370.25 24664.5 
 cat("determining tree traversal post-order... ")
 nodes.postorder = rep(NA,Ntip+guidetree$Nnode)
 guidetree$node.order = rep(NA,Ntip+guidetree$Nnode)
 nextnode = 1
 set.node.postordertraversal = function(nodeid){
  # external: nextnode, nodes.postorder, guidetree$node.order
  if (nodeid > Ntip) { # not a leaf
   ces = which(guidetree$edge[,1]==nodeid) # child edges
   for (ce in ces){
    set.node.postordertraversal(guidetree$edge[ce,2])
   }
  }
  guidetree$node.order[nodeid] <<- nextnode
  nodes.postorder[nextnode] <<- nodeid
  nextnode <<- nextnode+1
 }
 set.node.postordertraversal(Ntip+1) # ID for the root
 #plot(guidetree);tiplabels(adj=0);nodelabels(adj=1);edgelabels()
 cat("done.\n")

 cat("calculating matrix of descendant relationships... ")
 node2descendant = matrix(FALSE, Ntip+guidetree$Nnode, Ntip+guidetree$Nnode)
 # entry [i,j] is TRUE if j is a descendant of i.
 for (n in nodes.postorder){
  node2descendant[n,n] = T
  if (n > Ntip) { # not a leaf
   ces = which(guidetree$edge[,1]==n) # child edges     
   for (cn in guidetree$edge[ces,2]){ # child nodes
    node2descendant[n,] = node2descendant[cn,] | node2descendant[n,]
   }
  }
 }
 cat("done.\n")

 get.mrca = function(nodes){
  n = max(guidetree$node.order[nodes]) # first candidate node
  for (n in nodes.postorder){
   if (all(node2descendant[n,nodes])){ return(n) }
  }
 }

 #guidetree.dist = round(cophenetic(guidetree),digits=12)
 #get.dominant = function(x){
 # quartet.dist = guidetree.dist[x,x]
 # threedist = c(quartet.dist[1,2]+quartet.dist[3,4],quartet.dist[1,3]+quartet.dist[2,4],
 #               quartet.dist[1,4]+quartet.dist[2,3])
 # return(order(threedist)[1]) #  index for which treedist is min
 #} #dominant = apply(dat.names,1,get.dominant)

 cat("calculating matrix of edges spanned by each quartet... ")
 dominant.table  = c(1:3,1:3) # ie tree displays 12|34 if dominant = 1
                              #                  13|24               2
                              #                  14|23               3
 dominant = rep(NA,M)
 quartet2edge = matrix(FALSE,M,length(guidetree$edge.length)) # big matrix!

 for (i in 1:M){
  x = dat.leaves[i,]
  ancestor = sapply(list(x[1:2],x[c(1,3)],x[c(1,4)],x[3:4],x[c(2,4)],x[2:3]),get.mrca)
  ind = which.min(guidetree$node.order[ancestor])
  a1 = ancestor[ind]
  dominant[i] = dominant.table[ind]
  mrca = ancestor[which.max(guidetree$node.order[ancestor])]
  ancestor[ind] = Ntip+1 # replacing by the root to find 2nd best, next
  a2 = ancestor[which.min(guidetree$node.order[ancestor])]
  mynode = a1
  while (mynode!=a2 & mynode!=mrca){
   parentE = which(guidetree$edge[,2]==mynode)
   quartet2edge[i,parentE]=TRUE
   mynode = guidetree$edge[parentE,1]
  }
  if (mynode == mrca){ # the mrca of all 4 taxa is on path between both ancestors
   mynode = a2
   while (mynode!=mrca){
    parentE = which(guidetree$edge[,2]==mynode)
    quartet2edge[i,parentE]=TRUE
    mynode = guidetree$edge[parentE,1]
   }
  }
 } # i=4288; dat.names[i,]; quartet2edge[i,]
 colnames(quartet2edge)=as.character(1:length(guidetree$edge.length))
 cat("done.\n")

 total.ew = rep(1,M) %*% quartet2edge[,internal.edge] # number quartet per edge
 qw =  1/quartet2edge %*% rep(1,dim(quartet2edge)[2])
    # quartet weight: 1/Ne for each edge, where Ne= # edges for that quartet
 qw.ew = t(qw) %*% quartet2edge[,internal.edge]
 score.edge = function(quartetID){
 # !Warning! external variables: cf, guidetree
 # gives score to each branch in guidetree: number/weight of 4-taxon sets 
 #       in 'quartetID' whose internal edge includes that branch
  ues = rep(1,length(quartetID)) %*% quartet2edge[quartetID,internal.edge,drop=F] / total.ew
  wes = qw[quartetID,drop=F] %*% quartet2edge[quartetID,internal.edge,drop=F] / qw.ew
  return(rbind(ues,wes))
 }

 test.species.tree = function(edge.keep,plot=FALSE,fullOutput=FALSE){
  # !Warning! external variables:  guidetree, tab0, M, cf.obs, dominant
  #----------- get expected concordance factors -----------------------#
  if (length(edge.keep)==0) { tk = numeric(M) } else { 
   tk = quartet2edge[,edge.keep] %*% cbind(guidetree$edge.length[edge.keep]) }
  cf.exp = matrix(exp(-tk)/3, M,3) # minor CFs only, so far
  # cf.exp[i,j] also = cf.exp[i+(j-1)*M] 
  cf.exp[1:M + (dominant-1)*M] = 1-exp(-tk)*2/3
  # expected CFs: q12.34, q13.24, q14.23
  # dat.names[1336,]; cf.exp[1336,]; tk[1336] # A_Lyr Bik_1 Stw_0 Ting_1, 0.484 0.269 0.247
  #--------------------- estimate alpha -------------------------#
  logCFtilde = sum(log(cf.obs)*cf.exp)
  tmp = rle(sort(tk))
  tk.unique = tmp$values
  nk = tmp$lengths # number of quartets with internal edge of length tk.
  pk.minor = exp(-tk.unique)/3
  pk.domin = 1-2*pk.minor
  logCFsum.bytk = as.vector(tapply(log(cf.obs[,1])+log(cf.obs[,2])+log(cf.obs[,3]), tk, sum))
  f1 = function(alpha){ # -logPL(alpha) = negative pseudo log-likelihood
    if (add2shape.alpha == "adaptive"){ # shape_j = expCF_j*alpha + b >=1 for all j=1,2,3.
      # adaptive: when b depends on the 4-taxon sets.
      b <- pmax(1 - alpha*pk.minor, rep(0,length(pk.minor))) # b of length < M.
    } else {b <- rep(add2shape.alpha,length(pk.minor))}
    -sum(nk*lgamma(alpha + 3*b)) - alpha * logCFtilde +
     sum(nk*lgamma(alpha*pk.domin+b)) +
   2*sum(nk*lgamma(alpha*pk.minor+b)) + sum((1-b)*logCFsum.bytk)
  }
  #if (length(edge.keep)==0 & add2shape.alpha == "adaptive"){ alpha.lower.bound <- 3 }
  res=optimize(f1,interval=c(0,10^5))
  alpha=res$minimum
  minus.pll = res$objective # - pseudo-log-likelihood at best alpha
  #--------------------- get a p-value for each quartet -------------------#
  p.exp = pmax(cf.exp[,1], cf.exp[,2],cf.exp[,3])
  ind.tree = which(tk>0)  
  ind.panmixia = which(tk==0)
  p.obs = pmax(cf.obs[,1], cf.obs[,2],cf.obs[,3]) # max for panmixia test
  # cf.obs[i+(j-1)*M] is same as cf.obs[i,j]
  p.obs[ind.tree] = cf.obs[ind.tree+(dominant[ind.tree]-1)*M] 
  if (add2shape.outlier == "adaptive"){ # to get expCF*alpha + shapeAdd >=1 for all quartets.
    shapeAdd <- pmax(1 - (1-p.exp)*alpha/2, rep(0,M))
    # minor CF = (1-p.exp)/2
  } else { shapeAdd <- rep(add2shape.outlier,M) }

  pval = rep(NA,M) # p-values
  pval[ind.panmixia] = 3*pbeta(p.obs[ind.panmixia],lower.tail=F,
                               shape1=alpha/3+shapeAdd[ind.panmixia],
                               shape2=(alpha/3+shapeAdd[ind.panmixia])*2)
  dev = abs(p.exp[ind.tree]-p.obs[ind.tree]) # deviations
  pe = p.exp[ind.tree] # expected proportions
  pval[ind.tree] =
    pbeta(pe+dev,shape1=alpha*pe+shapeAdd[ind.tree],shape2=alpha*(1-pe)+2*shapeAdd[ind.tree],lower.tail=F) +
    pbeta(pe-dev,shape1=alpha*pe+shapeAdd[ind.tree],shape2=alpha*(1-pe)+2*shapeAdd[ind.tree],lower.tail=T)
  pval[pval>1] = 1

  pcat = cut(pval, breaks=c(0,.01,.05,.10,1),
                   labels=c(".01",".05",".10","large"),include.lowest=T)
  tab = table(pcat)
  res = chisq.test(tab, p=tab0) # p-value from X-squared and df=3
  tmp <- list(alpha=alpha,minus.pll=minus.pll,X2=as.numeric(res$statistic),
              chisq.pval=res$p.value,p.01=as.numeric(tab[".01"]),p.05=as.numeric(tab[".05"]),
              p.10=as.numeric(tab[".10"]),p.large=as.numeric(tab["large"]),pval=pval)
  if (fullOutput){
    tmp$outlier.table <- rbind(observed=tab,expected=tab0*M)
    tmp$cf.exp <- cf.exp
  }
  if (plot){
   if (length(edge.keep)>0) {
    layout(matrix(c(1:5,5),2,3,byrow=T))
    plot(cf.exp,cf.obs,ylab="Observed CFs",xlab="Expected CFs");
    abline(a=0,b=1,col="red")
    cf2plot <- cf.exp[1:M + (dominant-1)*M]
   } else {
    layout(matrix(1:4,2,2,byrow=T))
    hist(cf.obs,xlab="Observed CFs",breaks=50,col="tan",xlim=0:1,main="",ylab="Number of quartets")
    mtext("Expected CFs=1/3",side=3,at=1/3,adj=0, line=0.1)
    abline(v=1/3)
    cf2plot <- cf.obs[1:M + (dominant-1)*M]
   }
   aa = seq(.01,3*alpha,length.out=500)
   plot(aa,sapply(aa,f1),type="l",xlab="alpha",ylab="neg. pseudo log-likelihood");
   points(alpha,minus.pll,col="red",pch=16)
   hist(pval, breaks=100, xlim=c(0,.5),col="tan", main="",ylab="Number of 4-taxon sets",
       xlab="p-value (showing those < 0.5 only)")
   plot(pval~jitter(cf2plot,amount=0.005),cex=0.5,ylim=c(0,0.01),
        xlab=paste(ifelse(length(edge.keep)>0,"Expected","Observed"),
                   "CF of dominant quartet\n(jittered)"),
        ylab="p-value of 4-taxon set")
   mtext("(zooming to p-values<0.01 only)",side=3,line=0.01,cex=0.7)
   if (length(edge.keep)>0)
     .plot.species.tree(guidetree,edge.keep)
  }
  return(tmp)
 }

 # main: heuristic search through partial trees
 if (search=="heuristic"){
  mystat = data.frame(number.edges=0:length(internal.edge))
  mystat$newedge=NA; mystat$alpha=NA; mystat$negPseudoLoglik=NA; mystat$X2=NA;
  mystat$p.01=NA; mystat$p.05=NA; mystat$p.10=NA;
  for (Nedge in mystat$number.edges){ 
   if (Nedge==0){ edge2keep = c(); # starting from panmixia 
   } else { 
    if (Nedge==1){ nextedge = as.numeric(names(which.max(res2[2,]))) }
    if (Nedge> 1){
     edge2keep.2 = which(colnames(res2) %in% as.character(edge2keep))
     nextedge = as.numeric(names(which.max(res2[2,-edge2keep.2])))
    }
    mystat$newedge[Nedge+1] = nextedge
    edge2keep = c(edge2keep, nextedge)
    cat("new edge to keep:",nextedge,"\n")
   }
   res1 = test.species.tree(edge2keep)
   mystat$negPseudoLoglik[Nedge+1] = res1$minus.pll
   mystat$alpha[Nedge+1] = res1$alpha
   mystat$X2[Nedge+1]    = res1$X2
   mystat$p.01[Nedge+1]  = res1$p.01
   mystat$p.05[Nedge+1]  = res1$p.05
   mystat$p.10[Nedge+1]  = res1$p.10
   res2 = score.edge(which(res1$pval<.01));
   # print(res2[,order(res2[2,],decreasing=T)]) # to see edge scores
  }
  return(mystat)
 }
 # main: stepwise, forward+backward search through partial trees
 if (search=="both"){
  if (all(startT=="panmixia")){ edge2keep.current = c() }         else {
    if (all(startT=="fulltree")){edge2keep.current=internal.edge} else {
      if (is.numeric(startT)){ edge2keep.current = startT }       else {
        stop('problem with bad startT. Options: "panmixia", "fulltree", or numeric vector of edges')
      }
    }
  }
  Nedge = length(edge2keep.current)
  bestRes1 = test.species.tree(edge2keep.current)
  res2 = score.edge(which(bestRes1$pval<.01))
  if (method=="PLL"){ bestCriterion = bestRes1$minus.pll } # criterion
  iter = 0
  cat("iter=",iter,"\n\tNedge=",Nedge,"\n\tedges=",edge2keep.current)
  cat("\n\talpha=",bestRes1$alpha,"\n\t-PLL =",bestRes1$minus.pll)
  cat("\n\tX2   =",bestRes1$X2,"\n\tcrit.=",bestCriterion,"\n")
  flush.console()
  while (iter<maxiter){
   iter = iter+1
   bestAction = "stop"
   bestEdge = NULL
   # forward search
   if (Nedge<length(internal.edge)){
    candidates = setdiff(internal.edge,edge2keep.current)
    myk = min(kbest, length(candidates))
    candidates = res2[2,as.character(candidates)]
    edges2add = sort.int(candidates,index.return=T,decreasing=T)$ix[1:myk]
    edges2add = as.numeric(names(candidates[edges2add]))
    #cat("\tForward: trying edges",edges2add,"\n")
    for (ed in edges2add){
     edge2keep = c(edge2keep.current, ed)
     res1 = test.species.tree(edge2keep)
     #cat("\t edge:",ed,"-PLL=",res1$minus.pll,"\n")
     if (method=="PLL"){
      if (res1$minus.pll < bestCriterion) {
       bestCriterion = res1$minus.pll
       bestAction="add"
       bestEdge = ed
       bestRes1 = res1
     }}
    }
   }
   # backward search
   if (Nedge >0){
    #myk = min(kbest,Nedge)
    #candidates = res2[2,as.character(edge2keep.current)]
    #edges2remove = sort.int(candidates,index.return=T,decreasing=T)$ix[1:myk]
    #edges2remove = as.numeric(names(candidates[edges2remove]))
    #edges2removeI = which(edge2keep.current %in% edges2remove)
    #cat("\tBackward: trying edges",edges2remove,"\n")
    #cat("\tBackward:\n")
    for (i in 1:length(edge2keep.current)){ # or in edges2removeI
     edge2keep = edge2keep.current[-i]
     res1 = test.species.tree(edge2keep)
     #cat("\t edge:",edge2keep.current[i],"-PLL=",res1$minus.pll,"\n")
     if (method=="PLL"){
      if (res1$minus.pll < bestCriterion) {
       bestCriterion = res1$minus.pll
       bestAction="remove"
       bestEdge = i # index. the edge itself is edge2keep.current[i]
       bestRes1 = res1
     }}
    }
   }
   cat("iter=",iter," action=",bestAction,
       ifelse(bestAction=="stop","",", edge="),bestEdge,"\n", sep="")
   if (bestAction=="stop"){ break }
   if (bestAction=="add"){ edge2keep.current=c(edge2keep.current,bestEdge)}
   if (bestAction=="remove"){ edge2keep.current=edge2keep.current[-bestEdge]}
   res2 = score.edge(which(bestRes1$pval<.01))
   Nedge = length(edge2keep.current)
   #{cat("\tNedge=",Nedge,"\n\tedges=",edge2keep.current)
   #cat("\n\talpha=",bestRes1$alpha,"\n\t-PLL =",bestRes1$minus.pll)
   #cat("\n\tX2   =",bestRes1$X2,"\n\tcrit.=",bestCriterion,"\n")}
  }

  bestRes1 <- test.species.tree(edge2keep.current, plot=T, fullOutput=T)
  if (bestRes1$chisq.pval<0.05) {
    chisq.message <- "There is an excess of outlier quartets (with outlier p-value<=0.01).\nThis pattern could indicate an inadequate population tree or reticulate evolution.\n"
    if (bestRes1$outlier.table["observed",".01"] < bestRes1$outlier.table["expected",".01"])
      chisq.message <- "The chi-square test is significant at level .05,\nbut there is a deficit of outlier quartets (with outlier p-value<=0.01).\nThis pattern does not have a simple evolutionary explanation.\n"
  } else {
    chisq.message <- "The chi-square test is not significant:\nthe population tree fits the quartet concordance factors adequately\n"}
  cat(chisq.message)
  return(list(Nedge=Nedge,edges=sort(edge2keep.current),
              notincluded = sort(setdiff(internal.edge,edge2keep.current)),
              alpha=bestRes1$alpha,
              negPseudoLoglik=bestRes1$minus.pll, X2=bestRes1$X2,
              chisq.pval=bestRes1$chisq.pval, chisq.conclusion=chisq.message,
              outlier.table=bestRes1$outlier.table,
              outlier.pvalues=bestRes1$pval, cf.exp=bestRes1$cf.exp))
 }
}

