## [Currently unused completely]
get.signif.clusters <- function(tree,
                                c.test,
                                extinct)
{
  ## This is not required anymore
  ## if(length(extinct)>0)## check if the excting numbers are truly the leafs of extinct branches
  ##  stopifnot(all(unlist(lapply(tree[extinct],function(node) is.na(node$child1) && is.na(node$child2)))))
  signif.clusters <- extinct
  for(i in c.test) ## TODO w/o for()
    {
      signif.clusters <- c(signif.clusters,tree[[i]]$parent)
    }

  ## return
  unique(signif.clusters)
}

##' Only called from mssplit.hierarch.testing() :
signif.partialFtest <- function(tree,
                                c,## current cluster to test
                                x,
                                y,
                                split.out,
                                family)
{
  stopifnot((ns <- length(split.out)) >= 1)
  pvalue <- numeric(ns)
  for(i in 1:ns) {
      sel.model <- which(split.out[[i]]$sel.models)
      intersection <- intersect(tree[[c]]$vars, sel.model)
      if((hasI <- length(intersection) > 0)) {
          yt <- y[-split.out[[i]]$split ]
          xt <- x[-split.out[[i]]$split,]
      }
      pvalue[i] <-
        ## MM speed FIXME: use lm.fit() and  pf() {for p-val of F test}
        if(!hasI) {
            1
        } else if(length(intersection)==length(sel.model)) {
            ## the intersection equals the screened set
            anova(lm(yt~1),
                  lm(yt~xt[,sel.model,drop=FALSE]),## RIGHT? we should only look at what
                  test="F")$P[2]
        } else {
              Shat.withoutc <- setdiff(sel.model, tree[[c]]$vars)
              anova(lm(yt~xt[, Shat.withoutc, drop=FALSE]),
                    lm(yt~xt[, sel.model, drop=FALSE]),## RIGHT? we should only look at what
                    test="F")$P[2]
        }
    }
  ## return
  pvalue
}

find.onesidedparents <- function(c,
                                 tree,
                                 extinct,
                                 inheritance.enable=TRUE) {
  onesidedparents <- c()
  if(inheritance.enable && length(extinct)>=1) {
    ## we assume that all nodes that are part of an extinct branch are in the vector extinct
    currentnode <- c
    while(TRUE) {
      oldnode <- currentnode
      currentnode <- tree[[currentnode]]$parent
      if(is.na(currentnode)) {
        ## we reached the top
        ## no more one sided parents to be found
        break
      }
      if((tree[[currentnode]]$child1 == oldnode || tree[[currentnode]]$child1 %in% extinct) &&
         (tree[[currentnode]]$child2 == oldnode || tree[[currentnode]]$child2 %in% extinct)) {
        ## the current node is a one-sided parent!
        onesidedparents <- c(onesidedparents,currentnode)
      }
    }
  }
  return(onesidedparents)
}

getmultadj.fortree <- function(tree,
                               c,
                               split.out,
                               extinct,
                               c.test)
{## we assume here that the c has a parent that is significant
  ## this should be implicity true in the way we loop over the nodes in the R code

  ## This approach is not necessarily efficient
  one.sided.parents <- find.onesidedparents(c=c,
                                            tree=tree,
                                            extinct=extinct)## find the ancestor nodes where the other branch side is extinct

  mC <- numeric(length(split.out))
  for(i in 1:length(split.out)) {
      sel.model <- which(split.out[[i]]$sel.models)
      if(any(tree[[c]]$vars %in% sel.model)) {
          mC[i] <- length(sel.model)/length(intersect(tree[[c]]$vars,sel.model))
          ## Inheritance procedure
          for(p in one.sided.parents) {
            notextinct.child <- setdiff(c(tree[[p]]$child1,
                                          tree[[p]]$child2),
                                        extinct)
            mC[i] <- mC[i] * length(intersect(tree[[notextinct.child]]$vars,sel.model))/
              length(intersect(tree[[p]]$vars,sel.model))
          }

          ## Shaeffer improvement
          if(!is.na(tree[[c]]$parent)) {
            sibling <- setdiff(c(tree[[tree[[c]]$parent]]$child1,
                                 tree[[tree[[c]]$parent]]$child2),
                               c)
            if(sibling %in% c.test) {
              ## the sibling is a leaf and is not rejected yet
              mC[i] <- mC[i] * length(intersect(tree[[c]]$vars,sel.model))/
                length(intersect(tree[[tree[[c]]$parent]]$vars,sel.model))
            }
          }
        }else{
          mC[i] <- 1
        }
    }
  mC
}

update.extinct <- function(nowextinct,
                           tree,
                           c.test,
                           oldextinct)
{
  ## with the new extinct cluster, update all parts of the tree that are now extinct

  newextinct <- c(nowextinct,oldextinct)

  current.node <- nowextinct
  while(TRUE) {
    ## go up in the tree till you find a node with a second child that is not extinct yet
    current.node <- tree[[current.node]]$parent
    found.notextinct <- (!tree[[current.node]]$child1 %in% newextinct ||
                         !tree[[current.node]]$child2 %in% newextinct)
    if(found.notextinct)
      break## stop, this node is not extinct and therefore we are done
    ## else
    ##  print("We found extra extincts!")
    newextinct <- c(newextinct,current.node)
  }
  if(any(duplicated(newextinct)))
    message("we have some duplicates in newextinct")
  ## return
  unique(newextinct)
}

mssplit.hierarch.testing <- function(tree,
                                     hh,
                                     x,
                                     y,
                                     gamma,
                                     split.out,
                                     alpha = 0.05,
                                     family,
                                     verbose=FALSE)
{

  ## signif.cluster## this is implicitly known since we only go through the cluster points
  ## of which the parents were possible to reject

  ## But what to do with the branches that were extinct? need to save them too someway
  extinct <- c()
  ## the current clusters we are about to test
  c.test <- 1## the root node

  tree.pvals <- matrix(NA,
                       length(tree),## number of cluster nodes
                       length(split.out))## number of splits
  clusters.tested <- c(1)
  pvals.tested <- c(-1)
  while(TRUE)
    {
      found.signif <- FALSE
      ## loop through the c.test depth first
      c.test.new <- c.test## we aren't supposed to change the set we are looping over
      for(c in c.test)
        {
          if(is.na(tree.pvals[c,1]))
            {## haven't calculated the p-values yet for this node
              pvalues <- signif.partialFtest(tree=tree,
                                             c=c,
                                             y=y,
                                             x=x,
                                             split.out=split.out,
                                             family=family)
              tree.pvals[c,] <- pvalues## save them for later, if we get back to the same nodes
            }else{
              pvalues <- tree.pvals[c,]
            }
          ## Now determine the weights to use to determine significance
          ## 4.2 INFERENCE OF HIERARCHICALLY ORDERED CLUSTERS OF VARIABLES
          mult.adj <- getmultadj.fortree(tree=tree,
                                         c=c,
                                         split.out=split.out,
                                         extinct=extinct,
                                         c.test=c.test.new)
          pvalues.adj <- pvalues*mult.adj
          ## Then do the aggregation as for multisplit
          ##(Based on multi-split.R code)
          ##(should extract this as a separate function)
          quant.gamma <- quantile(pvalues.adj, gamma) / gamma
          if(length(gamma) > 1)
            penalty <- (1 - log(min(gamma)))
          else
            penalty <- 1

          pvals.pre <- min(quant.gamma) * penalty
          agg.pvalue <- pmin(pvals.pre, 1)

          ## track the pval
          pvals.tested[match(c,clusters.tested)] <- agg.pvalue

          if(agg.pvalue <= alpha)## add the children
            {## depth first search
              if(is.na(tree[[c]]$child1))
                {
                  stopifnot(is.na(tree[[c]]$child2))## the second child should be NA too
                  ## OLD## extinct <- c(extinct,c)
                  ## TODO add all parents of an extinct branch that are extinct too
                  extinct <- update.extinct(nowextinct=c,
                                            tree=tree,
                                            c.test=c.test.new,
                                            oldextinct=extinct)
                }else{
                  clusters.tested <- c(clusters.tested,
                                       tree[[c]]$child1,
                                       tree[[c]]$child2)
                  c.test.new <- c(tree[[c]]$child1,
                                  tree[[c]]$child2,
                                  c.test.new)
                }
              c.test.new <- setdiff(c.test.new,c)## remove the current tested cluster that was found significant from the to test
              found.signif <- TRUE
              ## Note: can optimize the c.test further with the
            }
        }
      ## never remove clusters you did not found significant from c.test,
      ## we will revisit it later potentially

      ## stop looping through the c.test if we have looped once through c.test and nothing was found significant
      c.test <- c.test.new
      if(verbose) {
        cat("length of c.test:", length(c.test), "\n")
        cat("c.test itself:", c.test, "\n")
      }
      if(!found.signif)
        break
    }
  message("doing final post processing")
  ## will still have to transform the output to be able to match the one from lowerbound etc.

  ## out <- list()
  ## out$c.test <- c.test
  ## out$extinct <- extinct
  ## out$tree <- tree

  ## TODO convert to the lowerbound output thingie and you are done!
  ## ---- check out the tree plot and verify correctness
  ## out <- list()
  ## signif.clusters <- get.signif.clusters(tree=tree,
  ##                                        c.test=c.test,
  ##                                        extinct=extinct)
  ## can we also only return the significant clusters?

  clusters <- list()
  pval <- list()
  leftChild <- list()
  rightChild <- list()
  k <- 0
  ## we somehow need to all clusters we have tested
  ## OLDOLDOLD
  ##  for(i in 1:length(signif.clusters))
  ## {
  ##   ct <- signif.clusters[i]
  ## OLDOLDOLD
  for(i in 1:length(clusters.tested))
    {
      ct <- clusters.tested[i]
      k <- k+1
      clusters[[k]] <- tree[[ct]]$vars

      pval[[k]] <- pvals.tested[i]##
      ## if(ct %in% signif.clusters)
      ##   pval[[k]] <- 0
      ## else
      ##   pval[[k]] <- 1
      leftChild[[k]] <- match(tree[[ct]]$child1,clusters.tested)
      rightChild[[k]] <- match(tree[[ct]]$child2,clusters.tested)

      ## lc <- leftChild[[k]]
      ## rc <- rightChild[[k]]

      ## if(!is.na(lc)) {
      ##   k <- k+1
      ##   clusters[[k]] <- tree[[lc]]$vars
      ##   pval[[k]] <- 1## not significant
      ##   leftChild[[k]] <- tree[[lc]]$child1
      ##   rightChild[[k]] <- tree[[lc]]$child2
      ## }

      ## if(!is.na(rc)) {
      ##   k <- k+1
      ##   clusters[[k]] <- tree[[rc]]$vars
      ##   pval[[k]] <- 1## not significant
      ##   leftChild[[k]] <- tree[[rc]]$child1
      ##   rightChild[[k]] <- tree[[rc]]$child2
      ## }
    }
  lCh <- unlist( leftChild); lCh[is.na(lCh)] <- -1
  rCh <- unlist(rightChild); rCh[is.na(rCh)] <- -1
  ## return
  list(
    clusters = clusters,
    pval = unlist(pval),
    leftChild  = lCh,
    rightChild = rCh,
    alpha = alpha,
    hh = hh)
  ## TODO next steps are doing 4.2.1 inheritance procedure
  ## TODO afterwards try 4.1 exploiting logical relationships Shaffer improvements
}

createtree.from.hclust <- function(hh,
                                   verbose=TRUE)
{## Note: an object oriented approach probably has better performance
  p <- length(hh$order)

  ## the root node
  tree <- list()## the implicit numbering is the node number
  tree[[1]] <- list(child1=NA,
                    child2=NA,
                    parent=NA,
                    vars=which(cutree(hh,1)==1))## this is equivalent to 1:p but conserves variable names

  parent <- rep(1,p)

  if(verbose)
      cat("Converting the hclust output to a tree data structure\n",
          paste0("building the tree is ",0,"% done.\n"))
  nn <- 2 ## this is the running node number to add to the tree
  for(i in 2:nrow(hh$merge)) {
      if(verbose && (i %% round(nrow(hh$merge)/4)==0))
          cat("building the tree is ", round(i/round(nrow(hh$merge))*100),"% done.\n")
      old.cut <- cutree(hh,i-1)
      new.cut <- cutree(hh,i)

      tparent <- NA

      ## find out what cluster has been split up FAST WAY
      tparent <- parent[min(which(old.cut!=new.cut))]
      c <- (parent==tparent)

      ## ## Identification approach, clusters themselves
      ## old.clusters <- split(1:p,old.cut)
      ## new.clusters <- split(1:p,new.cut)

      ## ## I actually only have to do this for the clusters of cluster size that is
      ## ## not in the new clusters :)
      ## old.table <- table(old.cut)## table(parent)
      ## new.table <- table(new.cut)

      ## ## alternative OLD
      ## old.disappeared <- as.numeric(names(old.table)[old.table==setdiff(old.table,intersect(old.table,new.table))])
      ## if(length(old.disappeared)==0)
      ##   {## there are duplicate cluster sizes
      ##     ## just check out the duplicate ones
      ##     duplicate.c.sizes <- as.numeric(names(table(old.table))[table(old.table)>1])
      ##     old.disappeared <- NA
      ##     for(d in duplicate.c.sizes)
      ##       {
      ##         if(sum(old.table==d)>sum(new.table==d))
      ##         {
      ##           old.disappeared <- which(old.table==d)
      ##           break
      ##         }
      ##       }
      ##   }
      ##
      ## if(length(old.disappeared)==1)
      ##   {
      ##     c <- old.clusters[[old.disappeared]]
      ##     tparent <- unique(parent[c])
      ##   }else{
      ##     for(c in old.clusters[old.disappeared])
      ##       {
      ##         if(any(sapply(sapply(new.clusters,identical,c),isTRUE)))
      ##           {## there is a new cluster identical to the c
      ##             ##==> the cluster c was not split up
      ##           }
      ##         else{
      ##           tparent <- unique(parent[c])
      ##           stopifnot(length(tparent)==1)
      ##           break
      ##         }
      ##       }
      ##   }

      browser(expr=is.na(tparent))
      browser(expr=(length(tparent)>1))
      stopifnot(!is.na(tparent))

      ## update the parent
      stopifnot(is.na(tree[[tparent]]$child1) && is.na(tree[[tparent]]$child2))
      tree[[tparent]]$child1 <- nn
      tree[[tparent]]$child2 <- nn+1


      ## create child 1
      child1.c <- which(new.cut == unique(new.cut[c])[1])
      tree[[nn]] <- list(child1=NA,
                         child2=NA,
                         parent=tparent,
                         vars=child1.c)
      parent[child1.c] <- nn
      nn <- nn+1
      ## create child 2
      child2.c <- which(new.cut == unique(new.cut[c])[2])
      tree[[nn]] <- list(child1=NA,
                         child2=NA,
                         parent=tparent,
                         vars=child2.c)
      parent[child2.c] <- nn
      nn <- nn+1
  }
  if(verbose)
    cat("building the tree is ",100,"% done.\n", sep="")
  tree
} ## {createtree.from.hclust}
