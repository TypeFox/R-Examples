#******************************************************************************* 
#
# Bayesian Regression and Adaptive Sampling with Gaussian Process Trees
# Copyright (C) 2005, University of California
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (rbgramacy@ams.ucsc.edu)
#
#*******************************************************************************


"tgp.postprocess" <-
function(ll, Xnames, response, pred.n, zcov, Ds2x, improv, sens.p, Zm0r1, params, rmfiles=TRUE)
{
  ## deal with X, and names of X, as well as Xsplit
  ll$X <- framify.X(ll$X, Xnames, ll$d)
  ll$Xsplit <- framify.X(ll$Xsplit, Xnames, ll$d)
  ll$nsplit <- NULL
  
  ## deal with Z, and names of Z
  if(is.null(response)) ll$response <- "z"
  else ll$response <- response

  ## remove from the list if not requested
  if(Ds2x == FALSE) { ll$Ds2x <- NULL; }
  if(improv == FALSE || is.null(improv)) { ll$improv <- NULL; }
  
  ## deal with predictive data locations (ZZ)
  if(ll$nn == 0 || (ll$BTE[2]-ll$BTE[1])==0 || !is.null(sens.p)) { 
    ll$XX <- ll$ZZ.mean <- ll$ZZ.s2 <-  ll$ZZ.q <- ll$ZZ.km <- ll$ZZ.ks2 <- ll$ZZ.vark <- NULL
    ll$ZZ.q1 <- ll$ZZ.med <- ll$ZZ.q2 <- ll$ZpZZ.s2 <- ll$Ds2x <- ll$improv <- NULL
  } else {
    ## do predictive input/output processing
    
    ## replace NaN's in improv with zeros
    ## shouldn't happen because check have been moved to C code
    if((!is.null(ll$improv)) && sum(is.nan(ll$improv) > 0)) {
      warning(paste("encountered", sum(is.nan(ll$improv)),
                    "NaN in Improv, replaced with zeros"), call.=FALSE)
      ll$improv[is.nan(ll$improv)] <- 0
    }
    
    ## make sure XX has the correct output format
    ll$XX <- framify.X(ll$XX, Xnames, ll$d)
  }

  ## turn improv into a data.frame where the second column is the rankings
  if(!is.null(improv)){
    ll$irank[ll$irank == 0] <- ll$nn
    ll$improv <- data.frame(improv=ll$improv, rank=ll$irank)
  }
  ll$irank <- NULL

  ## NULL-out data-predictive output if unused
  if(pred.n == FALSE || ll$BTE[2]-ll$BTE[1] == 0) {
    ll$Zp.mean <- ll$Zp.q <- ll$Zp.q1 <- ll$Zp.q2 <- NULL;
    ll$Zp.s2 <- ll$ZpZZ.s2 <-  ll$Zp.km <-  ll$Zp.vark <- ll$Zp.ks2 <- ll$Zp.med <- NULL
  }
  
  ## gather information about partitions
  if(file.exists(paste("./", "best_parts_1.out", sep=""))) {
    ll$parts <- as.matrix(read.table("best_parts_1.out"))
    if(rmfiles) unlink("best_parts_1.out")
  } else { ll$parts <- NULL }
  
  ## gather information about MAP trees as a function of height
  ll$trees <- tgp.get.trees(ll$Xsplit, rmfiles)
  ll$posts <- read.table("tree_m0_posts.out", header=TRUE)
  if(ll$BTE[2] - ll$BTE[1] == 0) ll$posts <- NULL
  if(rmfiles) unlink("tree_m0_posts.out")

  ## read the trace in the output files, and then delete them
  if(ll$trace) ll$trace <- tgp.read.traces(ll$n, ll$nn, ll$d, params$corr, ll$verb, rmfiles)
  else ll$trace <- NULL
  
  ## store params
  ll$params <- params

  ## clear the verb, state, tree and MAP fields for output
  ll$verb <- NULL; ll$state <- NULL; ll$tree <- NULL; ll$MAP <- NULL; ll$nt <- NULL
  ll$ncol <- NULL; ll$hier <- NULL;

  ## clear output dimensions
  ll$pred.n <- ll$nnprime <- ll$krige <- ll$bDs2x <- NULL

  ## consolidate itemps
  nt <- as.integer(ll$itemps[1])
  lambda <- ll$itemps[length(ll$itemps)]
  if(lambda == 1) lambda <- "opt"
  else if(lambda == 2) lambda <- "naive"
  else if(lambda == 3) lambda <- "st"
  else stop(paste("bad lambda = ", lambda, sep=""))
  ll$itemps <- list(c0n0=as.integer(ll$itemps[2:3]), k=ll$itemps[4:(nt+3)],
                    pk=ll$itemps[(nt+4):(2*nt+3)], 
                    counts=as.integer(ll$itemps[(2*nt+4):(3*nt+3)]),
                    lambda=lambda)

  ## consolidate ess
  if(nt == 1) ll$ess <- ll$ess[1]
  else {
    ll$ess=list(combined=ll$ess[1],
      each=data.frame(k=ll$itemps$k, count=ll$ess[2:(nt+1)], ess=ll$ess[(nt+2):(2*nt+1)]))
  }
  

  ## change {0,1} to {TRUE,FALSE}
  if(ll$linburn) ll$linburn <- TRUE
  else ll$linburn <- FALSE

  ## pretty-up the grow, prune, change and swap stats
  ll$gpcs[is.nan(ll$gpcs)] <- NA
  ll$gpcs <- data.frame(t(ll$gpcs))
  names(ll$gpcs) <- c("grow", "prune", "change", "swap")

  ## deal with sensitivity analysis outputs
  if(!is.null(sens.p)){
    names(sens.p) <- NULL
    sens.par <- list(nn.lhs=sens.p[1], rect=matrix(sens.p[2:(ll$d*2+1)], nrow=2),
                 shape=sens.p[(ll$d*2+2):(ll$d*3+1)], mode=sens.p[(ll$d*3+2):(ll$d*4+1)],
                 ngrid=ll$sens.ngrid, span=ll$sens.span)
    sens <- list()
    sens$par <- sens.par
    sens$ngrid <- NULL
    sens$span <- NULL
    sens$Xgrid <- matrix(ll$sens.Xgrid, ncol=ll$d)
    sens$ZZ.mean <- matrix(ll$sens.ZZ.mean, ncol=ll$d)
    sens$ZZ.q1 <- matrix(ll$sens.ZZ.q1, ncol=ll$d)
    sens$ZZ.q2 <- matrix(ll$sens.ZZ.q2, ncol=ll$d)
    sens$S <- matrix(ll$sens.S, ncol=ll$d, byrow=TRUE)
    sens$T <- matrix(ll$sens.T, ncol=ll$d, byrow=TRUE)
  } else{ sens <- NULL }

  ## clear ll$sens.* and replace with single list
  ll$sens.Xgrid <- ll$sens.ZZ.mean <- ll$sens.ZZ.q1 <- ll$sens.ZZ.q2 <- NULL
  ll$sens.ngrid <- ll$sens.span <- ll$sens.S <- ll$sens.T <- NULL
  ll$sens <- sens
  
  ## undo mean0.range1
  if(!is.null(Zm0r1)) {
    ll$Z <- undo.mean0.range1(ll$Z,Zm0r1$undo)
    ll$Zp.mean <- undo.mean0.range1(ll$Zp.mean,Zm0r1$undo)
    ll$ZZ.mean <- undo.mean0.range1(ll$ZZ.mean,Zm0r1$undo)
    ll$Zp.km <- undo.mean0.range1(ll$Zp.km,Zm0r1$undo)
    ll$ZZ.km <- undo.mean0.range1(ll$ZZ.km,Zm0r1$undo)
    ll$Zp.vark <- undo.mean0.range1(ll$Zp.vark,Zm0r1$undo, nomean=TRUE, s2=TRUE)
    ll$ZZ.vark <- undo.mean0.range1(ll$ZZ.vark,Zm0r1$undo, nomean=TRUE, s2=TRUE)
    ll$Zp.ks2 <- undo.mean0.range1(ll$Zp.ks2,Zm0r1$undo, nomean=TRUE, s2=TRUE)
    ll$ZZ.ks2 <- undo.mean0.range1(ll$ZZ.ks2,Zm0r1$undo, nomean=TRUE, s2=TRUE)
    ll$ZpZZ.ks2 <- undo.mean0.range1(ll$ZpZZ.ks2,Zm0r1$undo, nomean=TRUE, s2=TRUE)
    ll$Zp.q <- undo.mean0.range1(ll$Zp.q,Zm0r1$undo, nomean=TRUE)
    ll$ZZ.q <- undo.mean0.range1(ll$ZZ.q,Zm0r1$undo, nomean=TRUE)
    ll$Zp.s2 <- undo.mean0.range1(ll$Zp.s2,Zm0r1$undo, nomean=TRUE, s2=TRUE)
    ll$ZZ.s2 <- undo.mean0.range1(ll$ZZ.s2,Zm0r1$undo, nomean=TRUE, s2=TRUE)
    ll$Zp.q1 <- undo.mean0.range1(ll$Zp.q1,Zm0r1$undo)
    ll$Zp.med <- undo.mean0.range1(ll$Zp.med,Zm0r1$undo)
    ll$Zp.q2 <- undo.mean0.range1(ll$Zp.q2,Zm0r1$undo)
    ll$ZZ.q1 <- undo.mean0.range1(ll$ZZ.q1,Zm0r1$undo)
    ll$ZZ.med <- undo.mean0.range1(ll$ZZ.med,Zm0r1$undo)
    ll$ZZ.q2 <- undo.mean0.range1(ll$ZZ.q2,Zm0r1$undo)
    for(j in 1:ll$d){
      ll$sens.ZZ.mean[,j] <- undo.mean0.range1(ll$sens.ZZ.mean[,j],Zm0r1$undo)
      ll$sens.ZZ.q1[,j] <- undo.mean0.range1(ll$sens.ZZ.q1[,j],Zm0r1$undo)
      ll$sens.ZZ.q2[,j] <- undo.mean0.range1(ll$sens.ZZ.q2[,j],Zm0r1$undo)
    }
    ll$m0r1 <- TRUE
  } else { ll$m0r1 <- FALSE }

  ## turn Z*.s2 into a matrix (covariance matrix)
  if(!is.null(ll$Zp.s2) && ll$zcov) ll$Zp.s2 <- matrix(ll$Zp.s2, ncol=ll$n)
  if(!is.null(ll$ZZ.s2) && ll$zcov) ll$ZZ.s2 <- matrix(ll$ZZ.s2, ncol=ll$nn)
  if(!is.null(ll$ZpZZ.s2) && ll$zcov) ll$ZpZZ.s2 <- t(matrix(ll$ZpZZ.s2, ncol=ll$n))
  else ll$ZpZZ.s2 <- NULL
  ll$zcov <- NULL

  ## set class information and return
  class(ll) <- "tgp"
  return(ll)
}


"tgp.get.trees" <-
function(X, rmfiles=TRUE)
{
  trees <- list()

  ## get all of the names of the tree files
  tree.files <- list.files(pattern="tree_m0_[0-9]+.out")

  ## return no trees if the run was only burn-in
  if(length(tree.files) == 0) return(NULL)

  ## for each tree file
  for(i in 1:length(tree.files)) {

    ## grab the height from the filename
    h <- as.numeric(strsplit(tree.files[i], "[_.]")[[1]][3])
    
    ## read it in, then remove it
    trees[[h]] <- read.table(tree.files[i], header=TRUE)
    if(rmfiles) unlink(tree.files[i])

    ## correct the precision of the val (split) locations
    ## by replacing them with the closest X[,var] location
    if(nrow(trees[[h]]) == 1) next;
    nodes <- (1:length(trees[[h]]$var))[trees[[h]]$var != "<leaf>"]
    for(j in 1:length(nodes)) {
	col <- as.numeric(as.character(trees[[h]]$var[nodes[j]])) + 1
      m <- which.min(abs(X[,col] - trees[[h]]$val[nodes[j]]))
      trees[[h]]$val[nodes[j]] <- X[m,col]
    }          
  }
  
  return(trees)
}
