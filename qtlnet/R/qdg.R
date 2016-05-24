##############################################################################
##
## $Id: codeQDG.R,v 2007/11/28 byandell Exp $
##
##     Copyright (C) 2007 Elias Chaibub Neto and Brian S. Yandell
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the
## Free Software Foundation; either version 2, or (at your option) any
## later version.
##
## These functions are distributed in the hope that they will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## The text of the GNU General Public License, version 2, is available
## as http://www.gnu.org/copyleft or by writing to the Free Software
## Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
##
## Routines: qdg, summary.qdg, print.qdg
##           qdg.perm.test
##           summary.qdg.sem, print.qdg.sem, summary.qdg.perm.test, print.qdg.perm.test
##############################################################################

## Results of qdg can be plotted using graph.qdg (using igraph).

qdg <- function(cross, phenotype.names, marker.names, QTL, alpha, 
                n.qdg.random.starts, addcov = NULL, intcov = NULL,
                skel.method=c("pcskel","udgskel"), udg.order = 2)
{
  if(!inherits(cross, "cross"))
    stop("cross must be an object of class cross")

  pheno.data <- cross$pheno[,phenotype.names]
  skel.method <- match.arg(skel.method)
  if(skel.method == "pcskel"){
    ## Create skeleton using R/pcalg package.
    suffStat <- list(C = cor(pheno.data), n = nrow(pheno.data))
    pcskeleton <- skeleton(suffStat, gaussCItest, p = ncol(pheno.data), alpha = alpha)

    ## Transform to UDG.
    UDG <- transformPCtoUDG(pcskeleton)
    UDG <- renameUDG(selpheno = phenotype.names, UDG = UDG)
  }
  else if (skel.method == "udgskel") 
    UDG <- approximate.UDG(Data = pheno.data, alpha = alpha, fixed.order = udg.order)

  DG <- orient.graph.edges(cross=cross,UDG=UDG,QTLs=QTL,addcov=addcov,intcov=intcov)
  rc <- recheck.directions(cross=cross,QTLs=QTL,oldDG=DG,addcov=addcov,intcov=intcov)
  aux.cross <- argmax.geno(cross)
  genotypes <- pull.geno.argmax(aux.cross)
  as <- get.all.solutions(DG=DG, rc=rc, n.shuffles=n.qdg.random.starts, cross=cross,
                          QTLs=QTL, markers=marker.names, phenotypes=phenotype.names,
                          genotypes=genotypes, addcov=addcov, intcov=intcov)
  best <- which(as$BIC == min(as$BIC))
  mylist <- list(UDG, DG, best, as)
  names(mylist) <- c("UDG","DG","best.lm","Solutions")
  mylist$marker.names <- marker.names
  mylist$phenotype.names <- phenotype.names
  mylist$addcov <- addcov
  class(mylist) <- c("qdg", "list")
  
  mylist
}

#################################################
summary.qdg <- function(object, ...)
{
  cat("\n Number of solutions:\n")
  print(length(object$Solutions$BIC))
  cat("\nBest solution:\n")
  print(object$Solutions$solutions[[object$best.lm]])
  bic.lm <- object$Solutions$BIC[object$best.lm]
  cat("\nBIC:\n")
  print(c(lm = bic.lm))
  cat("\nBest solution is solution number:\n")
  print(object$best.lm)
  cat("\nCaution:\n")
  print("If one of the solutions is a cyclic graph you should run qdg.sem in order to score the networks using SEM.")
  invisible()
}

#################################################
print.qdg <- function(x, ...) summary(x, ...)

#################################################
transformPCtoUDG <- function(PC)
{
  edges <- graph::edges(PC@graph)  
  tmp1 <- rep(names(edges), sapply(edges, length))
  tmp2 <- unlist(edges)
  UDG <- data.frame(matrix(NA,length(tmp1), 2))
  UDG[,1] <- as.numeric(tmp1)
  UDG[,2] <- as.numeric(tmp2)
  UDG <- UDG[as.numeric(tmp1) < as.numeric(tmp2), ]
  names(UDG) <- paste("node", 1:2, sep = "")
  UDG$edge <- 1
  class(UDG) <- c("qdg", "data.frame")
  attr(UDG, "edgemode") <- "undirected"
  attr(UDG, "message") <- ""
  attr(UDG, "cont") <- 0
  UDG
}
###################################
renameUDG <- function(selpheno,UDG)
{
  rUDG <- UDG
  n <- length(UDG[,1])
  for(i in 1:n){
    rUDG[i,1] <- selpheno[UDG[i,1]] 
    rUDG[i,2] <- selpheno[UDG[i,2]]
  }
  rUDG
}
##############################################################
approximate.UDG <- function(Data, alpha, fixed.order = 2)
{
  partial.correlation <- function(i, j, k, comb, R){
    RR <- R[c(i, j, comb[, k]), c(i, j, comb[, k])]
    RRinv <- solve(RR)
    D <- diag(1/sqrt(diag(RRinv)))
    return(-D%*%RRinv%*%D)
  }
  pvalue <- function(correlation, n, np){
    tobs <- correlation/sqrt((1-correlation^2)/(n-2-np))
    return(2*pt(abs(tobs), df = n-2-np, lower.tail = FALSE))	
  }
  n <- length(Data[, 1])
  nv <- length(Data[1, ])
  aux.comb <- c(1:nv)
  R <- cor(Data, method = "spearman")
  UDG <- data.frame(matrix(1, nv*(nv-1)/2, 3))
  names(UDG) <- c("node1", "node2", "edge")
  cp <- 1
  for(i in 1:(nv-1)){
    for(j in (i+1):nv){
      UDG[cp, "node1"] <- names(Data)[i]
      UDG[cp, "node2"] <- names(Data)[j]
      order0.cor <- R[i, j]
      pv <- pvalue(correlation = order0.cor, n = n, np = 0)
      order0.ht <- ifelse(pv > alpha, 0, 1)
      if(order0.ht == 0) UDG[cp, "edge"] <- 0
      else{
        if(nv > 2) end <- 0 
        else end <- 1
        order <- 1
        while((order <= fixed.order) & (end == 0)){
          comb <- combn(x = aux.comb[-c(i, j)], m = order)
          nc <- length(comb[1, ])
          k <- 1
          while((k <= nc) & (end == 0)){
            PR <- partial.correlation(i = i, j = j, k = k, comb = comb, R = R)
            pv <- pvalue(correlation = PR[1, 2], n = n, np = order)
            ht <- ifelse(pv > alpha, 0, 1)
            if(ht == 0){
              UDG[cp, "edge"] <- 0
              end <- 1
            }
            k <- k+1
          }
          order <- order+1
          if(nv <= order+1){end <-1}				
        }
      }
      cp <- cp + 1
    }
  }
  class(UDG) <- c("qdg", "data.frame")
  attr(UDG, "edgemode") <- "undirected"
  attr(UDG, "message") <- ""
  attr(UDG, "cont") <- 0
  
  UDG
}
