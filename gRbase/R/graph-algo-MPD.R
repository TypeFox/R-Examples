# R function returning a junction tree representation of the MPD of an undirected graph
# Author: Clive Bowsher

# Inputs:
# uG: graphNEL representation of the undirected graph

# Output: returns a junction tree representation of the MPD of uG using the Recursive Thinning and Aggregate Cliques
#         algorithms of Olesen & Madsen 2002. Names(MPDTree[[r]]) retains the original clique numbers from the RIP
#         ordering of cliques of the minimal triangulated graph returned by rip(TuG) below (for r=2,3,4)

# date: 23.03.09
# checked: using the DAG of Figs 2(Asia) & 9 of Olesen & Madsen 2002, Fig1s and Fig2 of Leimer93 (CGB 09.03.09).
# 		   also using the known MPD of ugraph(kLIGnelChap), and line-by-line on 19.03.09
# issues: none known


mpd <- function(object, tobject=minimalTriang(object), details=0) {
    UseMethod("mpd")
}

mpd.default <- function(object, tobject=triangulate(object), details=0){
    cls <- match.arg(class( object ),
                     c("graphNEL","matrix","dgCMatrix"))
    switch(cls,
           "graphNEL" ={.mpd(object, TuG=tobject, details=details) },
           "dgCMatrix"=,
           "matrix"   ={ #FIXME: minimalTriang: Not sure if this is correct...
               object2 <- as(object,  "graphNEL")
               tobject2<- as(tobject, "graphNEL")
               .mpd(object2, TuG=tobject2, details=details)
           })
}




mpdMAT <- function(amat, tamat=minimalTriangMAT(amat), details=0){
  .mpd(as(amat, "graphNEL"), TuG=as(tamat, "graphNEL"), details=details)
}

.mpd <- function(uG, TuG=minimalTriang(uG), details=0) {

  ##TuG <- MinimalTriang(uG)

  Tree <- rip(TuG) ## Tree arranged as [1]nodes,[2]cliques,[3]separators,[4]parents
  #print(Tree)

  MPDTree <- vector("list",length=4)

  if (length(Tree[[2]]) == 1)
    {
      MPDTree <- Tree		## Tree contains only 1 clique so cannot be decomposed
    }
  else
    {
      i <- length(Tree[[3]])	## indexes separators of Tree, starting with the last one (although the order of aggregation is immaterial)
      while (i > 1)		## 1st clique never has a parent
        {
          if(!(is.complete(graph::subGraph(c(Tree[[3]][i],recursive=TRUE),uG))))
            {
              Tree[[3]][i] <- 0		# set separator i to 'empty'
              parent.i <- Tree[[4]][[i]]
              #cat(sprintf("i=%i parent.i=%i\n", i, parent.i))
              sel <- !(Tree[[2]][[i]] %in% Tree[[2]][[parent.i]])
              Tree[[2]][[parent.i]] <- c(Tree[[2]][parent.i],Tree[[2]][[i]][c(sel)],recursive=TRUE)	# merge clique i into its parent clique

              if (details>0)
                cat(sprintf("Clique %i merged into clique %i\n", i, parent.i))

              Tree[[2]][[i]] <- 0		# set clique i to 'empty'
              Tree[[4]][[i]] <- 0		# i has been emptied, hence has no parent
              if(i < length(Tree[[3]]))	# last cluster can have no children
		{
                  for (j in (i+1):length(Tree[[4]]))  		# determine previous children of i, and make them children of parent.i
                    {
                      if(Tree[[4]][[j]] == i)
                        {
                          Tree[[4]][[j]] <- parent.i
                        }
                    }
		}
            }
          i <- i - 1
        }

      MPDTree <- Tree

      sel <- matrix(0,length(MPDTree[[2]]),2)
      for (i in 1:length(MPDTree[[2]]))
        {
          if(0 == MPDTree[[2]][[i]][1])
            {
              sel[i,1] <- FALSE # no change to sel
            }
          else
            {
              sel[i,1] <- TRUE	## to select only those cluster numbers that are now not empty
              sel[i,2] <- i	## record original cluster number i
            }
        }

      ch    <- as.logical(sel[,1])
      Appel <- sel[,2][ch]	## original cluster numbers of non-empty clusters
      MPDTree[[2]] <- MPDTree[[2]][ch]
      names(MPDTree[[2]]) <- Appel
      MPDTree[[3]] <- MPDTree[[3]][ch]
      names(MPDTree[[3]]) <- Appel
      MPDTree[[4]] <- MPDTree[[4]][ch]
      names(MPDTree[[4]]) <- Appel

      return(MPDTree)
    }
}
