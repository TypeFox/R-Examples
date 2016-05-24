genTrees <-
function(x, N = 200, filebase = 'trial', method = c('nni', 'random'), maxmoves = 2, perms = 'DEFAULT', software = c('raxml'), ...) {
  ## Arguments:
  ## x = phylo tree
  ## N = total number of trees to generate
  ## filebase = file name base; a tree file (.tre) and paup command file (.nex) will be created for both
  ## method = method for generating trees
  ## maxmoves = maximum number of rearrangements per tree for nni or spr
  ## perms = number of permutations per maxmoves
  ## ... = additional arguments to pass along to rtree or rNNI
  ## works with nni, 12 nov 10
  ## January 2014: as written, this doesn't unroot the tree. It ought to, unless you are evaluating trees in a rooted framework (e.g., not using GTR)
  ## 20 January 2014: updated to make sure all trees are unique
  if(perms == 'DEFAULT') perms <- c(length(nni(x)), max(1, as.integer(N-(length(nni(x))))))
  if(class(x) != 'phylo') stop('This function requires a phylo object as its first argument')
  if(method[1] == 'nni') {
	for(i in seq(maxmoves)) {
	  message(paste('doing maxmoves', i))
	  if(i == 1) {
	    originalTree <- list(x)
		treeset <- c(originalTree, lapply(nni(x), function(x) x))
		}
	  # else treeset <- c(treeset, rNNI(x, i, perms[i]))
	  else treeset <- c(treeset, lapply(unique(rNNI(x, i, perms[i] * 1.5)), function(x) x))
	  treeset <- lapply(treeset, unroot)
	  class(treeset) <- 'multiPhylo'
	  if(length(treeset) >= sum(perms[1:i], 1)) treeset <- unique(treeset)[1:sum(perms[1:i], 1)]
	  else(warning(paste('Treeset only includes', length(treeset), 'trees of the', sum(perms[1:i], 1), 'expected')))
	  treeset <- treeset[!sapply(treeset, is.null)]
	  # just takes the first set of uniques... chops off non-uniques presented so far
      }	# end i
	} # end if	  
  else if(method[1] == 'random') treeset = c(x, rtreePhylo(x, N, ...))
  class(treeset) <- 'multiPhylo'
  if(software[1] == 'raxml') {
	message('writing raxml')
	write.tree(treeset, file = paste(filebase, '.trees.tre', sep = ''))
    # write.tree(x, file = paste(filebase, '.optimal.tre', sep = '')) ## no longer separating optimal from full trees
    }
  return(treeset)
  }
