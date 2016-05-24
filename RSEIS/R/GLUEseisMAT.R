`GLUEseisMAT` <-
function(GFIL)
      {
	###  find duplicated stations in a matrix and
        ###   fill in the traces that are continuations
	### return the new matrix and the vector duplicates
	dot = which(duplicated(GFIL$KNOTES))
	G = GFIL$JMAT
	for(i in 1:length(dot))
	  {
	    w = which(!is.na(match(GFIL$KNOTES, GFIL$KNOTES[dot[i]])))
	    
	    a = G[,w[1]]
	    a[!is.na(G[,w[2]])] = G[!is.na(G[,w[2]]), w[2]]
	    G[,w[1]]  = a 
	    
	  }
	invisible(list(JMAT=G, dpl=dot) )
	
      }

