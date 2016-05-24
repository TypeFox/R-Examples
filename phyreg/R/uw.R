uw <-
function(phy,spu,nspec) { # called by getphytxp

  lphy<-length(phy);
  spc<-c(rep(1,nspec),rep(0,lphy-nspec)) 
  spn<-c(spu,rep(0,lphy-nspec))
  
  isd<-rep(0,lphy+1)
  for (ii in 1:lphy) isd[phy[ii]]<-isd[phy[ii]]+isd[ii]+spn[ii]  # This gives the total number of included species descendants
  isd[1:nspec]<-spu                                                              #  and now including itself in the descendants

 # A node is to be called "alive" if it has included species descendants, thus
 alive<-(isd>=1)
 
 # The number of daughters alive of a given higher node is important
 ndalive<-rep(0,lphy+1)
 for (ii in 1:lphy) ndalive[phy[ii]]<-ndalive[phy[ii]]+alive[ii]
 
 # A node is to be kept in the new phylogeny if it (1) is a species node OR (2) has two alive daughters, thus
 
 keep<-(ndalive>=2L) | c(spc,0) 
  
 # We're going to assign the closest ancestor who's to be kept. We require an ancestor of the root, so its given the next number up
 
 bigphy<-c(phy,1+max(phy))
 
 for (ii in lphy:1) if (!keep[bigphy[ii]]) bigphy[ii]<-bigphy[bigphy[ii]]
 
 # We proceed by simplifying bigphy by multiplying it by alive and keep. Then it contains non-zero entries at only the kept nodes
 #  of the old phylogeny, and it has the (old) names of the parent. I think this *just* makes it easier to read.
 
 bigphy<-bigphy*alive*keep

 # bigphy now has the MRA of each kept and alive node. So every kept and alive node now has a kept parent.
 # We now want to have cumulative keep, which is the new name for a given node and zero elsewhere, and the length of the 
 #  new phylogeny, which is one less than the total number of nodes
 
 # Now, keep indicates for each of the original nodes whether its to be kept, but not for the extra top node -- however that is
 #  never kept, as it cannot have two alive daughters. But it may be the (old) name of the highest node if the original top
 #  is not being kept.
 
 # To define the names of the new nodes as indexed by the old nodes, we define cumkeep, and we need to add a final
 #  element for the extra top. The new phylogeny has length sum(keep) minus one if the old root remains, but minus two if
 #  the old root has gone, for then keep contains not only the old root but also its parent.  # 2014_01_31: still "true"?
 
 cumkeep<-cumsum(keep)*keep
 cumkeep<-c(cumkeep,1+max(cumkeep))
 # if (keep[length(keep)]) {newlphy<-sum(keep)-1} else {newlphy<-sum(keep)-2} # Seems not to be necessary:
 newlphy<-sum(keep)-1
 # We can now create a new phylogeny with zero for non-included species, and each higher node having at least two daughters
 
 newphy<-rep(-1,newlphy)
 
 for (ii in 1:lphy) if (keep[ii] & alive[ii] & (cumkeep[ii]<=newlphy)) newphy[cumkeep[ii]]<-cumkeep[bigphy[ii]] 
                                                                                                        # newphy now has the new name of new parent of new name of ii
 newphy[1:nspec]<-newphy[1:nspec]*spu                                                                        # but non-included species are set to zero
 
 # That's the new phylogeny made, but we also need the original names of each node
 #  So we invert cumkeep...
 
 on<-rep(0,newlphy+1)
 for (ii in 1:(lphy+1)) if (cumkeep[ii] != 0) on[cumkeep[ii]]<-ii
  
 return(list(txp=newphy,on=on))
 
 }
