"getPrel" <-
function (model) 
{
	# dependent params. are specified as list or vector
	# in prel; get their indices in the parameter group 
        # in its vectorized form 

  prelspec <- model@prelspec
  ppars <- intersect(slotNames(theta()), slotNames(model))
  pinde <- vector("list", length(ppars))
  for(p in 1:length(pinde)) pinde[[p]] <- vector()
  names(pinde) <- ppars
  for(rel in prelspec) {
    if(length(rel$ind1) == 1) 
      pinde[[rel$what1]] <- append(pinde[[rel$what1]], rel$ind1)
    else {
      indf <- rel$ind1[2]
      if(rel$ind1[1] > 1) {
        slW<-slot(model, rel$what1)
        indf <- indf + sum(unlist(lapply(slW, length))[1:(rel$ind1[1]-1)]) 
      }
      pinde[[ rel$what1 ]] <-
        append(pinde[[ rel$what1 ]], indf)
    }
  }
  pinde
}

