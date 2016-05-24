"addPrel" <-
function (model) 
{
	## prel has structure 
	## list(list(what1, ind1, 
	##           what2, ind2, 
	##           rel, start), ...) 
	
  prelspec <- model@prelspec
  model@prel <- vector()
  for(diffs in prelspec){
    model@prel <- append(model@prel, diffs$start)  
  }
  for(diffs in prelspec){
    #model@prel <- append(model@prel, diffs$start)  
    if(length(diffs$rel) == 0 || diffs$rel == "lin"){
      newpar <- multiLin(model, diffs, diffs$start[1]) + diffs$start[2]
      if(length(diffs$ind1)==1)
        slot(model, diffs$what1)[diffs$ind1] <- newpar 
      if(length(diffs$ind1)==2) 
        slot(model, diffs$what1)[[diffs$ind1[1]]][diffs$ind1[2]] <- newpar      
    }
    else {
      if(diffs$rel == "multilin"){
        newpar <- diffs$start[[1]] + multiLin(model, diffs, diffs$start[2:length(diffs$start)]) # was diffs$start[1]
        if(length(diffs$ind1)==1)
          slot(model, diffs$what1)[diffs$ind1] <- newpar 
        if(length(diffs$ind1)==2) 
                slot(model, diffs$what1)[[diffs$ind1[1]]][diffs$ind1[2]] <- newpar 
      }
    }
  }
  model
}


