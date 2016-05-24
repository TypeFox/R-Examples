"getPrelTheta" <-
function (th, modellist, diffs, d, parorder) 
{
  model <- modellist[[d]]
  fixed <- model@fvecind
  removepar <- fixed[["prel"]]
  if(length(unlist(slot(model, "prel"))) - length(removepar) != 0) {
    if(diffs$what %in%  modellist[[d]]@positivepar)
      parapp <- log(unlist(slot(model, "prel"))[-removepar])
    else {
      if(length(modellist[[d]]@clinde[[diffs$what]]) > 0) 
        for(i in 1:length(modellist[[d]]@clinde[[diffs$what]]))  
          parapp <- log(parapp[i])
      if(length(modellist[[d]]@chinde[[diffs$what]]) > 0) 
        for(i in 1:length(modellist[[d]]@chinde[[diffs$what]]))  
          parapp <- log(parapp[i])
      
    }
    parapp <- parapp[-removepar] 
  }
  else 
    parapp <- vector() 
  if(length(parapp) != 0) 
    ind <- (length(th) + 1):(length(th) + length(parapp))
  else ind <- vector()
  parorder[[length(parorder)+1]] <- list(name="prel",
                                         ind=ind, dataset = d, rm=removepar)
  th <- append(th, parapp)
  list(th=th, parorder=parorder)
}

