"getTheta" <- function (mod) 
{
  modelspec <- mod@modelspec
  modellist <- mod@modellist
  datasetind <- mod@datasetind 
  th <- vector()
  parorder <- list()
  for(i in 1:length(modelspec)) {
    ds <- which(datasetind == i)
    model <- modellist[[ds[1]]]
    fixed <- model@fvecind
    prel <- model@pvecind
    mrel <- model@mvecind 
    ppars <- append(model@parnames, "prel")
    for (p in ppars) {
      if(length(slot(modelspec[[datasetind[i]]], p)) != 0) {  
        removepar <- sort(unique(c(fixed[[p]], prel[[p]])))
        rmpar <- sort(unique(c(fixed[[p]], prel[[p]], mrel[[p]])))
        parapp <- unlist(slot(model, p))
        if(length(model@clinde[[p]]) > 0) {
          for(j in 1:length(model@clinde[[p]]))  {
            ind <- model@clinde[[p]][j]
            parapp[ind] <- model@lowcon[[p]][j] - parapp[ind]
            if(! p %in% model@positivepar)
              parapp[ind] <- log(parapp[ind])
          }
        }
        if(length(model@chinde[[p]]) > 0) { 
          for(j in 1:length(model@chinde[[p]]))  {
            ind <- model@chinde[[p]][j]
            parapp[ind] <- parapp[ind] + model@highcon[[p]][j]
            if(! p %in% model@positivepar)
              parapp[ind] <- log(parapp[ind])
          }
        }
        if (length(rmpar) != 0)
          parapp <- parapp[-rmpar]
        if(length(parapp) != 0) 
          ind <- (length(th) + 1):(length(th) + length(parapp))
        else 
          ind <- vector()
        parorder[[length(parorder)+1]] <- list(name=p, ind=ind, 
                                               dataset = ds, rm=rmpar)
        if(p %in% model@positivepar && length(parapp) != 0) 
          parapp <- log(parapp)
        th <- append(th, parapp)
      }
    }  
  }
  mod@parorder <- processOrder(parorder, mod) 
  df <- getDiffTheta(th, mod)
  list(theta = df$theta, mod = df$mod)
}
    
