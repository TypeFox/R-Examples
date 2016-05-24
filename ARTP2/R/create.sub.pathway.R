
create.sub.pathway <- function(super.pathway){
  
  sub.pathway <- list()
  grp <- sort(unique(super.pathway$Chr))
  ngrp <- length(grp)
  for(i in 1:ngrp){
    sub.pathway[[i]] <- super.pathway[super.pathway$Chr == grp[i], ]
    rownames(sub.pathway[[i]]) <- NULL
  }
  
  sub.pathway
  
}

