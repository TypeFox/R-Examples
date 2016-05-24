
recreate.pathway <- function(setup, pathway){
  
  setup$super.pathway <- setup$pathway
  npath <- length(pathway)
  path <- list()
  for(i in 1:npath){
    path[[i]] <- setup$super.pathway[setup$super.pathway$Gene %in% pathway[[i]][, 'Gene'], ]
  }
  setup$pathway <- path
  setup
  
}

