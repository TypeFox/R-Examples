
set.modules <- function(modules,DIC){

  #Load/unload appropriate modules (besides dic)
  called.set <- c('basemod','bugs',modules)
  current.set <- list.modules()

  load.set <- called.set[!called.set%in%current.set]
  unload.set <- current.set[!current.set%in%called.set]

  if(length(load.set)>0){
    for (i in 1:length(load.set)){
      load.module(load.set[i],quiet=TRUE)
    }
  }
  if(length(unload.set)>0){
    for (i in 1:length(unload.set)){
      unload.module(unload.set[i],quiet=TRUE)
    }
  }
  if(DIC){
    load.module("dic",quiet=TRUE)
  }
}