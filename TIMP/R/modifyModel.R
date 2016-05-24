"modifyModel" <- function(model = list(), newest = list(), 
   exceptslots = vector()){
     if(length(newest) == 0) 
       theta <- get(".currTheta", envir = .GlobalEnv)[[1]]
     else 
        theta <- newest$currTheta[[1]]
     if(length(model) == 0)
        m <- get(".currModel", envir = .GlobalEnv)@model
     else 
        m <- model
     if(length(exceptslots) > 0) 
        nm <- slotNames(theta)[- which(slotNames(theta) %in% exceptslots)]
     else nm <- slotNames(theta)

     nm <- intersect(nm, slotNames(m)) 

     arglist <- as.list(m@datCall[[length(m@datCall)]])
     aa <- list()
     for(i in 2:length(arglist)) 
                aa[[ names(arglist)[i] ]] <- arglist[[i]]
     for (i in 1:length(nm)) 
           aa[[nm[i]]] <- slot(theta, nm[i] )
  
      do.call("initModel", aa)
}
