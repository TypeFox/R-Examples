assg2xBalance <- function(assg.obj, data, id.var, bal.vars, to.report = "all"){
  
  ## require("RItools")
  
  tr.vec <- rep(NA, nrow(data))
  fff <- formula(paste("Tr ~ ", paste(bal.vars, collapse = "+")))
  xbal.list <- list()
  n.groups <- length(assg.obj$assg)
  
  for(i in 1:n.groups){
    assg.gp <- assg.obj$assg[[i]]
    
    tr.idx <- unfactor(assg.gp[, 1])
    co.idx <- unfactor(assg.gp[, 2])
    
    tr.vec[data[[id.var]] %in% tr.idx] <- 1
    tr.vec[data[[id.var]] %in% co.idx] <- 0
    wh.gp <- data[[id.var]] %in% c(tr.idx, co.idx)
    
    data.tr.gp <- data.frame(cbind(data[wh.gp, ], Tr=tr.vec[wh.gp]))
    
    xbal.out <- RItools::xBalance(fff, data = data.tr.gp, report = c(to.report))
    xbal.list[[paste("Group", i, sep="")]] <- xbal.out
  }
  
  data.tr <- data.frame(cbind(data, Tr = tr.vec))
  data.tr <- data.tr[!(is.na(data.tr$Tr)), ]
  xbal.list[["Overall"]] <- RItools::xBalance(fff, data = data.tr, report = c(to.report))
  return(xbal.list)	
}