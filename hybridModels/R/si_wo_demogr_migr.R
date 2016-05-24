buildModelClass.siWoDemogrMigr <- function(x, var.names, init.cond, model.parms,
                                           prop.func, state.var, state.change.matrix){
  
  #### args creation ####
  nodes.ID <- sort(unique(c(x$network[, var.names$from],x$network[, var.names$to])))
  nodes.info <- x$nodes.info[which(x$nodes.info[, 1] %in% nodes.ID), ]
  mov.dates <- sort(unique(x$network[, var.names$Time]))
  time.diff <- c(as.numeric(mov.dates[2:length(mov.dates)] - mov.dates[1:(length(mov.dates)-1)]),1)
  number.nodes <- length(nodes.ID)
  
  #### first model ####
  
  # building vectors
  propFunc  <-  paste("Beta * (N", nodes.ID, " - I", nodes.ID, ") * I",
                      nodes.ID, " / N", nodes.ID, sep = "")
  x0 <- vector(mode='integer', number.nodes)
  model.parms <- c(model.parms, vector(mode='numeric', number.nodes))
  scM <- diag(1, nrow = number.nodes, ncol = number.nodes)
  
  # naming vectors
  nomes <- paste("I", nodes.ID, sep = "")
  
  nomes.results <- c(paste("S", nodes.ID, sep = ""),
                     paste("I", nodes.ID, sep = ""),
                     paste("N", nodes.ID, sep = ""))
  names(x0) <- nomes
  names(model.parms)[2:(number.nodes + 1)] <- paste("N", nodes.ID, sep = "")
  
  results <- as.data.frame(stats::setNames(replicate(length(nomes.results),
                                                     integer(length(mov.dates)), 
                                                     simplify = FALSE), nomes.results))
  
  results <- cbind.data.frame(sim = integer(length(mov.dates)), Time = mov.dates, results)
  colnames(results)[2] <- var.names$Time
  
  ### intial contidtions ###
  nodes.info[, 'S.ID'] <- paste("S", nodes.info[, 1], sep = '')
  nodes.info[, 'I.ID'] <- paste("I", nodes.info[, 1], sep = '')
  nodes.info[, 'N.ID'] <- paste("N", nodes.info[, 1], sep = '')
  # Infected
  results[1, names(init.cond)] <- init.cond
  # N
  results[1, nodes.info$N.ID]  <- nodes.info[, 2]
  # Suscptible
  results[1, nodes.info$S.ID]  <- (results[1, nodes.info$N.ID] - results[1, nodes.info$I.ID])
  
  #### adding last day ####
  results <- rbind(results, results[length(mov.dates),])
  results[(length(mov.dates) + 1), var.names$Time] <- mov.dates[length(mov.dates)] + 1
  
  
  return(structure(list(ssaObjet = list(propFunction = propFunc, x0 = x0, sCMatrix = scM,
                                        parms = model.parms, mov.dates = mov.dates,
                                        time.diff = time.diff, number.nodes = number.nodes,
                                        var.names = var.names, ssa.method =  x$ssa.method,
                                        IDs = nodes.info[, c(1,3:5)], pop.correc = x$pop.correc),
                        results = results), class = c('siWoDemogrMigr', 'HM')))
}