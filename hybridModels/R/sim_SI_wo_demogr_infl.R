#' @import foreach
#' 
simHM.siWoDemogrInfl <- function(x, network, sim.number, num.cores){
  
  siWoDemogrInfl <- function(){
    
    # making it readable
    sim.result <- x$results
    ssaObject <- x$ssaObjet
    from <- x$ssaObjet$var.names$from
    arc <- x$ssaObjet$var.names$arc
    to <- x$ssaObjet$var.names$to
    Time <- x$ssaObjet$var.names$Time
    
    # starting
    sim.result$sim <- sims
    
    for(tempo in 1:length(ssaObject[['mov.dates']])){
      
      ssaObject$parms <- x$ssaObjet$parms
      
      # Extracting the total number of migrants per node per tempo
      # emigrants
      emigrants <- 
        stats::aggregate(network[which(network[, Time] == ssaObject[['mov.dates']][tempo]),
                          c(from,arc)][, arc],
                  by = list(network[which(network[, Time] == ssaObject[['mov.dates']][tempo]),
                                    c(from,arc)][, from]),
                  FUN = sum)
      colnames(emigrants) <- c(from, arc)
      
      # imigrants
      imigrants <-
        stats::aggregate(network[which(network[, Time] == ssaObject[['mov.dates']][tempo]),
                          c(to,arc)][, arc],
                  by = list(network[which(network[, Time] == ssaObject[['mov.dates']][tempo]),
                                    c(to,arc)][, to]),
                  FUN = sum)
      colnames(imigrants) <- c(to, arc)
      
      # Creating a Vector with current infected/Suscetibles/Efetivo/Traded animals
      s.individuals <- as.vector(sim.result[which(sim.result[, Time] == ssaObject[['mov.dates']][tempo]),
                                            ssaObject$IDs$S.ID], mode = 'integer')
      names(s.individuals) <- ssaObject$IDs$nodes.ID
      i.individuals <- as.vector(sim.result[which(sim.result[, Time] == ssaObject[['mov.dates']][tempo]),
                                            ssaObject$IDs$I.ID], mode = 'integer')
      names(i.individuals) <- ssaObject$IDs$nodes.ID
      n.individuals <- as.vector(sim.result[which(sim.result[, Time] == ssaObject[['mov.dates']][tempo]),
                                            ssaObject$IDs$N.ID], mode = 'integer')
      names(n.individuals) <- ssaObject$IDs$nodes.ID
      
      
      # ---------------- Hypergeom Wiki--------------------
      # mean da hypergeom  <- n*(m/N}, N is poputiation size
      # m, number of successes states in the population
      # n, number of draws
      # --------------------------------------------------
      infected.emigrants <- apply(cbind(i.individuals[as.character(emigrants[, from])],
                                        s.individuals[as.character(emigrants[, from])],
                                        emigrants[, arc]),
                                  1, function(x) {stats::rhyper(1,x[[1]],x[[2]],x[[3]])})
      
      susceptible.emigrants <- emigrants[, arc] - infected.emigrants
      
      
      # -------------        Expensive Option        --------------
      # ------------- Randomly distributing infected --------------
      # connected.nodes is a vector data frame with the connected nodes in the time tempo
      connected.nodes <-
        stats::aggregate(network[which(network[, Time] == ssaObject[['mov.dates']][tempo]),
                          c(from,arc)][, arc],
                  by = list(network[which(network[, Time] == ssaObject[['mov.dates']][tempo]),
                                    c(from,arc)][, from],
                            network[which(network[, Time] == ssaObject[['mov.dates']][tempo]),
                                    c(to,arc)][, to]),
                  FUN = sum)
      colnames(connected.nodes) <- c(from, to, arc)
      
      # individual.emigrants is a list of traded individuals (susceptible = zero and infected = 1 in a random 
      # order in each node)
      individual.emigrants <- apply(cbind(rbind(susceptible.emigrants, infected.emigrants,
                                                emigrants = emigrants[, arc]), foo = 0), 2,
                                    function(x){sample(c(rep.int(0, x[1]), rep.int(1, x[2])), x[3],
                                                       replace=FALSE)})
      
      connected.nodes <- cbind.data.frame(connected.nodes, iin = integer(length(connected.nodes[, to])))
      
      # distribution of infected and susceptible
      for(donor.node in emigrants[, from]){
        
        first <- 1
        last <- 0
        
        for(reciever.node in which(connected.nodes[ , from] == donor.node)){
          
          last <- last + connected.nodes[reciever.node, arc]
          
          connected.nodes[reciever.node, 'iin'] <- sum(individual.emigrants[[as.character(donor.node)]][first:last])
          
          first <- first + connected.nodes[reciever.node, arc]
          
        }
      }
      connected.nodes <- cbind.data.frame(connected.nodes,
                                          sin = connected.nodes[, arc] - connected.nodes[, 'iin'])
      # ----------------------------------------------------------
      
      # first balance
      emigrants <- cbind.data.frame(emigrants, iout = sapply(individual.emigrants,sum)[-length(individual.emigrants)])
      emigrants <- cbind.data.frame(emigrants, sout = emigrants[, arc] - emigrants$iout)
      
      imigrants <- cbind.data.frame(imigrants, 
                                    stats::aggregate(connected.nodes[, c('iin', 'sin')],
                                              by = list(connected.nodes[, to]), sum))
      imigrants <- imigrants[, -which(names(imigrants) == 'Group.1')]
          
      # Alocating parameters      
      names(i.individuals) <- paste('I', names(i.individuals), sep = '')
      names(n.individuals) <- paste('N', names(n.individuals), sep = '')      
      
      ssaObject$x0 <- i.individuals[names(ssaObject$x0)]
      
      ssaObject$parms[paste("iin", imigrants[, to], sep = "")] <- imigrants$iin
      ssaObject$parms[paste("nin", imigrants[, to], sep = "")] <- imigrants[, arc]   
      ssaObject$parms[names(n.individuals)] <- n.individuals
      
      
      
      # checking whether will be another trade
      out.sim <- GillespieSSA::ssa(x0 = ssaObject$x0, a = ssaObject$propFunction, nu = ssaObject$sCMatrix,
                                   parms = ssaObject$parms, tf = ssaObject$time.diff[tempo],
                                   censusInterval = 0, verbose = FALSE, method = "D")$data
      
      out.sim <- out.sim[(length(out.sim[,1])-2), -1]
      
      sim.result[which(sim.result[, Time] == ssaObject[['mov.dates']][tempo]) + 1,
                 names(out.sim)] <- out.sim
      
      sim.result[which(sim.result[, Time] == ssaObject[['mov.dates']][tempo]) + 1,
                 paste('S', names(s.individuals), sep = '')] <-
        ssaObject$parms[gsub('I','N',names(out.sim))] - out.sim
      
      sim.result[which(sim.result[, Time] == ssaObject[['mov.dates']][tempo]) + 1,
                 names(n.individuals)] <- ssaObject$parms[names(n.individuals)]
    }
    
    return(sim.result)
  }
  
  if(num.cores == 'max') {
    num.cores <- parallel::detectCores() 
  }
  
  cl <- parallel::makeCluster(num.cores, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  sims <- NULL
  sim.result <- foreach(sims = 1:sim.number, .verbose=FALSE, .inorder=FALSE,
                        .packages = 'GillespieSSA') %dopar% (siWoDemogrInfl())
  
  parallel::stopCluster(cl)
  sim.result <- do.call(rbind.data.frame, sim.result)
  
  return(sim.result)
}