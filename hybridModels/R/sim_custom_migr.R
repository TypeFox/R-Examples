#' @import foreach
#' 
simHM.customMigr <- function(x, network, sim.number, num.cores, fill.time){
  
  if (fill.time == F){
    parallelCustomMigr <- function(){
      
      # making it readable
      sim.result <- x$results
      ssaObject <- x$ssaObjet
      from <- x$ssaObjet$var.names$from
      arc <- x$ssaObjet$var.names$arc
      to <- x$ssaObjet$var.names$to
      Time <- x$ssaObjet$var.names$Time
      state.var <- x$ssaObjet$state.var
      
      # starting
      sim.result$sim <- sims
      
      for(tempo in 1:length(ssaObject[['mov.dates']])){
        
        # Extracting the total number of migrants per node per tempo
        # emigrants
        emigrants <- 
          stats::aggregate(network[which(network[, Time] == ssaObject[['mov.dates']][tempo]),
                            c(from,arc)][, arc],
                    by = list(network[which(network[, Time] == ssaObject[['mov.dates']][tempo]),
                                      c(from,arc)][, from]),
                    FUN = sum)
        colnames(emigrants) <- c(from, arc)
        
        ### sampling from nodes ###
        sampled <- apply(emigrants, 1,
                         function(x){
                           sampled <- sample(rep(state.var,
                                                 sim.result[tempo,
                                                            paste(state.var, x[1],
                                                                  sep ='')]),
                                             x[2], replace = F)})
        if (is.matrix(sampled) == T)
          sampled <- as.list(data.frame(sampled, stringsAsFactors = F))
        names(sampled) <- emigrants[,1]
        
        # -------------        Expensive Option        --------------
        # ------------- Randomly distributing infected --------------
        # connected.nodes is a data frame with the connected nodes in the time tempo
        connected.nodes <-
          stats::aggregate(network[which(network[, Time] == ssaObject[['mov.dates']][tempo]),
                            c(from,arc)][, arc],
                    by = list(network[which(network[, Time] == ssaObject[['mov.dates']][tempo]),
                                      c(from,arc)][, from],
                              network[which(network[, Time] == ssaObject[['mov.dates']][tempo]),
                                      c(to,arc)][, to]),
                    FUN = sum)
        colnames(connected.nodes) <- c(from, to, arc)
        connected.nodes[,state.var] <- 0
        
        
        # distribution of individuals
        for(donor.node in emigrants[, from]){
          
          first <- 1
          last <- 0
          
          for(reciever.node in which(connected.nodes[ , from] == donor.node)){
            
            last <- last + connected.nodes[reciever.node, arc]          
            
            connected.nodes[reciever.node, state.var] <- 
              apply(as.matrix(state.var), 1,
                    function(x){
                      length(which(sampled[as.character(donor.node)][[1]][first:last] == x))})
            
            first <- first + connected.nodes[reciever.node, arc]
            
          }
        }
        
        # balancing
        connected.emigrants <- stats::aggregate(connected.nodes[, state.var], by = list(connected.nodes[, from]), sum)
        connected.imigrants <- stats::aggregate(connected.nodes[, state.var], by = list(connected.nodes[, to]), sum)
        
        ssaObject$x0[as.vector(apply(as.matrix(state.var), 1, function(x)
          paste(x, connected.emigrants[,'Group.1'], sep = '')))] <- as.vector(t(apply(connected.emigrants, 1, function(x){
            ssaObject$x0[paste(state.var, x['Group.1'], sep = '')] - as.numeric(x[state.var])})))
        
        ssaObject$x0[as.vector(apply(as.matrix(state.var), 1, function(x)
          paste(x, connected.imigrants[,'Group.1'], sep = '')))] <- as.vector(t(apply(connected.imigrants, 1, function(x){
            ssaObject$x0[paste(state.var, x['Group.1'], sep = '')] + as.numeric(x[state.var])})))
        
        # checking whether will be another trade
        out.sim <- GillespieSSA::ssa(x0 = ssaObject$x0, a = ssaObject$propFunction, nu = ssaObject$sCMatrix,
                                     parms = ssaObject$parms, tf = ssaObject$time.diff[tempo],
                                     censusInterval = 0, verbose = FALSE, method = ssaObject$ssa.method$method,
                                     epsilon = ssaObject$ssa.method$epsilon, nc = ssaObject$ssa.method$nc,
                                     dtf = ssaObject$ssa.method$dtf, nd = ssaObject$ssa.method$nd)$data
        
        out.sim <- out.sim[(length(out.sim[,1])-2), -1]
        
        sim.result[which(sim.result[, Time] == ssaObject[['mov.dates']][tempo]) + 1,
                   names(out.sim)] <- out.sim
        
        ssaObject$x0[names(out.sim)] <- out.sim
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
                          .packages = 'GillespieSSA') %dopar% (parallelCustomMigr())
    
    parallel::stopCluster(cl)
    
    sim.result <- do.call(rbind.data.frame, sim.result)
    
    return(sim.result) 
  }
}