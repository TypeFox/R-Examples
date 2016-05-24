#' Hybrid model simulation.
#' 
#' @description  \code{hybridModel} function runs hybrid models simulations.
#' 
#' @param network a \code{\link{data.frame}} with variables that describe
#'        the donor node, the reciever node, the time when each connection between
#'        donor to the reciever happened and the weight of these connection.
#' 
#' @param var.names a \code{\link{list}} with variable names of the network:
#'        the donor node, the reciever node, the time when each connection between
#'        donor to the reciever happened and the weight of these connection.
#'        The variables names must be "from", "to", "Time" and "arc", respectively.
#' 
#' @param link.type a \code{\link{character}} describing the link type between 
#'        nodes. There are two types: 'migration' and 'influence'. In the migration
#'        link type there are actual migration between nodes. In the influence 
#'        link type individuals does not migrate, just influences another node.
#'        
#' @param model a \code{\link{character}} describing model's name.
#' 
#' @param init.cond a named \code{\link{vector}} with initial conditions.
#' 
#' @param fill.time It indicates whether to return all dates or just the dates
#'        when nodes get connected.
#'
#' @param model.parms a named \code{\link{vector}} with model's parameters.
#' 
#' @param prop.func a character \code{\link{vector}} with propensity functions
#'        of a generic node. See references for more details
#' 
#' @param state.change.matrix is a state-change \code{\link{matrix}}. See references
#'        for more details
#' 
#' @param state.var a character \code{\link{vector}} with the state varialbes of
#'        the propensity functions.
#' 
#' @param ssa.method a \code{\link{list}} with SSA parameters. The default method
#'        is the direct method. See references for more details
#'
#' @param nodesCensus a \code{\link{data.frame}} with the first column describing
#'        nodes' ID, the second column with the number of individuals and the third
#'        describing the day of the census.
#'
#' @param sim.number Number of repetitions.The default value is 1
#'
#' @param pop.correc Whether \code{hybridModel} function tries to balance the number
#'        of individuals or not. The default value is TRUE.
#'        
#' @param num.cores  number of  threads/cores that the simulation will use. the
#'        default value is num.cores = 'max', the algothim will use all
#'        threads/cores available.
#'
#' @return Object containing a \code{\link{data.frame}} (results) with the number
#'         of individuals through time per node and per state.
#'
#' @references
#' [1] Pineda-krch, M. (2008). GillespieSSA : Implementing the Stochastic
#'     Simulation Algorithm in R. Journal of Statistical Software, 25(12).
#'
#' @seealso \link{GillespieSSA}.
#' @export
#' @examples 
#' # Parameters and initial conditions for an SIS model
#' # loading the data set 
#' data(networkSample) # help("networkSample"), for more info
#' networkSample <- networkSample[which(networkSample$Dia < "2012-03-20"),]
#' 
#' var.names <- list(from = 'originID', to = 'destinyID', Time = 'Dia',
#'                   arc = 'num.animais')
#'                   
#' prop.func <- c('beta * S * I / (S + I)', 'gamma * I')
#' state.var <- c('S', 'I')
#' state.change.matrix <- matrix(c(-1,  1,  # S
#'                                  1, -1), # I
#'                               nrow = 2, ncol = 2, byrow = TRUE)
#'                               
#' model.parms <- c(beta = 0.1, gamma = 0.01)
#'
#' init.cond <- rep(100, length(unique(c(networkSample$originID,
#'                                       networkSample$destinyID))))
#' names(init.cond) <- paste('S', unique(c(networkSample$originID,
#'                                         networkSample$destinyID)), sep = '')
#' init.cond <- c(init.cond, c(I36811 = 10, I36812 = 10)) # adding infection
#'                   
#' # running simulations, check num of cores available (num.cores)
#' sim.results <- hybridModel(network = networkSample, var.names = var.names,
#'                            model.parms = model.parms, state.var = state.var,
#'                            prop.func = prop.func, init.cond = init.cond,
#'                            state.change.matrix = state.change.matrix,
#'                            sim.number = 2, num.cores = 2)
#' 
#' # default plot layout (plot.types: 'pop.mean', 'subpop', or 'subpop.mean')
#' plot(sim.results, plot.type = 'subpop.mean')
#' 
#' # changing plot layout with ggplot2 (example)
#' # uncomment the lines below to test new layout exemple
#' #library(ggplot2)
#' #plot(sim.results, plot.type = 'subpop') + ggtitle('New Layout') + 
#' #  theme_bw() + theme(axis.title = element_text(size = 14, face = "italic"))
#'
hybridModel <-   function(network, var.names, link.type = 'migration',
                          model = 'custom', init.cond, fill.time = F,
                          model.parms, prop.func = NULL, state.change.matrix = NULL, 
                          state.var = NULL,
                          ssa.method = list(method = "D", epsilon = 0.03,
                                            nc = 10, dtf = 10, nd = 100),
                          nodesCensus = NULL, sim.number = 1, pop.correc = TRUE,
                          num.cores = 'max'){
  
  
  #### Extracting, trasforming and loading the dynamic network #####
  network <- network[, c(var.names$from, var.names$to, var.names$Time, var.names$arc)]
  
  #### building classes ####
  if(model == 'SI model without demographics' & link.type == 'migration'){
    model1 <- 'siWoDemogrMigr'
  } else if(model == 'SI model without demographics' & link.type == 'influence'){
    model1 <- 'siWoDemogrInfl'
  } else if(model == 'custom' & link.type == 'migration'){
    model1 <- 'customMigr'
  }
  
  model2simulate <- buildModelClass(structure(list(network = network,
                                                   ssa.method = ssa.method,
                                                   pop.correc = pop.correc,
                                                   nodes.info = nodesCensus),
                                              class = c(model1, 'HM')), var.names,
                                    init.cond, model.parms, prop.func, state.var,
                                    state.change.matrix)
                                    
                            
  #### running the simulation ####
  model2simulate$results <- simHM(model2simulate, network, sim.number, num.cores, fill.time = F)

  return(model2simulate)
}