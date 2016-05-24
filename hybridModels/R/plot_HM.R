#' @name plot
#' @title Plot for SI without demographics model
#' 
#' @description  \code{plot.HM} is a method to plot hybrid models from this
#'               package
#'               
#' 
#' @param x \code{HM} object
#' 
#' @param sim points to which simulation to plot.
#' 
#' @param facet.scales should scales be fixed ("free_y", the default), free ("free"), or free
#'        in one dimension ("free_x", "free_y"). See ggplot2 package for more
#'        details.
#'        
#' @param plot.type plots the mean number of each state variable for the whole
#'        population ('pop.mean'), or the subpopulations of a particular
#'        simulation ('subpop'), or the mean of each subpopulation ('subpop.mean').
#'        
#' @param ... arguments to be passed to methods.
#' 
#' @export
#' 
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
plot.HM <- function(x, sim = 1, plot.type = 'subpop', facet.scales = 'free_y', ...){
  
  plot.type <- match.arg(plot.type, c('subpop', 'pop.mean', 'subpop.mean'))
  
  Time <- Number <- variable <- State <- Subpop <- NULL
  
  if(plot.type == 'subpop'){
    sim.result <- x$results[which(x$results$sim == sim), ]
    
    sim.result.plot <- reshape2::melt(sim.result, id.vars = c('sim',
                                                              x$ssaObjet$var.names$Time))
    sim.result.plot[, 'State'] <- substring(sim.result.plot$variable, 1, 1)
    colnames(sim.result.plot)[c(2,4)] <- c('Time','Number')
    
    sim.result.plot$State <- factor(sim.result.plot$State, levels = x$ssaObjet$state.var)
    
    return(ggplot2::ggplot(sim.result.plot, ggplot2::aes(x = Time, y = Number,
                                                         group = variable, color = State)) + 
             ggplot2::geom_line(alpha = 0.4, size = 0.3) + ggplot2::ggtitle(paste('Simulation', sim)) +
             ggplot2::ylab('Number Of Individuals') + ggplot2::guides(color=FALSE) +
             ggplot2::facet_wrap(~State, ncol = 1, scales = facet.scales) +      
             ggplot2::theme(panel.border = ggplot2::element_rect(colour = "grey", fill = NA),
                            strip.background = ggplot2::element_rect(fill = NA, colour = "grey", size = 0.1),
                            strip.text = ggplot2::element_text(face = "bold", size = 12),
                            panel.margin = grid::unit(0.6, "lines"),
                            plot.title = ggplot2::element_text(size = 14, face = "bold"),
                            axis.text = ggplot2::element_text(size = 10),
                            axis.title = ggplot2::element_text(size = 12, face = "bold")))
  } else if(plot.type == 'pop.mean'){
    sim.result <- x$results
    sim.result.plot <- reshape2::melt(sim.result, id.vars = c('sim',
                                                              x$ssaObjet$var.names$Time))
    sim.result.plot[, 'State'] <- substring(sim.result.plot$variable, 1, 1)
    sim.result.plot <- stats::aggregate(sim.result.plot$value,
                                 by = list(State = sim.result.plot$State,
                                           Sim = sim.result.plot$sim,
                                           Time = sim.result.plot[ , x$ssaObjet$var.names$Time]), sum)
    sim.result.plot <- stats::aggregate(sim.result.plot$x,
                                 by = list(State = sim.result.plot$State,
                                           Time = sim.result.plot$Time), mean)
    sim.result.plot$State <- factor(sim.result.plot$State, levels = x$ssaObjet$state.var)
    
    return(ggplot2::ggplot(sim.result.plot, ggplot2::aes(x = Time, y = x,
                                                         group = State, color = State)) + 
             ggplot2::geom_line(size = 0.3) + ggplot2::ggtitle('Population') +
             ggplot2::ylab('Mean Number Of Individuals') + ggplot2::guides(color=FALSE) +
             ggplot2::facet_wrap(~State, ncol = 1, scales = facet.scales) +      
             ggplot2::theme(panel.border = ggplot2::element_rect(colour = "grey", fill = NA),
                            strip.background = ggplot2::element_rect(fill = NA, colour = "grey", size = 0.1),
                            strip.text = ggplot2::element_text(face = "bold", size = 12),
                            panel.margin = grid::unit(0.6, "lines"),
                            plot.title = ggplot2::element_text(size = 14, face = "bold"),
                            axis.text = ggplot2::element_text(size = 10),
                            axis.title = ggplot2::element_text(size = 12, face = "bold")))
  } else if(plot.type == 'subpop.mean'){
    sim.result <- x$results
    sim.result.plot <- reshape2::melt(sim.result, id.vars = c('sim',
                                                              x$ssaObjet$var.names$Time))
    sim.result.plot <- stats::aggregate(sim.result.plot$value,
                                 by = list(Subpop = sim.result.plot$variable,
                                           Time = sim.result.plot[ , x$ssaObjet$var.names$Time]), mean)
    
    sim.result.plot[, 'State'] <- substring(sim.result.plot$Subpop, 1, 1)
    sim.result.plot$State <- factor(sim.result.plot$State, levels = x$ssaObjet$state.var)
    
    return(ggplot2::ggplot(sim.result.plot, ggplot2::aes(x = Time, y = x,
                                                         group = Subpop, color = State)) + 
             ggplot2::geom_line(alpha = 0.4, size = 0.3) + ggplot2::ggtitle('Subpopulations') +
             ggplot2::ylab('Mean Number Of Individuals') + ggplot2::guides(color=FALSE) +
             ggplot2::facet_wrap(~State, ncol = 1, scales = facet.scales) +      
             ggplot2::theme(panel.border = ggplot2::element_rect(colour = "grey", fill = NA),
                            strip.background = ggplot2::element_rect(fill = NA, colour = "grey", size = 0.1),
                            strip.text = ggplot2::element_text(face = "bold", size = 12),
                            panel.margin = grid::unit(0.6, "lines"),
                            plot.title = ggplot2::element_text(size = 14, face = "bold"),
                            axis.text = ggplot2::element_text(size = 10),
                            axis.title = ggplot2::element_text(size = 12, face = "bold")))
  }
}