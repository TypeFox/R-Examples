#' @name Hyde-package
#' @aliases Hyde
#' 
#' @title Hydbrid Decision Networks
#' @description Facilities for easy implementation of hybrid Bayesian networks 
#' using R. Bayesian networks are directed acyclic graphs representing joint 
#' probability distributions, where each node represents a random variable and 
#' each edge represents conditionality. The full joint distribution is 
#' therefore factorized as a product of conditional densities, where each node 
#' is assumed to be independent of its non-desendents given information on its 
#' parent nodes. Since exact, closed-form algorithms are computationally 
#' burdensome for inference within hybrid networks that contain a combination 
#' of continuous and discrete nodes, particle-based approximation techniques 
#' like Markov Chain Monte Carlo are popular. We provide a user-friendly 
#' interface to constructing these networks and running inference using rjags. 
#' Econometric analyses (maximum expected utility under competing policies, 
#' value of information) involving decision and utility nodes are also 
#' supported.
#'   

NULL