#' Overview over the parameters used in the \code{unitedR} package
#' 
#' This list of parameters yields a comprehensive overview of the parameters
#' used in the \code{unitedR} package.
#'
#' @param away away team (an object of the \code{S4}class \code{formation})
#' @param DF \code{numeric} vector for the strengths of the players in the 
#' defense
#' @param formation object of the \code{S4}class \code{formation}
#' @param GK \code{integer} for the strength goalkeeper
#' @param hardness \code{numeric} vector of length five with integers for the used hardness
#' @param home home team (an object of the \code{S4}class \code{formation})
#' @param homeAdv \code{numeric} vector of length five with integers for the used hardness
#' @param MF \code{numeric} vector for the strengths of the players in the 
#' midfield
#' @param penaltyGoalProb probability of a goal by a singular penalty
#' @param posPenalties number of possible penalties in a game
#' @param preventGoalGK factor multiplicied with the strength of the GK for computing the 
#' probability of preventing a goal by the goalkeeper
#' @param preventGoalSW factor multiplicied with the strength of the SW for computing the 
#' probability of preventing a goal by the sweeper
#' @param r number of replications for the simulation of hardness and penalties, can 
#' be \code{missing} (exact results will be computed)
#' @param ST \code{numeric} vector of integers for the strenghts of the strikers
#' @param SW \code{vector} for the strength of the sweeper, can be 
#' \code{NA} or a \code{numeric}
#' @param x a variable \code{x}.
#'
#' @name overview
NULL