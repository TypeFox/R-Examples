#' Overview over the parameters used in the \code{randomizeR} package
#' 
#' This list of parameters yields a comprehensive overview of the parameters
#' used in the \code{randomizeR} package.
#'
#' @param add integer representing the number of balls that are added to the
#' urn in each step.
#' @param alpha the level of the t.test in each simulation.
#' @param bc vector which contains the lengths \code{k_1,...,k_l} of each block. 
#' This means that the vector \code{bc} will have one entry for each block.
#' @param b numeric vector of length at most 2 specifying the weight(s) for the punishment of 
# deviations from the target value.
#' @param ini integer representing the initial urn composition.
#' @param compr factor of compression for the sigmoid-time trend.
#' @param d effect size.
#' @param delta first noncentrality parameter of the doubly noncentral t-distribution.
#' @param df degrees of freedom (i.a. \code{N-2}).
#' @param eta numeric specifying the magnitude of selection bias.
#' @param file A connection, or a character string naming the file to write to.
#' @param filledBlock \code{logical} whether the last block should be filled or not. 
#' @param FTI final tolerated imbalance. This is the difference in number of
#' patients of groups A and B that is permitted at the end of a trial. Usually
#' this is set to zero.
#' @param groups character vector of labels for the different treatments.
#' @param k length of the block to be permuted. \code{k} should be divisible by
#' the number of treatment arms.
#' @param K number of treatment groups (e.g. K=2 if we compare one experimental
#' against one control treatment).
#' @param lb lower bound for the starting value of the poisson distribution.
#' @param lambda (second) noncentrality parameter of the doubly noncentral t-distribution.
#' @param method  method that is used to generate the (random) allocation
#' sequence. It can take values \code{PBR}, \code{RAR}, \code{HAD}, \code{PWR},
#' \code{EBC}, \code{BSD}, \code{CR}, \code{TBD}, \code{UD}, and \code{MP}.
#' @param mti maximum tolerated imbalance in patient numbers during the trial.
#' @param N integer for the total sample size of the trial.
#' @param name name of a variable.
#' @param mu vector of the expected responses of the treatment groups, should have
#' length \code{K} 
#' (i.e. one entry for each treatment group).
#' @param obj object specifying the randomization procedure, i.e. an object of a
#' class.
#' @param object any R object.
#' @param oject any R object.
#' Inheriting from \linkS4class{randPar}. See also \code{\link{createParam}}.
#' @param p success probability of the biased coin (e.g. in Efron's Biased Coin
#' Design).
#' @param pr vector with patient responses, i.e. each patients resulting value
#' after the treatment.
#' @param q "cut-off" value in \code{[0.5,1]}. This is the ratio of patients up 
#' from which the experimenter imposes selection bias on the data.
#' @param r numeric indicating the number of random sequences to be generated at
#' random or missing.
#' @param ratio vector of length \code{K}. The total sample number \code{N} and
#' all used block lengths (\code{bc}) have to be divisible by \code{sum(ratio)}.
#' @param rb block lengths of the blocks that can be selected equiprobable at random.
#' @param rsob randomization sequence (of one block).
#' @param rs randomization sequence (of all blocks).
#' @param S matrix for the computation of the probabilities in the maximal
#' procedure.
#' @param saltus integer or  \code{missing} specifying the patient index (i.e. position)
#' of the step in case of step time trend.
#' @param seed seed for the random number generation
#' @param sigma vector of the standard deviations in each the treatment group, 
#' should have length \code{K} (i.e. one entry for each treatment group).
#' @param SLs numeric vector of length at most 2 specifying the lower and/or upper 
#' specified border.
#' @param theta factor of the time trend for further details see \code{type}.
#' @param type character vector indicating which biasing strategy the
#' experimenter is using (selection bias) and which other bias is present in the
#' clinical trial (e.g. time trend). All biases included in the vector are
#' combined (i.e. added up) to form the total bias. Possible values are
#' \code{"none"} (if no bias occurs), \code{"CS"} (resp. \code{"DS"}) (if the
#' experimenter uses the convergence (resp. divergence) strategy to invoke
#' selection bias), \code{LinT} for linear time trend, \code{LogT} for
#' log-linear time trend, \code{StepT} for step time trend, \code{SigT} for
#' sigmoid time trend, \code{PWR} for knowledge of all up to the first
#' observation in each block, \code{MTI} the next observation after reaching the
#' maximal tolerated imbalance is reached will be known to the physican.
#' @param TV numeric specifying the optimal desired value called the target value.
#' @param varEq \code{logical} parameter for the t.test: Shall the variances be treated 
#' as equal (TRUE= t.test) or different (FALSE= Welch.test).
#' @param ub upper bound for the last value of the poisson distribution.
#' @param x a variable \code{x}.
#' @param rho nonnegative parameter which my be adjusted according to how strongly it is 
#' desired to balance the experiment. If \code{rho = 1}, we have Wei's urn design with
#' \code{alpha = 0}. If \code{rho = 0}, we have complete randomization. 
#' @param a nonnegative parameter which my be adjusted according to how strongly it is 
#' desired to balance the experiment. a = 0 gives the complete randomization, while the
#' assignments become more deterministic as \code{a} increases.
#' @param a nonnegative parameter which controls the degree of randomness: 
#' For decreasing \code{a} the allocations become deterministic, while for increasing
#' \code{a} this procedure tends to complete randomization. 
#'
#' @name overview
NULL

