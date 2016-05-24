#' @title CUB package
#' @description 
#' The analysis of human perceptions is often carried out by resorting to questionnaires, 
#' where respondents are asked to express ratings about the items being evaluated. The standard goal of the 
#' statistical framework proposed for this kind of data (e.g. cumulative models) is to explicitly characterize
#' the respondents' perceptions about a latent trait, by taking into account, at the same time, 
#' the ordinal categorical scale of measurement of the   involved statistical variables. The new class of models 
#' starts from a particular assumption about the unconscious mechanism leading individuals' responses  
#'  to choose an ordinal category on a rating scale. The basic idea derives from the awareness that two latent
#'   components move the psychological process of   selection among discrete alternatives: attractiveness
#'    towards the item and uncertainty in the response. Both components of models concern the stochastic 
#'  mechanism in term of feeling, which is an internal/personal movement of the subject towards the item,
#'   and uncertainty pertaining to the final choice.\cr
#'   Thus, on the basis of experimental data and statistical motivations, the response distribution is modelled 
#'   as the convex Combination of a discrete Uniform and a shifted Binomial random variable (denoted as CUB models) 
#'   whose parameters may be consistently estimated and validated by maximum likelihood inference. 
#'   In addition, subjects' and objects' covariates can be included in the model in order to assess how the 
#'   characteristics of the respondents may affect the ordinal score. \cr
#'   CUB models have been firstly introduced by Piccolo (2003) and implemented on real datasets by D'Elia 
#'   and Piccolo (2005), Iannario and Piccolo (2012).\cr
#'   The CUB package allows the user to estimate and test CUB models and their extensions by using maximum 
#'   likelihood methods.\cr
#'  ACKNOWLEDGEMENTS: The Authors are grateful to Maria Antonietta Del Ferraro, Francesco Miranda and
#'   Giuseppe Porpora for their support in the implementation of the package.
#' @details 
#'   \tabular{ll}{
#' Package: \tab CUB\cr
#' Type: \tab Package\cr
#' Version: \tab 0.0\cr
#' Date: \tab 2015-10-30\cr
#' License: GPL-2 | GPL-3
#'  }
#' @source  \url{http://www.labstat.it/home/research/resourses/cub-data-sets-2/}
#' @author  Maria Iannario, Domenico Piccolo 
#' @references   Piccolo D. (2003). On the moments of a mixture of uniform and shifted binomial random variables,
#'  \emph{Quaderni di Statistica}, \bold{5}, 85--104 \cr  D'Elia A. and Piccolo D. (2005). 
#'  A mixture model for preferences data analysis, \emph{Computational Statistics & Data Analysis}, 
#'  \bold{49}, 917--937 \cr
#'  Iannario M. and Piccolo D. (2012). CUB models: Statistical methods and empirical evidence, in: 
#'  Kenett  R. S. and Salini S. (eds.), \emph{Modern Analysis of Customer Surveys: with applications using R},
#'   J. Wiley and   Sons, Chichester, 231--258\cr
#'   Iannario M. and Piccolo D. (2014). Inference for CUB models: a program in R, 
#'   \emph{Statistica & Applicazioni}, \bold{XII} n.2, 177--204
#' @name CUB_package 
#' @keywords package
NULL
#> NULL