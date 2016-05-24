#' @name RRTCS-package
#' @aliases RRTCS-package
#' @aliases RRTCS
#' @docType package
#' @title Randomized Response Techniques for Complex Surveys
#' 
#' @description The aim of this package is to calculate point and interval estimation for linear parameters with data obtained from randomized response surveys. 
#' Twenty one RR methods are implemented for complex surveys:
#' 
#' - Randomized response procedures to estimate parameters of a qualitative stigmatizing characteristic: Christofides model, Devore model, Forced-Response model, 
#' Horvitz model, Horvitz model with unknown B, Kuk model, Mangat model, Mangat model with unknown B, Mangat-Singh model, Mangat-Singh-Singh model, 
#' Mangat-Singh-Singh model with unknown B, Singh-Joarder model, SoberanisCruz model and Warner model.
#' 
#' - Randomized response procedures to estimate parameters of a quantitative stigmatizing characteristic: BarLev model, Chaudhuri-Christofides model, 
#' Diana-Perri-1 model, Diana-Perri-2 model, Eichhorn-Hayre model, Eriksson model and Saha model.
#' 
#' Using the usual notation in survey sampling, we consider a finite population \eqn{U=\{1,\ldots,i,\ldots,N\}}, consisting of \eqn{N} different elements. 
#' Let \eqn{y_i} be the value of the sensitive aspect under study for the \eqn{i}th population element. Our aim is to estimate the finite population total 
#' \eqn{Y=\sum_{i=1}^N y_i} of the variable of interest \eqn{y} or the population mean \eqn{\bar{Y}=\frac{1}{N}\sum_{i=1}^N y_i}. If we can estimate the 
#' proportion of the population presenting a certain stigmatized behaviour \eqn{A}, the variable \eqn{y_i} takes the value 1 if \eqn{i\in G_A} 
#' (the group with the stigmatized behaviour) and the value zero otherwise. Some qualitative models use an innocuous or related attribute \eqn{B} whose 
#' population proportion can be known or unknown.
#'  
#' Assume that a sample \eqn{s} is chosen according to a general design \eqn{p} with inclusion probabilities \eqn{\pi_i=\sum_{s\ni i}p(s),i\in U}.
#'  
#' In order to include a wide variety of RR procedures, we consider the unified approach given by Arnab (1994). The interviews of individuals in the sample \eqn{s}
#' are conducted in accordance with the RR model. For each \eqn{i\in s} the RR induces a random response \eqn{z_i} (denoted scrambled response) so that the revised 
#' randomized response \eqn{r_i} (Chaudhuri and Christofides, 2013) is an unbiased estimation of \eqn{y_i}. Then, an unbiased estimator for the population total of 
#' the sensitive characteristic \eqn{y} is given by
#' \deqn{\widehat{Y}_R=\sum_{i\in s}\frac{r_i}{\pi_i}}
#' The variance of this estimator is given by:
#' \deqn{V(\widehat{Y}_R)=\sum_{i\in U}\frac{V_R(r_i)}{\pi_i}+V_{HT}(r)} 
#' where \eqn{V_R(r_i)} is the variance of \eqn{r_i} under the randomized device and \eqn{V_{HT}(r)} is the design-variance of the Horvitz Thompson estimator 
#' of \eqn{r_i} values.
#'  
#' This variance is estimated by:
#' \deqn{\widehat{V}(\widehat{Y}_R)=\sum_{i\in s}\frac{\widehat{V}_R(r_i)}{\pi_i}+\widehat{V}(r)}
#' where \eqn{\widehat{V}_R(r_i)} varies with the RR device and the estimation of the design-variance, \eqn{\widehat{V}(r)}, is obtained using Deville's method 
#' (Deville, 1993).
#'  
#' The confidence interval at \eqn{(1-\alpha)} \% level is given by
#' \deqn{ci=\left(\widehat{Y}_R-z_{1-\frac{\alpha}{2}}\sqrt{\widehat{V}(\widehat{Y}_R)},\widehat{Y}_R+z_{1-\frac{\alpha}{2}}\sqrt{\widehat{V}(\widehat{Y}_R)}\right)}
#' where \eqn{z_{1-\frac{\alpha}{2}}} denotes the \eqn{(1-\alpha)} \% quantile of a standard normal distribution.
#'  
#' Similarly, an unbiased estimator for the population mean \eqn{\bar{Y}} is given by
#' \deqn{\widehat{\bar{Y}}_R= \frac{1}{N}\sum_{i\in s}\frac{r_i}{\pi_i}}
#' and an unbiased estimator for its variance is calculated as:
#' \deqn{\widehat{V}(\widehat{\bar{Y}}_R)=\frac{1}{N^2}\left(\sum_{i\in s}\frac{\widehat{V}_R(r_i)}{\pi_i}+\widehat{V}(r)\right)} 
#' In cases where the population size \eqn{N} is unknown, we consider Hàjek-type estimators for the mean:
#' \deqn{\widehat{\bar{Y}}_{RH}=\frac{\sum_{i\in s}r_i}{\sum_{i\in s}\frac{1}{\pi_i}}}
#' and Taylor-series linearization variance estimation of the ratio (Wolter, 2007) is used.
#'  
#' In qualitative models, the values \eqn{r_i} and \eqn{\widehat{V}_R(r_i)} for \eqn{i\in s} are described in each model.
#'  
#' In some quantitative models, the values \eqn{r_i} and \eqn{\widehat{V}_R(r_i)} for \eqn{i\in s} are calculated in a general form (Arcos et al, 2015) as follows:
#'  
#' The randomized response given by the person \eqn{i} is 
#' \deqn{z_i=\left\{\begin{array}{lccc}
#' y_i & \textrm{with probability } p_1\\
#' y_iS_1+S_2 & \textrm{with probability } p_2\\
#' S_3 & \textrm{with probability } p_3
#' \end{array}
#' \right.}
#' with \eqn{p_1+p_2+p_3=1} and where \eqn{S_1,S_2} and \eqn{S_3} are scramble variables whose distributions are assumed to be known. We denote by \eqn{\mu_i} 
#' and \eqn{\sigma_i} respectively the mean and standard deviation of the variable \eqn{S_i,(i=1,2,3)}.
#'  
#' The transformed variable is
#' \deqn{r_i=\frac{z_i-p_2\mu_2-p_3\mu_3}{p_1+p_2\mu_1},} 
#' its variance is
#' \deqn{V_R(r_i)=\frac{1}{(p_1+p_2\mu_1)^2}(y_i^2A+y_iB+C)}
#' where 
#' \deqn{A=p_1(1-p_1)+\sigma_1^2p_2+\mu_1^2p_2-\mu_1^2p_2^2-2p_1p_2\mu_1}
#' \deqn{B=2p_2\mu_1\mu_2-2\mu_1\mu_2p_2^2-2p_1p_2\mu_2-2\mu_3p_1p_3-2\mu_1\mu_3p_2p_3}
#' \deqn{C=(\sigma_2^2+\mu_2^2)p_2+(\sigma_3^2+\mu_3^2)p_3-(\mu_2p_2+\mu_3p_3)^2}
#' and the estimated variance is
#' \deqn{\widehat{V}_R(r_i)=\frac{1}{(p_1+p_2\mu_1)^2}(r_i^2A+r_iB+C).} 
#' Some of the quantitative techniques considered can be viewed as particular cases of the above described procedure. Other models are described in the respective function.   
#' 
#' Alternatively, the variance can be estimated using certain resampling methods.  
#'    
#'    
#' @author Beatriz Cobo Rodríguez, Department of Statistics and Operations Research. University of Granada \email{beacr@@ugr.es}
#' 
#'  María del Mar Rueda García, Department of Statistics and Operations Research. University of Granada \email{mrueda@@ugr.es} 
#'  
#'  Antonio Arcos Cebrián, Department of Statistics and Operations Research. University of Granada \email{arcos@@ugr.es}
#'  
#'  Maintainer: Beatriz Cobo Rodríguez \email{beacr@@ugr.es}
#'
#' @references Arcos, A., Rueda, M., Singh, S. (2015). 
#' \emph{A generalized approach to randomised response for quantitative variables.}
#'  Quality and Quantity 49, 1239-1256.
#'
#' @references Arnab, R. (1994). 
#' \emph{Non-negative variance estimator in randomized response surveys.}
#'  Comm. Stat. Theo. Math. 23, 1743-1752.
#'  
#' @references Chaudhuri, A., Christofides, T.C. (2013).
#' \emph{Indirect Questioning in Sample Surveys}
#'  Springer-Verlag Berlin Heidelberg.
#'  
#' @references Deville, J.C. (1993).
#' \emph{Estimation de la variance pour les enquêtes en deux phases.}
#'  Manuscript, INSEE, Paris.
#'
#' @references Wolter, K.M. (2007).
#' \emph{Introduction to Variance Estimation.}
#'  2nd Edition. Springer.
#'  
#' @keywords package
NULL