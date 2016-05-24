#' @name HorvitzUB
#' @aliases HorvitzUB
#' @title Horvitz-UB model
#' 
#' @description Computes the randomized response estimation, its variance estimation and its confidence interval through the Horvitz model (Horvitz et al., 1967, and Greenberg et al., 1969)
#' when the proportion of people bearing the innocuous attribute is unknown. 
#' The function can also return the transformed variable. 
#' The Horvitz-UB model can be seen in Chaudhuri (2011, page 42). 
#' 
#' @usage HorvitzUB(I,J,p1,p2,pi,type=c("total","mean"),cl,N=NULL,pij=NULL)
#' @param I first vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param J second vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param p1 proportion of marked cards with the sensitive attribute in the first box
#' @param p2 proportion of marked cards with the sensitive attribute in the second box
#' @param pi vector of the first-order inclusion probabilities
#' @param type the estimator type: total or mean 
#' @param cl confidence level 
#' @param N size of the population. By default it is NULL
#' @param pij matrix of the second-order inclusion probabilities. By default it is NULL
#' 
#' @details 
#' In the Horvitz model, when the population proportion \eqn{\alpha} is not known, two independent samples are taken. Two boxes are filled with a large number of similar 
#' cards except that in the first box a proportion \eqn{p_1(0<p_1<1)} of them is marked \eqn{A} and the complementary proportion \eqn{(1-p_1)} each bearing the mark \eqn{B}, 
#' while in the second box these proportions are \eqn{p_2} and \eqn{1-p_2}, maintaining \eqn{p_2} different from \eqn{p_1}. A sample is chosen and every person sampled is requested 
#' to draw one card randomly from the first box and to repeat this independently with the second box. In the first case, a randomized response should be given, as
#' \deqn{I_i=\left\{\begin{array}{lcc}
#' 1 & \textrm{if card type draws "matches" the sensitive trait } A \textrm{ or the innocuous trait } B \\
#' 0 & \textrm{if there is "no match" with the first box }
#' \end{array}
#' \right.}  
#' and the second case given a randomized response as
#' \deqn{J_i=\left\{\begin{array}{lcc}
#' 1 & \textrm{if there is "match" for the second box} \\
#' 0 & \textrm{if there is "no match" for the second box}
#' \end{array}
#' \right.}  
#' The transformed variable is \eqn{r_i=\frac{(1-p_2)I_i-(1-p_1)J_i}{p_1-p_2}} and the estimated variance is \eqn{\widehat{V}_R(r_i)=r_i(r_i-1)}.
#'
#' @return Point and confidence estimates of the sensitive characteristics using the Horvitz-UB model. The transformed variable is also reported, if required.
#' 
#' @references Chaudhuri, A. (2011). 
#' \emph{Randomized response and indirect questioning techniques in surveys.}
#' Boca Raton: Chapman and Hall, CRC Press. 
#' 
#' @references Greenberg, B.G., Abul-Ela, A.L., Simmons, W.R., Horvitz, D.G. (1969).
#' \emph{The unrelated question RR model: Theoretical framework.}
#' Journal of the American Statistical Association, 64, 520-539.
#' 
#' @references Horvitz, D.G., Shah, B.V., Simmons, W.R. (1967).
#' \emph{The unrelated question RR model.}
#'  Proceedings of the Social Statistics Section of the American Statistical Association. 65-72. Alexandria, VA: ASA.
#'  
#' @seealso \code{\link{HorvitzUBData}}
#' @seealso \code{\link{Horvitz}}
#' @seealso \code{\link{ResamplingVariance}}
#' 
#' @keywords Randomized_response Qualitative Horvitz Greenberg Estimation Variance Transformed_variable Confidence_interval
#' 
#' @examples
#' N=802
#' data(HorvitzUBData)
#' dat=with(HorvitzUBData,data.frame(I,J,Pi))
#' p1=0.6
#' p2=0.7
#' cl=0.95
#' HorvitzUB(dat$I,dat$J,p1,p2,dat$Pi,"mean",cl,N)
#' 
#' @export
HorvitzUB=function(I,J,p1,p2,pi,type=c("total","mean"),cl,N=NULL,pij=NULL){
  call=match.call()
  
  out1=Qualitative(I,12,p1,2,p2,z2=J)
  out2=Estimator(out1,pi,type,cl,N,pij)
  
  out=c(Call=call,out1,out2,Name="Horvitz unknown B",Model="Qualitative",Type=type,Param=list(c("p1"=p1,"p2"=p2)),ConfidenceLevel=cl)
  class(out)="RRT"
  return(out)
}