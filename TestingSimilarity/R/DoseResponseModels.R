################################################################################
#' Implementation of EMAX models
#'
#' Emax model: \deqn{m(d,\theta)=E_0+E_{max}\frac{d}{ED_{50}+d}}
#' 
#' @param d real-valued argument to the function (dose variable)
#' @param e model parameter
#' @return Response value.
################################################################################
emax <- function(d,e){
  (e[1]+(e[2]*d)/(e[3]+d))
}
################################################################################
#' Implementation of linear models
#'
#' Linear model: \deqn{m(d,\theta)=E_0+\delta d}
#' 
#' @param d real-valued argument to the function (dose variable)
#' @param e model parameter
#' @return Response value.
################################################################################
linear <- function(d,e){
  e[1]+d*e[2]
}
################################################################################
#' Implementation of exponential models
#'
#' Exponential model: \deqn{m(d,\theta)=E_0+E_1(exp(d/\delta)-1)}
#' 
#' @param d real-valued argument to the function (dose variable)
#' @param e model parameter
#' @return Response value.
################################################################################
exponential <- function(d,e){
  e[1]+e[2]*(exp(d/e[3])-1)
}
################################################################################
#' Implementation of quadratic models
#'
#' Quadratic model: \deqn{m(d,\theta)=E_0+\beta_1 d+\beta_2 d^2}
#' 
#' @param d real-valued argument to the function (dose variable)
#' @param e model parameter
#' @return Response value.
################################################################################
quadratic <- function(d,e){
  e[1]+d*e[2]+d^2*e[3]
}
################################################################################
#' Implementation of logistic models
#'
#' Logistic model: \deqn{m(d,\theta)=E_0+\frac{E_{max}}{1+exp[(ED_{50}-d)/\delta]}}
#' 
#' @param d real-valued argument to the function (dose variable)
#' @param e model parameter
#' @return Response value.
################################################################################
logistic <- function(d,e){
  e[1]+e[2]/(1+exp((e[3]-d)/e[4]))
}
################################################################################
#' Implementation of Sigmoid Emax models
#'
#' Sigmoid Emax Model model: \deqn{m(d,\theta)=E_0+E_{max} \frac{d^h}{ED_{50}^h+d^h}}
#' 
#' @param d real-valued argument to the function (dose variable)
#' @param e model parameter
#' @return Response value
################################################################################
sigEmax <- function(d,e){
  e[1]+e[2]*d^e[4]/(e[3]^e[4]+d^e[4])
}
################################################################################
#' Implementation of Beta models
#'
#' Beta model: \deqn{m(d,\theta)=E_0+E_{max}B(\delta_1,\delta_2)(d/scal)^{\delta_1}(1-d/scal)^{\delta_2}}
#' with \deqn{B(\delta_1,\delta_2)=(\delta_1+\delta_2)^{\delta_1+\delta_2}/(\delta_1^{\delta_1} \delta_2^{\delta_2})}
#' and \eqn{scal} is a fixed dose scaling parameter.
#' 
#' @param d real-valued argument to the function (dose variable)
#' @param e model parameter
#' @param scal fixed dose scaling parameter
#' @return Response value.
################################################################################
betaMod <- function(d,e,scal){
  e[1]+e[2]*(e[3]+e[4])^(e[3]+e[4])/(e[3]^e[3]*e[4]^e[4])*(d/scal)^e[3]*(1-d/scal)^e[4]
}
################################################################################
#' Implementation of linear in log models
#'
#' Linear in log Model model: \deqn{m(d,\theta)=E_0+\delta\ log(d+off)}
#' and \eqn{off} is a fixed offset parameter.
#' 
#' @param d real-valued argument to the function (dose variable)
#' @param e model parameter
#' @param off fixed offset parameter
#' @return Response value.
################################################################################
linlog <- function(d,e,off){
  
  e[1]+e[2]*log(d+off)
}