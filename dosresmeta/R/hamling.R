#' Approximating effective-counts as proposed by Hamling
#' 
#' @description The function \code{hamling} reconstructs the set of pseudo-numbers (or "effective" numbers) of cases and non-cases consistent
#' with the input data (log relative risks) for either a single or multiple summarized data.
#' The method was proposed in 2008 by Hamling.
#' 
#' @param logrr a vector to specify the reported log relative risks.
#' @param v a vector to specify the variances of the reported log relative risks. Alternatively the user can provide the standard error in the \code{se} argument, or the confidence interval for the reported relative risks 
#' in the \code{lb} and \code{ub} arguments.
#' @param cases a vector to specify the number of cases for each exposure level.
#' @param n a vector to specify the total number of subject for each exposure level. For incidence-rate data \code{n} indicates the amount of person-time for each exposure level.
#' @param type a vector (or a string) to specify the study-specific design. The values for case-control, incidence-rate, and cumulative incidence data are \code{cc},
#' \code{ir}, and \code{ci} (or 1, 2, and 3), respectively.
#' @param id an optional vector to specify the id variable for the studies included in the analysis.
#' @param data an optional data frame (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the previous arguments.
#' @param se an optional vector to specify the standard error of the reported log relative risks; needed if \code{v} is not provided.
#' @param lb an optional vector to specify the lower bound of the confidence interval for the reported relative risks; needed if \code{v} and \code{se} are not provided.
#' @param ub an optional vector to specify the upper bound of the confidence interval for the reported relative risks; needed if \code{v} and \code{se} are not provided.
#' @param order a logical value to specify if the vectors need to be sorted. See details.
#' @param alpha a scalar to specify the alpha nominal value used in the published data, by defaul equal to .05.
#' 
#' @details The function reconstructs the effective counts corresponding to the multivariable adjusted log relative risks as well as their standard errors. A unique solution is guaranteed by keeping the ratio 
#' non-cases to cases and the fraction of unexposed subjects equal to the unadjusted data (Hamling). The function requires the data to be sorted by id and in such a way that the referent values correspond to the first record for each study. 
#' This is automatically done by the function when \code{order = TRUE} (default).
#'
#' @return The results are returned structured in a data frame.
#'
#' @examples
#' data("alcohol_cvd")
#' hamling(logrr = logrr, se = se, cases = cases, n = n, type = type,
#'  id = id, data = alcohol_cvd)
#'  
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' 
#' 
#' @references 
#' Hamling, J., Lee, P., Weitkunat, R., Ambuhl, M. (2008). Facilitating meta-analyses by deriving relative effect and precision estimates for alternative comparisons from a set of estimates presented
#' by exposure level or disease category. Statistics in medicine, 27(7), 954-970.
#' 
#' Orsini, N., Li, R., Wolk, A., Khudyakov, P., Spiegelman, D. (2012). Meta-analysis for linear and nonlinear dose-response relations: examples, an evaluation of approximations, and software. 
#' American journal of epidemiology, 175(1), 66-73.
#' 
#' @export hamling



hamling <- function (logrr, v, cases, n, type, id, data, 
                     se, lb, ub, order = TRUE, alpha = .05)
{
  if (missing(data)) data <- NULL
  if (is.null(data)) {
    data <- sys.frame(sys.parent())
  } else {
    if (!is.data.frame(data)) {
      data <- data.frame(data)
    }
  }
  mf <-  match.call(expand.dots = FALSE)
  mf.logrr <- mf[[match("logrr", names(mf))]]
  mf.v <- mf[[match("v", names(mf))]]
  mf.se <- mf[[match("se", names(mf))]]
  mf.cases <- mf[[match("cases", names(mf))]]
  mf.n <- mf[[match("n", names(mf))]]
  mf.type <- mf[[match("type", names(mf))]]
  mf.id <- mf[[match("id", names(mf))]]
  mf.lb <- mf[[match("lb", names(mf))]]
  mf.ub <- mf[[match("ub", names(mf))]]
  logrr <- eval(mf.logrr, data, enclos = sys.frame(sys.parent()))
  if (is.null(mf.v)){
    se <- eval(mf.se, data, enclos = sys.frame(sys.parent()))
    if (is.null(mf.se)){
      lb <- eval(mf.lb, data, enclos = sys.frame(sys.parent()))
      ub <- eval(mf.ub, data, enclos = sys.frame(sys.parent()))
      v <- (log(ub) - log(lb)) / (2 * qnorm(1 - alpha/2))^2
    } else { v <- se^2
    }
  } else {
    v <- eval(mf.v, data, enclos = sys.frame(sys.parent()))
  }
  v[is.na(v)] <- 0
  cases <- eval(mf.cases, data, enclos = sys.frame(sys.parent()))
  n <- eval(mf.n, data, enclos = sys.frame(sys.parent()))
  type <- eval(mf.type, data, enclos = sys.frame(sys.parent()))
  if (is.null(type)) type <- as.vector(mf.type)
  id <- as.factor(eval(mf.id, data, enclos = sys.frame(sys.parent())))
  if (is.null(mf.id)) id <- as.factor(rep(1, length(v)))
  index <- order(id, v)
  in1 <- row.names(data.frame(id, v))
  index <- order(id, v)
  if(order == T){
    logrr <- logrr[index]
    cases <- cases[index]
    n <- n[index]
    type <- type[index]
    v <- v[index]
    id <- id[order(id)]
  }
  fun <- function(cc, ris=F){
    A <- cc[1]
    B <- cc[2]
    if (type[id==j][1] == "cc" | type[id==j][1] == 1) {
      cases1 <- (1 + (A/B) * exp(logrr[id==j & v!=0])) / (v[id==j & v!=0] - 1/A -1/B)
      n1 <- (1 + B / (A * exp(logrr[id==j & v!=0]))) / (v[id==j & v!=0] - 1/A - 1/B)
    }
    if ((type[id==j][1] == "ir" | type[id==j][1] == 2) |
          (type[id==j][1] == "ci" | type[id==j][1] == 3)){ 
      cases1 <- 1 / (v[id==j & v!=0] - 1/A)
      n1 <- (B / ( A * exp(logrr[id==j & v!=0]) ) ) / (v[id==j & v!=0] - 1/A )
    }
    if (ris==F) {
      p1 <- B / (B + sum(n1))
      z1<- (B + sum(n1)) / (A + sum(cases1))
      (((p - p1) / p)^2 + ((z - z1) / z)^2)
    } else {
      return(data.frame(cases1 = c(A, cases1), n1 = c(B, n1)))
    }
  }
  pscounts <- data.frame()
  for (j in unique(id)){
    if (type[id==j][1] == "cc" | type[id==j][1] == 1){
      controls <- n[id==j] - cases[id==j]
      p <- controls[v[id==j]==0]/sum(controls)
      z <- sum(controls)/sum(cases[id==j])
      st.in <- c(cases[id==j & v==0], controls[v[id==j]==0])
    }
    if ( (type[id==j][1] == "ir" | type[id==j][1] == "ci") |
           (type[id==j][1] == 2 | type[id==j][1] == 3) ){
      p <- n[id==j & v==0]/sum(n[id==j])
      z <- sum(n[id==j])/sum(cases[id==j])
      st.in <- c(cases[id==j & v==0], n[id==j & v==0])
    }
    st.fin <- optim(st.in, fun)$par
    pscounts <- rbind(pscounts, fun(st.fin, ris=T))
  }
  rownames(pscounts) <- index
  if (order == T) pscounts <- pscounts[in1, ]
  data.frame(id, pscounts)
}
