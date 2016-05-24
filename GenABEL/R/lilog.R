#' LiLog (Linear to Logistic Beta estimator)
#'
#' Estimating Logistic Betas from a linear regression.
#'
#' The function transforms the betas from the linear regression of the
#' binomial trait into its logistic counterpart.
#' It is based on the following relationship: if we denote with
#' \eqn{\alpha} the coefficients of the logistic regression and with
#' \eqn{\beta} those coming from the linear regression we have that:
#'
#' \deqn{z=\alpha_0 + \alpha_{SNP} + \alpha_{COV1} \ldots +
#'  \alpha_{COVi} = \ln \left( \frac{p}{1-p} \right),}{%
#'  z = alpha0 + alphaSNP + alphaCOV1 + ... + alphaCOVi = ln(p/(1-p)),}
#' with
#' \deqn{p = \beta_0 + \beta_{SNP}}{%
#'       p = beta0 + betaSNP}
#' Covariates are disregarded and the relationship
#' \deqn{\alpha_0 + \alpha_{SNP} = \mathrm{logit}(\beta_0 + \beta_{SNP})}{%
#' alpha0 + \alphaSNP = logit(beta0 + betaSNP)}
#' holds.
#'
#' \eqn{\beta_0} is not the coefficient coming from the linear
#' regression but it is estimated as:
#' \deqn{\mathrm{prev} - 2  q  (1-q)  \beta - q^2  \beta}{%
#'       prev - 2 * q * (1-q) * beta - q^2 * beta}
#'
#' and \eqn{\alpha_0} is estimated as
#' \eqn{\mathrm{logit}(\beta_0)}{logit(beta0)}.
#' So \eqn{\alpha_{SNP}}, which is the parameter of interest, will be:
#' \deqn{\alpha_{SNP} = \mathrm{logit}(\beta_0 + \beta_{SNP}) -
#'            \mathrm{logit}(\beta_0).}{%
#'       alphaSNP = logit(beta0 + betaSNP) - logit(beta0).}
#'
#' \eqn{SE(\alpha_{SNP})}{SE(\alphaSNP)} can be estimated considering
#' that the following relationship must be true:
#' \deqn{\frac{\alpha_{SNP}}{SE(\alpha_{SNP})} =
#'       \frac{\beta_{SNP}}{SE(\beta_{SNP})},}{%
#'       alphaSNP/SEalphaSNP = betaSNP/SEbetaSNP,} which is equivalent to
#' \deqn{SE(\alpha_{SNP}) = \alpha_{SNP}
#'         \frac{SE(\beta_{SNP})}{\beta_{SNP}}.}{%
#'       SE(alphaSNP) = alphaSNP * SE(betaSNP)/betaSNP.}
#'
#' @param prev: Prevalence of the disease as observed from the data
#' @param beta: \eqn{\beta_{SNP}}{beta_SNP} as estimated from a linear
#'  regression of the binomial trait
#' @param SEbeta: Standard Error of \eqn{\beta}{beta} coming from the
#'  linear regression
#' @param q: Frequency of the coded allele.
#'
#' @return Dataframe containing \eqn{\beta}{beta} and
#' \eqn{SE(\beta)}{SE(beta)}.
#'
#' @author Nicola Pirastu, Lennart Karssen, Yurii Aulchenko
#'
#' @references
#'
#' \emph{Genome-wide feasible and meta-analysis oriented methods for
#' analysis of binary traits using mixed models},
#' Nicola Pirastu, Lennart C. Karssen, Pio D'adamo, Paolo Gasparini,
#' Yurii Aulchenko (Submitted).
#'
#' @seealso
#' For transforming results from \code{\link{formetascore}} see
#' \code{\link{formetascore2lilog}}
#' while for tranforming the output from ProbABEL use \code{\link{palinear2LiLog}}
#'
#' @examples
#'
#' # simulate a snp with q = 0.3
#' snp <- rbinom(n=1000, size=2, prob=0.3)
#' # estimate the probability according to a logistic model with beta 0.4
#' probabil <- 1/(1 + exp( -(snp * 0.4 - 1) ))
#' # simulate the binary trait
#' bin.trait <- rbinom(n=1000, size=1, prob=probabil)
#' # run linear regression
#' res.lin <- glm(bin.trait~snp)
#' # run logistic regression
#' res.log <- glm(bin.trait~snp, family=binomial)
#'
#' # recover information from linear regression
#' snp.info <- summary(res.lin)$coefficients["snp",]
#'
#' # transform linear regression coefficients into logistic regression scale
#' LiLog(prev=mean(bin.trait), beta=snp.info[1], SEbeta=snp.info[2],
#'       q=mean(snp)/2)
#'
#' summary(res.log)$coefficients["snp", 1:2]
#'
#' # Example with qtscore
#' data(srdta)
#'
#' snp <- as.numeric(srdta[, 100])
#' probabil <- 1/(1 + exp( -(snp * 0.4 - 1) ))
#' bin.trait <- rbinom(n=nids(srdta), size=1, prob=probabil)
#' srdta@phdata$binary <- bin.trait
#' res <- qtscore(binary, data=srdta)
#' sum <- summary(srdta)
#'
#' lilog.res <- LiLog(prev=mean(bin.trait), beta=res@results$effB,
#'                    SEbeta=res@results$se_effB, q=sum$Q.2)
#' summary(glm(bin.trait~snp,
#'         family=binomial)$coefficients)["snp"1:2,] # logistic regression
#' lilog.res[100, ] # lilog result
#'
#' @keywords
#' LiLog
#'

LiLog <- function(prev, beta, SEbeta, q){
  Icept <- (prev - 2 * q * beta)
  alfa0 <- log(Icept / (1-Icept))
  T1    <- (Icept + beta)
  T1[which(T1<0 | T1>1)] <- NA

  log.beta   <- binomial()$linkfun(T1) - alfa0
  se.log     <- (log.beta * SEbeta) / beta
  res        <- data.frame(log.beta, se.log)
  names(res) <- c("beta", "sebeta")
  return(res)
}
