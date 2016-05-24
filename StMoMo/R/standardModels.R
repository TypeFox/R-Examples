#' Create a Lee-Carter model
#' 
#' Utility function to initialise a \code{StMoMo} object representing a 
#' Lee-Carter model.
#' 
#' The created model is either a log-Poisson (see Brouhns et al (2002)) or a 
#' logit-Binomial version of the Lee-Carter model which has predictor structure   
#' \deqn{\eta_{xt} = \alpha_x + \beta_x\kappa_t.}
#' To ensure identifiability one of the  following constraints is imposed
#' \deqn{\sum_t\kappa_t = 0,\,\kappa_1 = 0,\, \kappa_n = 0}
#' depending on the value of \code{const}, and
#' \deqn{\sum_x\beta_x = 1.}
#' 
#' @inheritParams StMoMo
#' @param const defines the constraint to impose to the period index of the
#'  model to ensure identifiability. The alternatives are 
#'  \code{"sum"}(default),  \code{"last"} and \code{"first"} which apply 
#'  constraints \eqn{\sum_t\kappa_t = 0}, \eqn{\kappa_n = 0} and 
#'  \eqn{\kappa_1 = 0}, respectively.
#' 
#' @return An object of class \code{"StMoMo"}.
#' 
#' @seealso \code{\link{StMoMo}}
#'  
#' @references
#' Brouhns, N., Denuit, M., & Vermunt, J. K. (2002). A Poisson log-bilinear 
#' regression approach to the construction of projected lifetables.
#'  Insurance: Mathematics and Economics, 31(3), 373-393.
#' 
#' Lee, R. D., & Carter, L. R. (1992). Modeling and forecasting U.S. mortality. 
#' Journal of the American Statistical Association, 87(419), 659-671. 
#' 
#' @examples
#' 
#' #sum(kt) = 0 and log link
#' LC1 <- lc()
#' LCfit1<-fit(LC1, Dxt = EWMaleData$Dxt,Ext = EWMaleData$Ext, 
#'             ages = EWMaleData$ages, years = EWMaleData$years,
#'             ages.fit = 55:89)
#' plot(LCfit1)
#' 
#' #kt[1] = 0 and log link
#' LC2 <- lc(const = "first")
#' LCfit2<-fit(LC2, Dxt = EWMaleData$Dxt,Ext = EWMaleData$Ext, 
#'             ages = EWMaleData$ages, years = EWMaleData$years,
#'             ages.fit = 55:89)
#' plot(LCfit2)
#' 
#' #kt[n] = 0 and logit link
#' LC3 <- lc("logit", "last")
#' LCfit3<-fit(LC3, Dxt = EWMaleData$Dxt,Ext = EWMaleData$Ext, 
#'             ages = EWMaleData$ages, years = EWMaleData$years,
#'             ages.fit = 55:89)
#' plot(LCfit3)
#' 
#' @export
lc <- function(link = c("log", "logit"), const = c("sum", "last", "first")) {
  link <- match.arg(link)
  const <- match.arg(const)  
  constLC <- function(ax, bx, kt, b0x, gc, wxt, ages) {    
    c1 <- switch(const, sum = mean(kt[1, ], na.rm = TRUE),
                 first = kt[1, 1], last = tail(kt[1, ], 1))
    ax <- ax + c1 * bx[, 1]
    kt[1, ] <- kt[1, ] - c1
    c2 <- sum(bx[, 1], na.rm = TRUE)
    bx[, 1] <- bx[, 1] / c2
    kt[1, ] <- kt[1, ] * c2
    list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
  }
  StMoMo(link = link, staticAgeFun = TRUE, periodAgeFun = "NP", 
         constFun = constLC)    
}

#' Create a Cairns-Blake-Dowd mortality model
#' 
#' Utility function to initialise a \code{StMoMo} object representing a 
#' Cairns-Blake-Dowd mortality model.
#' 
#' The created model is either a logit-Binomial or a log-Poisson version of 
#' the Cairns-Blake-Dowd mortality model which has predictor structure 
#' \deqn{\eta_{xt} = \kappa_t^{(1)} + (x-\bar{x})\kappa_t^{(2)},}
#' where \eqn{\bar{x}} is the average age in the data.
#' 
#' @param link defines the link function and random component associated with 
#'   the mortality model. \code{"log"} would assume that deaths follow a 
#'   Poisson distribution and use a log link while \code{"logit"} would 
#'   assume that deaths follow a Binomial distribution and a logit link. 
#'   Note that the default is the logit link.
#' @return An object of class \code{"StMoMo"}.
#' 
#' @seealso \code{\link{StMoMo}}, \code{\link{m6}}, \code{\link{m7}}, 
#' \code{\link{m8}}
#'  
#' @references
#' Cairns, A. J. G., Blake, D., & Dowd, K. (2006). A Two-Factor Model for 
#' Stochastic Mortality with Parameter Uncertainty: Theory and Calibration. 
#' Journal of Risk and Insurance, 73(4), 687-718.
#' 
#' @examples
#' 
#' CBD <- cbd()
#' Dxt <- EWMaleData$Dxt
#' Ext <- EWMaleData$Ext + 0.5 * EWMaleData$Dxt
#' CBDfit <- fit(CBD, Dxt = Dxt, Ext = Ext, 
#'             ages = EWMaleData$ages, years = EWMaleData$years,
#'             ages.fit = 55:89)
#' plot(CBDfit, parametricbx = FALSE)
#' 
#' @export
cbd <- function(link = c("logit", "log")) {
  link <- match.arg(link)
  f1 <- function(x,ages) x - mean(ages)
  StMoMo(link = link, staticAgeFun = FALSE, periodAgeFun=c("1", f1))
}


#' Create an Age-Period-Cohort mortality model
#' 
#' Utility function to initialise a \code{StMoMo} object representing an 
#' Age-Period-Cohort mortality model.
#' 
#' The created model is either a log-Poisson or a logit-Binomial version of 
#' the classical age-period-cohort mortality model which has predictor 
#' structure 
#' \deqn{\eta_{xt} = \alpha_x + \kappa_t + \gamma_{t-x}.}
#' 
#' To ensure identifiability we follow Cairns et al. (2009) and impose 
#' constraints \deqn{\sum_c \gamma_c = 0}  and  \deqn{\sum_c c\gamma_c = 0}
#' 
#' @inheritParams StMoMo
#' @return An object of class \code{"StMoMo"}.
#' 
#' @seealso \code{\link{StMoMo}}, \code{\link{rh}}
#' 
#' @references
#' 
#' Cairns, A. J. G., Blake, D., Dowd, K., Coughlan, G. D., Epstein, D., 
#' Ong, A., & Balevich, I. (2009). A quantitative comparison of stochastic 
#' mortality models using data from England and Wales and the United States. 
#' North American Actuarial Journal, 13(1), 1-35.
#' 
#' @examples
#' 
#' APC <- apc()
#' wxt <- genWeightMat(EWMaleData$ages,  EWMaleData$years, clip = 3)
#' APCfit <- fit(APC, Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'               ages = EWMaleData$ages, years = EWMaleData$years, 
#'               wxt = wxt)
#' plot(APCfit, parametricbx = FALSE, nCol = 3)
#' 
#' @export
apc <- function(link = c("log", "logit")) {
  link <- match.arg(link)
  constAPC <- function(ax, bx, kt, b0x, gc, wxt, ages) {    
    nYears <- dim(wxt)[2]  
    x <- ages  
    t <- 1:nYears
    c <- (1 - tail(ages, 1)):(nYears - ages[1])    
    #\sum g(c)=0  and  \sum cg(c)=0       
    phiReg <- lm(gc ~ 1 + c, na.action = na.omit)
    phi <- coef(phiReg)       
    gc <- gc - phi[1] - phi[2] * c
    ax <- ax + phi[1] - phi[2] *x 
    kt <- kt + phi[2] * t   
    #\sum k(t)=0 
    c1 <- mean(kt, na.rm = TRUE)
    kt <- kt - c1
    ax <- ax + c1    
    list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc) 
  }
  StMoMo(link = link, staticAgeFun = TRUE, periodAgeFun = "1",
                cohortAgeFun = "1", constFun = constAPC)

}



#' Create a Renshaw and Haberman (Lee-Carter with cohorts) mortality model
#' 
#' Utility function to initialise a \code{StMoMo} object representing a 
#' Renshaw and Haberman (Lee-Carter with cohorts) mortality model introduced
#' in Renshaw and Haberman (2006).
#' 
#' The created model is either a log-Poisson or a 
#' logit-Binomial version of the Renshaw and Haberman model which has 
#' predictor structure   
#' \deqn{\eta_{xt} = \alpha_x + \beta^{(1)}_x\kappa_t + \beta^{(0)} \gamma_{t-x}.}
#' or
#' \deqn{\eta_{xt} = \alpha_x + \beta^{(1)}_x\kappa_t + \gamma_{t-x}.}
#' depending on the value of argument \code{cohortAgeFun}.
#'   
#' To ensure identifiability the  following constraints are imposed
#' \deqn{\sum_t\kappa_t = 0, \sum_x\beta^{(1)}_x = 1, \sum_c\gamma_c = 0}
#' plus
#' \deqn{\sum_x\beta^{(0)}_x = 1}
#' if \code{cohortAgeFun = "NP"}
#' 
#' By default \eqn{\beta^{(0)}_x = 1} as this model has shown to be more
#' stable (see Haberman and Renshaw (2011) and Hunt and Villegas (2015)).
#' 
#' @inheritParams StMoMo
#' @param cohortAgeFun defines the cohort age modulating parameter 
#'   \eqn{\beta_x^{(0)}}. It can take values: \code{"NP"} for a non-parametric age
#'   term or \code{"1"} for \eqn{\beta_x^{(0)}=1} (the default). 
#' 
#' @return An object of class \code{"StMoMo"}.
#' 
#' @seealso \code{\link{StMoMo}}, \link{lc}, \link{apc}
#'  
#' @references
#'
#' Haberman, S., & Renshaw, A. (2011). A comparative study of parametric 
#' mortality projection models. Insurance: Mathematics and Economics, 
#' 48(1), 35-55. 
#' 
#' Hunt, A., & Villegas, A. M. (2015). Robustness and convergence in the 
#' Lee-Carter model with cohorts. Insurance: Mathematics and Economics, 
#' 64, 186-202. 
#' 
#' Renshaw, A. E., & Haberman, S. (2006). A cohort-based extension to the 
#' Lee-Carter model for mortality reduction factors. 
#' Insurance: Mathematics and Economics, 38(3), 556-570.
#' 
#' @examples
#' 
#' LCfit <-  fit(lc(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'                ages = EWMaleData$ages, years = EWMaleData$years, 
#'                ages.fit = 55:89)
#' wxt <- genWeightMat(55:89,  EWMaleData$years, clip = 3)
#' RHfit <- fit(rh(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'               ages = EWMaleData$ages, years = EWMaleData$years, 
#'               ages.fit = 55:89, wxt = wxt, start.ax = LCfit$ax, 
#'               start.bx = LCfit$bx, start.kt = LCfit$kt)
#' plot(RHfit)
#' 
#' @export
rh <- function(link = c("log", "logit"), cohortAgeFun = c("1", "NP")) {
  link <- match.arg(link)
  cohortAgeFun <- match.arg(cohortAgeFun)
  constRHgeneral <- function(ax, bx, kt, b0x, gc, wxt, ages, cohortAgeFun) {
    #\sum k[t] = 0
    c1 <- mean(kt[1, ], na.rm = TRUE) 
    ax <- ax + c1 * bx[, 1]
    kt[1, ] <- kt[1, ] - c1
    #\sum b[x, 1] = 0
    c2 <- sum(bx[, 1], na.rm = TRUE)
    bx[, 1] <- bx[, 1] / c2
    kt[1, ] <- kt[1, ] * c2
    #\sum g[c] = 0
    c3 <- mean(gc, na.rm = TRUE)
    ax <- ax + c3 * b0x
    gc <- gc - c3
    #\sum b0[x] = 1
    if (cohortAgeFun == "NP") {
      c4 <- sum(b0x, na.rm = TRUE)
      b0x <- b0x / c4
      gc <- gc * c4
    }
    list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
  }
  constRH <- function(ax, bx, kt, b0x, gc, wxt, ages) {
    constRHgeneral(ax, bx, kt, b0x, gc, wxt, ages, cohortAgeFun) 
  }
  StMoMo(link = link, staticAgeFun = TRUE, periodAgeFun = "NP",
         cohortAgeFun = cohortAgeFun, constFun = constRH)  
}


#' Create an M6 type extension of the Cairns-Blake-Dowd mortality model
#' 
#' Utility function to initialise a \code{StMoMo} object representing the 
#' M6 (CBD with cohorts) extension of the Cairns-Blake-Dowd mortality model 
#' introduced in Cairns et al (2009).
#' 
#' The created model is either a logit-Binomial or a log-Poisson version of the 
#' M6 model which has predictor structure 
#' \deqn{\eta_{xt} = \kappa_t^{(1)} + (x-\bar{x})\kappa_t^{(2)} + \gamma_{t-x},} 
#' where \eqn{\bar{x}} is the average age in the data.
#' 
#' Identifiability of the model is accomplished by applying parameters 
#' constraints
#' \deqn{\sum_c\gamma_c = 0, \sum_c c\gamma_c = 0}
#' which ensure that the cohort effect fluctuates around zero and has no 
#' linear trend. These constraints are applied using the strategy discussed 
#' in Appendix A of Haberman and Renshaw (2011).
#' 
#' @inheritParams cbd
#' @return An object of class \code{"StMoMo"}.
#' 
#' @seealso \code{\link{StMoMo}}, \code{\link{cbd}}, \code{\link{m7}}, 
#' \code{\link{m8}}
#' 
#' @references
#' 
#' Cairns, A. J. G., Blake, D., Dowd, K., Coughlan, G. D., Epstein, D., 
#' Ong, A., & Balevich, I. (2009). A quantitative comparison of stochastic 
#' mortality models using data from England and Wales and the United States. 
#' North American Actuarial Journal, 13(1), 1-35.
#' 
#' Haberman, S., & Renshaw, A. (2011). A comparative study of parametric 
#' mortality projection models. Insurance: Mathematics and Economics, 
#' 48(1), 35-55. 
#' 
#' @examples
#' 
#' M6 <- m6()
#' Dxt <- EWMaleData$Dxt
#' Ext <- EWMaleData$Ext + 0.5 * EWMaleData$Dxt
#' wxt <- genWeightMat(55:89,  EWMaleData$years, clip = 3)
#' M6fit <- fit(M6, Dxt = Dxt, Ext = Ext, ages = EWMaleData$ages, 
#'            years = EWMaleData$years, ages.fit = 55:89)
#' plot(M6fit, parametricbx = FALSE)
#' 
#' @export
m6 <- function(link = c("logit", "log")) {
  link <- match.arg(link)
  f1 <- function(x,ages) x - mean(ages)
  constM6 <- function(ax, bx, kt, b0x, gc, wxt, ages) {
    #See Appendix A in Haberman and Renshaw (2011)
    nYears <- dim(wxt)[2]
    x <- ages
    t <- 1:nYears
    c <- (1 - tail(ages, 1)):(nYears - ages[1])
    xbar <- mean(x)
    #\sum g(c)=0  and  \sum cg(c)=0
    phiReg <- lm(gc ~ 1 + c, na.action = na.omit)
    phi <- coef(phiReg)
    gc <- gc - phi[1] - phi[2] * c    
    kt[2, ] <- kt[2, ] - phi[2]
    kt[1, ] <- kt[1, ] + phi[1] + phi[2] * (t - xbar)
    list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
  }
  StMoMo(link = link, staticAgeFun = FALSE, periodAgeFun = c("1", f1),
               cohortAgeFun = "1", constFun = constM6)
  
}



#' Create an M7 type extension of the Cairns-Blake-Dowd mortality model
#' 
#' Utility function to initialise a \code{StMoMo} object representing the 
#' M7 extension of the Cairns-Blake-Dowd mortality model introduced
#' in Cairns et al (2009).
#' 
#' The created model is either a logit-Binomial or a log-Poisson version of 
#' the M7 model which has predictor structure 
#' \deqn{\eta_{xt} = \kappa_t^{(1)} + (x-\bar{x})\kappa_t^{(2)} + 
#'                            ((x-\bar{x})^2 - \hat{\sigma}^2_x)\kappa_t^{(2)} +  \gamma_{t-x},} 
#' where \eqn{\bar{x}} is the average age in the data and \eqn{\hat{\sigma}^2_x} 
#' is the average value of \eqn{(x-\bar{x})^2}.
#'                            
#' Identifiability of the model is accomplished by applying parameters 
#' constraints \deqn{\sum_c\gamma_c = 0, \sum_c c\gamma_c = 0, 
#' \sum_c c^2\gamma_c = 0} which ensure that the cohort effect fluctuates 
#' around zero and has no linear or quadratic trend. These constraints are 
#' applied using the strategy discussed in Appendix A of 
#' Haberman and Renshaw (2011).
#' 
#' @inheritParams cbd
#' @return An object of class \code{"StMoMo"}.
#' 
#' @seealso \code{\link{StMoMo}}, \code{\link{cbd}}, \code{\link{m6}}, 
#' \code{\link{m8}}
#' 
#' @references
#' 
#' Cairns, A. J. G., Blake, D., Dowd, K., Coughlan, G. D., Epstein, D.,
#' Ong, A., & Balevich, I. (2009). A quantitative comparison of stochastic 
#' mortality models using data from England and Wales and the United States. 
#' North American Actuarial Journal, 13(1), 1-35.
#' 
#' Haberman, S., & Renshaw, A. (2011). A comparative study of parametric 
#' mortality projection models. Insurance: Mathematics and Economics, 
#' 48(1), 35-55. 
#' 
#' @examples
#' 
#' M7 <- m7()
#' Dxt <- EWMaleData$Dxt
#' Ext <- EWMaleData$Ext + 0.5 * EWMaleData$Dxt
#' wxt <- genWeightMat(55:89,  EWMaleData$years, clip = 3)
#' M7fit <- fit(M7, Dxt = Dxt, Ext = Ext, ages = EWMaleData$ages, 
#'              years = EWMaleData$years, ages.fit = 55:89)
#' plot(M7fit, parametricbx = FALSE)
#' 
#' @export
m7 <- function(link = c("logit", "log")) {
  link <- match.arg(link)
  f1 <- function(x,ages) x - mean(ages)
  f2 <- function(x,ages) {
    xbar <- mean(ages)
    s2 <- mean((ages - xbar)^2)
    (x - xbar)^2 - s2
  }
  constM7 <- function(ax, bx, kt, b0x, gc, wxt, ages) {
    #See Appendix A in Haberman and Renshaw (2011)
    nYears <- dim(wxt)[2]
    x <- ages
    t <- 1:nYears
    c <- (1 - tail(ages, 1)):(nYears - ages[1])
    xbar <- mean(x)
    s2 <- mean((x - xbar)^2)
    #\sum g(c)=0, \sum cg(c)=0, \sum c^2g(c)=0
    phiReg <- lm(gc ~ 1 + c + I(c^2), na.action = na.omit)
    phi <- coef(phiReg)
    gc <- gc - phi[1] - phi[2] * c - phi[3] * c^2    
    kt[3, ] <- kt[3, ] + phi[3]
    kt[2, ] <- kt[2, ] - phi[2] - 2 * phi[3] * (t - xbar)  
    kt[1, ] <- kt[1, ] + phi[1] + phi[2] * (t - xbar) + 
      phi[3] * ((t - xbar)^2 + s2)
    list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
  }
  StMoMo(link = link, staticAgeFun = FALSE, periodAgeFun = c("1", f1, f2),
         cohortAgeFun = "1", constFun = constM7)
  
}


#' Create an M8 type extension of the Cairns-Blake-Dowd mortality model
#' 
#' Utility function to initialise a \code{StMoMo} object representing the 
#' M8 extension of the Cairns-Blake-Dowd mortality model introduced
#' in Cairns et al (2009).
#' 
#' The created model is either a logit-Binomial or a log-Poisson version of 
#' the M8 model which has predictor structure 
#' \deqn{\eta_{xt} = \kappa_t^{(1)} + (x-\bar{x})\kappa_t^{(2)} + (x_c-x)\gamma_{t-x}}
#' where \eqn{\bar{x}} is the average age in the data and \eqn{x_c} is a
#' predefined constant. 
#' Identifiability of the model is accomplished by applying parameters 
#' constraint
#' \deqn{\sum_c\gamma_c = 0.}
#' 
#' @inheritParams cbd
#' @param xc constant defining the cohort age-modulating parameter. 
#' 
#' @return An object of class \code{"StMoMo"}.
#' 
#' @seealso \code{\link{StMoMo}}, \code{\link{cbd}}, \code{\link{m6}}, 
#' \code{\link{m7}}
#' 
#' @references
#' 
#' Cairns, A. J. G., Blake, D., Dowd, K., Coughlan, G. D., Epstein, D., 
#' Ong, A., & Balevich, I. (2009). A quantitative comparison of stochastic
#' mortality models using data from England and Wales and the United States. 
#' North American Actuarial Journal, 13(1), 1-35.
#' 
#' @examples
#' 
#' M8 <- m8(xc = 89)
#' Dxt <- EWMaleData$Dxt
#' Ext <- EWMaleData$Ext + 0.5 * EWMaleData$Dxt
#' wxt <- genWeightMat(55:89,  EWMaleData$years, clip = 3)
#' M8fit <- fit(M8, Dxt = Dxt, Ext = Ext, ages = EWMaleData$ages, 
#'              years = EWMaleData$years, ages.fit = 55:89)
#' plot(M8fit, parametricbx = FALSE)
#' 
#' @export
m8 <- function(link = c("logit", "log"), xc) {
  link <- match.arg(link)
  f1 <- function(x,ages) x - mean(ages)
  f3 <- function(x,ages) xc - x
  constM8 <- function(ax, bx, kt, b0x, gc, wxt, ages) {
    #See Appenix A in Haberman and Renshaw (2011)
    x <- ages
    xbar <- mean(x)      
    c <- mean(gc,na.rm = TRUE)  
    gc <- gc - c
    kt[1, ] <- kt[1, ] + c * (xc - xbar)
    kt[2, ] <- kt[2, ] - c 
    list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
  }
  StMoMo(link = link, staticAgeFun = FALSE, periodAgeFun = c("1", f1),
         cohortAgeFun = f3, constFun = constM8)
  
}
