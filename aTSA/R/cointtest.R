#' Cointegration Test
#' @description Performs Engle-Granger(or EG) tests for the null hypothesis that two or more 
#' time series, each of which is I(1), are not cointegrated.
#' @param y the response
#' @param X the exogenous input variable of a numeric vector or a matrix.
#' @param d difference operator for both \code{y} and \code{X}. The default is 0.
#' @param nlag the lag order to calculate the test statistics. The default is \code{NULL}.
#' @param output a logical value indicating to print the results in R console. 
#' The default is \code{TRUE}.
#' 
#' @details To implement the original EG tests, one first has to fit the linear regression
#' \deqn{y[t] = \mu + B*X[t] + e[t],}
#' where \eqn{B} is the coefficient vector and \eqn{e[t]} is an error term.
#' With the fitted model, the residuals are obtained, i.e., \eqn{z[t] = y[t] - hat{y}[t]} 
#' and a Augmented Dickey-Fuller test is utilized to examine whether the sequence of 
#' residuals \eqn{z[t]} is white noise. The null hypothesis of non-cointegration 
#' is equivalent to the null hypothesis that \eqn{z[t]} is white noise. See \code{\link{adf.test}} 
#' for more details of Augmented Dickey-Fuller test, as well as the default \code{nlag}. 
#' 
#' @return A matrix for test results with three columns (\code{lag}, \code{EG}, \code{p.value})
#'  and three rows (\code{type1}, \code{type2}, \code{type3}). 
#'  Each row is the test results (including lag parameter, 
#' test statistic and p.value) for each type of linear regression models of residuals 
#' \eqn{z[t]}. See \code{\link{adf.test}} for more details of three types of linear models.
#' 
#' @author Debin Qiu
#' @references
#' MacKinnon, J. G. (1991). Critical values for cointegration tests, Ch. 13 in Long-run
#' Economic Relationships: Readings in Cointegration, eds. R. F. Engle and C. W. J.
#' Granger, Oxford, Oxford University Press.
#' @seealso \code{\link{adf.test}}
#' @examples X <- matrix(rnorm(200),100,2)
#' y <- 0.3*X[,1] + 1.2*X[,2] + rnorm(100)
#' # test for original y and X
#' coint.test(y,X)
#' 
#' # test for response = diff(y,differences = 1) and
#' # input = apply(X, diff, differences = 1)
#' coint.test(y,X,d = 1)  
#' @importFrom stats lm
#' @importFrom stats residuals
#' @importFrom stats embed
#' @importFrom stats as.formula
#' @importFrom stats update
#' @importFrom stats approx
#' @export 
coint.test <- function(y,X,d = 0,nlag = NULL,output = TRUE)
{
  dname.y <- deparse(substitute(y))
  dname.X <- deparse(substitute(X))
  if (NCOL(y) > 1)
    stop("'y' must be a numeric vector or univariate time series")
  n.y <- length(y)
  n.X <- NROW(X)
  m.X <- NCOL(X)
  if (n.y != n.X)
    stop("'y' and 'X' must have same number of observations")
  if (d > n.y || d < 0 || d%%1 != 0)
    stop("'d' must be an integer between [0,length(y)]")
  nlag <- ifelse(is.null(nlag),floor(4*(length(y)/100)^(2/9)),nlag)
  if (d > 0) {
    y <- diff(y,d)
    X <- if (m.X > 1) apply(X,2,"diff",differences = d) 
         else diff(X,differences = d)
    dname.y <- paste("diff(",dname.y,",",d,")",sep = "")
    dname.X <- paste("diff(",dname.X,",",d,")",sep = "")
  }
  fmla <- paste("y ~ X",if (d > 0) "- 1" else "+ 1") 
  fit.lm <- lm(as.formula(fmla))
  res <- residuals(fit.lm) 
  yy <- diff(res)
  n.yy <- length(yy)  
  z <- embed(yy,nlag)
  yyt <- z[,1]
  t <- nlag:n.yy
  rt1 <- res[t]
  m1 <- lm(yyt ~ rt1 - 1)
  m2 <- lm(yyt ~ rt1 + t - 1)
  m3 <- lm(yyt ~ rt1 + t + I(t^2) - 1)
  if (nlag > 1) {
    yyt1 <- z[,2:nlag]
    m1 <- update(m1,.~. + yyt1)
    m2 <- update(m2,.~. + yyt1)
    m3 <- update(m3,.~. + yyt1)
  }
  STAT <- c(summary(m1)$coefficients[1,1]/summary(m1)$coefficients[1,2],
            summary(m2)$coefficients[2,1]/summary(m2)$coefficients[2,2],
            summary(m3)$coefficients[2,1]/summary(m3)$coefficients[2,2])
  EG_Tab <- matrix(c(1,0.01,-3.4304,-6.5393,-16.786,-79.433,
                 1,0.05,-2.8615,-2.8903,-4.234,-40.04,
                 1,0.1,-2.5668,-1.5384,-2.809,0,
                 2,0.01,-3.8964,-10.9519,-22.527,0,
                 2,0.05,-3.3361,-6.1101,-6.823,0,
                 2,0.1,-3.0445,-4.2412,-2.72,0,
                 3,0.01,-4.2937,-14.4354,-33.195,47.433,
                 3,0.05,-3.7407,-8.5631,-10.852,27.982,
                 3,0.1,-3.4522,-6.2143,-3.718,0,
                 4,0.01,-4.6433,-18.1031,-37.972,0,
                 4,0.05,-4.096,-11.2349,-11.175,0,
                 4,0.1,-3.8102,-8.3931,-4.137,0,
                 5,0.01,-4.9576,-21.8883,-45.142,0,
                 5,0.05,-4.4152,-14.0406,-12.575,0,
                 5,0.1,-4.1316,-10.7417,-3.784,0,
                 6,0.01,-5.2457,-25.6688,-57.737,88.639,
                 6,0.05,-4.7069,-16.9178,-17.492,60.007,
                 6,0.1,-4.425,-13.1875,-5.104,27.877,
                 7,0.01,-5.5123,-29.576,-69.398,164.295,
                 7,0.05,-4.9768,-19.9021,-22.045,110.761,
                 7,0.1,-4.6965,-15.7315,-6.922,67.721,
                 8,0.01,-5.762,-33.5258,-82.189,256.289,
                 8,0.05,-5.2292,-23.0023,-24.646,144.479,
                 8,0.1,-4.9501,-18.3959,-7.344,94.872,
                 9,0.01,-5.9974,-37.6572,-87.365,248.316,
                 9,0.05,-5.467,-26.2057,-26.627,176.382,
                 9,0.1,-5.189,-21.1377,-9.484,172.704,
                 10,0.01,-6.221,-41.7154,-102.68,389.33,
                 10,0.05,-5.6924,-29.4521,-30.994,251.016,
                 10,0.1,-5.4153,-24.0006,-7.514,163.049,
                 11,0.01,-6.4338,-46.0084,-106.809,352.752,
                 11,0.05,-5.9071,-32.8336,-30.275,249.994,
                 11,0.1,-5.6309,-26.9693,-4.083,151.427,
                 12,0.01,-6.6379,-50.2095,-124.156,579.622,
                 12,0.05,-6.1128,-36.2681,-32.505,314.802,
                 12,0.1,-5.8372,-29.9864,-2.686,184.116,
                 1,0.01,-3.9588,-9.0531,-28.428,-134.155,
                 1,0.05,-3.4105,-4.3904,-9.036,-45.374,
                 1,0.1,-3.1271,-2.5856,-3.925,-22.38,
                 2,0.01,-4.3276,-15.4387,-35.679,0,
                 2,0.05,-3.7806,-9.5106,-12.074,0,
                 2,0.1,-3.4963,-7.0815,-7.538,21.892,
                 3,0.01,-4.6631,-18.7688,-49.793,104.244,
                 3,0.05,-4.1189,-11.8922,-19.031,77.332,
                 3,0.1,-3.8351,-9.0723,-8.504,35.403,
                 4,0.01,-4.9694,-22.4694,-52.599,51.314,
                 4,0.05,-4.4287,-14.5876,-18.228,39.647,
                 4,0.1,-4.1463,-11.25,-9.873,54.109,
                 5,0.01,-5.2528,-26.2183,-59.631,50.646,
                 5,0.05,-4.7154,-17.3569,-22.66,91.359,
                 5,0.1,-4.4342,-13.6078,-10.238,76.781,
                 6,0.01,-5.5173,-29.976,-75.222,202.253,
                 6,0.05,-4.9823,-20.305,-25.224,132.03,
                 6,0.1,-4.7023,-16.1253,-9.836,94.272,
                 7,0.01,-5.7654,-33.9165,-84.312,245.394,
                 7,0.05,-5.233,-23.3328,-28.955,182.342,
                 7,0.1,-4.9541,-18.7352,-10.168,120.575,
                 8,0.01,-6,-37.8892,-96.428,335.92,
                 8,0.05,-5.4697,-26.4771,-31.034,220.165,
                 8,0.1,-5.1918,-21.4328,-10.726,157.955,
                 9,0.01,-6.2229,-41.9496,-109.881,466.068,
                 9,0.05,-5.6945,-29.7152,-33.784,273.002,
                 9,0.1,-5.4174,-24.2882,-8.584,169.891,
                 10,0.01,-6.4355,-46.1151,-120.814,566.823,
                 10,0.05,-5.9089,-33.0251,-37.208,346.189,
                 10,0.1,-5.6326,-27.2042,-6.792,177.666,
                 11,0.01,-6.6389,-50.4287,-128.997,642.781,
                 11,0.05,-6.114,-36.461,-36.246,348.554,
                 11,0.1,-5.8385,-30.1995,-5.163,210.338,
                 12,0.01,-6.8349,-54.7119,-139.8,736.376,
                 12,0.05,-6.3113,-39.9676,-37.021,406.051,
                 12,0.1,-6.0365,-33.2381,-6.606,317.776,
                 1,0.01,-4.37113,-11.5882,-35.819,-334.047,
                 1,0.05,-3.83239,-5.9057,-12.49,-118.284,
                 1,0.1,-3.55326,-3.6596,-5.293,-63.559,
                 2,0.01,-4.69276,-20.2284,-64.919,88.884,
                 2,0.05,-4.15387,-13.3114,-28.402,72.741,
                 2,0.1,-3.87346,-10.4637,-17.408,66.313,
                 3,0.01,-4.99071,-23.5873,-76.924,184.782,
                 3,0.05,-4.45311,-15.7732,-32.316,122.705,
                 3,0.1,-4.1728,-12.4909,-17.912,83.285,
                 4,0.01,-5.2678,-27.2836,-78.971,137.871,
                 4,0.05,-4.73244,-18.4833,-31.875,111.817,
                 4,0.1,-4.45268,-14.7199,-17.969,101.92,
                 5,0.01,-5.52826,-30.9051,-92.49,248.096,
                 5,0.05,-4.99491,-21.236,-37.685,194.208,
                 5,0.1,-4.71587,-17.082,-18.631,136.672,
                 6,0.01,-5.77379,-34.701,-105.937,393.991,
                 6,0.05,-5.24217,-24.2177,-39.153,232.528,
                 6,0.1,-4.96397,-19.6064,-18.858,174.919,
                 7,0.01,-6.00609,-38.7383,-108.605,365.208,
                 7,0.05,-5.47664,-27.3005,-39.498,246.918,
                 7,0.1,-5.19921,-22.2617,-17.91,208.494,
                 8,0.01,-6.22758,-42.7154,-119.622,421.395,
                 8,0.05,-5.69983,-30.4365,-44.3,345.48,
                 8,0.1,-5.4232,-24.9686,-19.688,274.462,
                 9,0.01,-6.43933,-46.7581,-136.691,651.38,
                 9,0.05,-5.91298,-33.7584,-42.686,346.629,
                 9,0.1,-5.63704,-27.8965,-13.88,236.975,
                 10,0.01,-6.64235,-50.9783,-145.462,752.228,
                 10,0.05,-6.11753,-37.056,-48.719,473.905,
                 10,0.1,-5.84215,-30.8119,-14.938,316.006,
                 11,0.01,-6.83743,-55.2861,-152.651,792.577,
                 11,0.05,-6.31396,-40.5507,-46.771,487.185,
                 11,0.1,-6.03921,-33.895,-9.122,285.164,
                 12,0.01,-7.02582,-59.6037,-166.368,989.879,
                 12,0.05,-6.50353,-44.0797,-47.242,543.889,
                 12,0.1,-6.22941,-36.9673,-10.868,418.414),nrow = 108,byrow = TRUE)
  EG.critical.value <- function(N,n,type = c("nt","lt","qt")) {
    type <- match.arg(type)
    C.val <- function(n,BETA) BETA[,1] + BETA[,2]/n + 
                              BETA[,3]/n^2 + BETA[,4]/n^3
    T_EG <- switch(type,nt = EG_Tab[1:36,],lt = EG_Tab[37:72,],qt = EG_Tab[73:108,])
    Tab <- cbind(T_EG[,1:2],C.val(n,T_EG[,-c(1:2)]))
    t.N <- Tab[Tab[,1] == N,]
    names(t.N) <- c("N","P","C")
    return(t.N)            
  }
  table1 <- EG.critical.value(m.X,n.y,type = "nt")
  table2 <- EG.critical.value(m.X,n.y,type = "lt")
  table3 <- EG.critical.value(m.X,n.y,type = "qt")
  PVAL <- c(approx(table1[,3],table1[,2],STAT[1],rule = 2)$y,
            approx(table2[,3],table2[,2],STAT[2],rule = 2)$y,
            approx(table3[,3],table3[,2],STAT[3],rule = 2)$y) 
  result <- matrix(c(rep(nlag,3),STAT,PVAL),3,3,dimnames = 
                     list(paste("type",1:3),c("lag","EG","p.value")))
  if (output) {
    cat("Response:",dname.y,"\n")
    cat("Input:",dname.X,"\n")
    cat("Number of inputs:",m.X,"\n")
    cat("Model:",fmla,"\n")
    cat("-------------------------------","\n")
    cat("Engle-Granger Cointegration Test","\n")
    cat("alternative: cointegrated","\n\n") 
    cat("Type 1: no trend","\n")
    print(result[1,], digits = 3) 
    cat("-----","\n","Type 2: linear trend","\n")
    print(result[2,], digits = 3)
    cat("-----","\n","Type 3: quadratic trend","\n")
    print(result[3,], digits = 3)
    cat("-----------","\n")
    cat("Note: p.value = 0.01 means p.value <= 0.01","\n")
    cat("    : p.value = 0.10 means p.value >= 0.10","\n")
  }
  coint.test <- result
}