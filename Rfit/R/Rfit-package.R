

#' All Scores
#' 
#' An object of class scores which includes the score function and it's
#' derivative for rank-based regression inference.
#' 
#' Using Wilcoxon (linear) scores leads to inference which has ARE of 0.955 to
#' least squares (ML) when the data are normal. Wilcoxon scores are optimal
#' when the underlying error distribution is logistic. Normal scores are
#' optimal when the data are normally distributed. Log-rank scores are optimal
#' when the data are from an exponential distribution, e.g. in a proportional
#' hazards model. Log-Generalized F scores can also be used in the analysis of
#' survival data (see Hettmansperger and McKean p. 233).
#' 
#' @name allscores
#' @aliases wscores nscores bentscores1 bentscores2 bentscores3 bentscores4
#' logGFscores logrank.scores
#' @docType data
#' @format The format is: Formal class 'scores' [package ".GlobalEnv"] with 2
#' slots ..@ phi :function (u) ..@ Dphi:function (u)
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' @keywords datasets
#' @examples
#' 
#' data(wscores)
#' x<-runif(10)
#' y<-rlogis(10)
#' rfit(y~x,scores=wscores)
#' 
NULL





#' The World Famous Baseball Data
#' 
#' These data come from the back-side of 59 baseball cards that Carrie had.
#' 
#' 
#' @name baseball
#' @docType data
#' @format A data frame with 59 observations on the following 6 variables.
#' \describe{ \item{list("height")}{Height in inches}
#' \item{list("weight")}{Weight in pounds} \item{list("bat")}{a factor with
#' levels \code{L} \code{R} \code{S}} \item{list("throw")}{a factor with levels
#' \code{L} \code{R}} \item{list("field")}{a factor with levels \code{0}
#' \code{1}} \item{list("average")}{ERA if the player is a pitcher and his
#' batting average if the player is a fielder} }
#' @source Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' @keywords datasets
#' @examples
#' 
#' data(baseball)
#' wilcox.test(height~field,data=baseball)
#' rfit(weight~height,data=baseball)
#' 
NULL





#' Baseball Salaries
#' 
#' Salaries of 176 professional baseball players for the 1987 season.
#' 
#' 
#' @name bbsalaries
#' @docType data
#' @format A data frame with 176 observations on the following 8 variables.
#' \describe{ \item{list("logYears")}{Log of the number of years experience}
#' \item{list("aveWins")}{Average wins per year}
#' \item{list("aveLosses")}{Average losses per year} \item{list("era")}{Earned
#' Run Average} \item{list("aveGames")}{Average games pitched in per year}
#' \item{list("aveInnings")}{Average number of innings pitched per year}
#' \item{list("aveSaves")}{Average number of saves per year}
#' \item{list("logSalary")}{Log of the base salary in dollars} }
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' @source http://lib.stat.cmu.edu/datasets/baseball.data
#' @keywords datasets
#' @examples
#' 
#' data(bbsalaries)
#' summary(rfit(logSalary~logYears+aveWins+aveLosses+era+aveGames+aveInnings+aveSaves,data=bbsalaries))
#' 
NULL





#' Box and Cox (1964) data.
#' 
#' The data are the results of a 3 * 4 two-way design, where forty-eight
#' animals were exposed to three different poisons and four different
#' treatments. The design is balanced with four replications per cell. The
#' response was the log survival time of the animal.
#' 
#' 
#' @name BoxCox
#' @docType data
#' @format A data frame with 48 observations on the following 3 variables.
#' \describe{ \item{list("logSurv")}{log Survival Time} \item{list("Poison")}{a
#' factor indicating poison level} \item{list("Treatment")}{a factor indicating
#' treatment level} }
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' @source Box, G.E.P. and Cox, D.R. (1964), An analysis of transformations,
#' \emph{ Journal of the Royal Statistical Society, Series B, Methodological},
#' 26, 211-252.
#' @keywords datasets
#' @examples
#' 
#' data(BoxCox)
#' with(BoxCox,interaction.plot(Treatment,Poison,logSurv,median))
#' raov(logSurv~Poison+Treatment,data=BoxCox)
#' 
NULL





#' Cardiovascular risk factors
#' 
#' Data from a study to investigate assocation between uric acid and various
#' cardiovascular risk factors in developing countries (Heritier et. al. 2009).
#' There are 474 men and 524 women aged 25-64.
#' 
#' Data set and description taken from Heritier et. al. (2009) (c.f. Conen et.
#' al. 2004). Some not discussed (in Section 3.5) of the text and their
#' persummed meaning is listed followed by (?).
#' 
#' @name CardioRiskFactors
#' @docType data
#' @format A data frame with 998 observations on the following 14 variables.
#' \describe{ \item{list("age")}{Age of subject} \item{list("bmi")}{Body Mass
#' Index} \item{list("waisthip")}{waist/hip ratio(?)}
#' \item{list("smok")}{indicator for regular smoker}
#' \item{list("choles")}{total cholesterol} \item{list("trig")}{triglycerides
#' level in body fat} \item{list("hdl")}{high-density lipoprotien(?)}
#' \item{list("ldl")}{low-density lipoprotein} \item{list("sys")}{systolic
#' blood pressure} \item{list("dia")}{diastolic blood pressure(?)}
#' \item{list("Uric")}{serum uric} \item{list("sex")}{indicator for male}
#' \item{list("alco")}{alcohol intake (mL/day)} \item{list("apoa")}{apoprotein
#' A} }
#' @source Heritier, S., Cantoni, E., Copt, S., and Victoria-Feser, M. (2009),
#' \emph{Robust Methods in Biostatistics}, New York: John Wiley \& Sons.
#' 
#' Conen, D., Wietlisbach, V., Bovet, P., Shamlaye, C., Riesen, W., Paccaud,
#' F., and Burnier, M. (2004), Prevalence of hyperuricemia and relation of
#' serum uric acid with cardiovascular risk factors in a developing country.
#' \emph{BMC Public Health}, \url{http://www.biomedcentral.com/1471-2458/4/9}.
#' @keywords datasets
#' @examples
#' 
#' data(CardioRiskFactors)
#' fitF<-rfit(Uric~bmi+sys+choles+ldl+sex+smok+alco+apoa+trig+age,data=CardioRiskFactors)
#' fitR<-rfit(Uric~bmi+sys+choles+ldl+sex,data=CardioRiskFactors)
#' drop.test(fitF,fitR)
#' summary(fitR)
#' 
NULL





#' Free Fatty Acid Data
#' 
#' The response variable is level of free fatty acid in a sample of
#' prepubescent boys. The explanatory variables are age (in months), weight (in
#' lbs), and skin fold thickness.
#' 
#' 
#' @name ffa
#' @docType data
#' @format A data frame with 41 rows and 4 columns. \describe{
#' \item{list("age")}{ age in years } \item{list("weight")}{ weight in lbs }
#' \item{list("skin")}{ skin fold thinkness } \item{list("ffa")}{ free fatty
#' acid } }
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' @source Morrison, D.F. (1983), \emph{Applied Linear Statistical Models},
#' Englewood Cliffs, NJ:Prentice Hall.
#' @keywords datasets
#' @examples
#' 
#' data(ffa)
#' summary(rfit(ffa~age+weight+skin,data=ffa))  #using the default (Wilcoxon scores)
#' summary(rfit(ffa~age+weight+skin,data=ffa,scores=bentscores1))
#' 
NULL





#' ~~ Methods for Function getScores ~~
#' 
#' ~~ Methods for function \code{getScores} ~~
#' 
#' 
#' @name getScores-methods
#' @aliases getScores-methods getScores,scores-method getScores
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(object = \"scores\")")}{ %% ~~describe this method
#' here~~ } }
#' @keywords methods ~~ other possible keyword(s) ~~
NULL





#' ~~ Methods for Function getScoresDeriv ~~
#' 
#' ~~ Methods for function \code{getScoresDeriv} ~~
#' 
#' 
#' @name getScoresDeriv-methods
#' @aliases getScoresDeriv-methods getScoresDeriv,scores-method getScoresDeriv
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(object = \"scores\")")}{ %% ~~describe this method
#' here~~ } }
#' @keywords methods ~~ other possible keyword(s) ~~
NULL





#' Class "param"
#' 
#' Internal class for use with score functions.
#' 
#' 
#' @name param-class
#' @aliases param-class param
#' @docType class
#' @section Objects from the Class: A virtual Class: No objects may be created
#' from it.
#' @author John Kloke
#' @seealso \code{\linkS4class{scores}}
#' @keywords classes
#' @examples
#' 
#' showClass("param")
#' 
NULL





#' Quail Data
#' 
#' Thirty-nine quail were randomized to one of for treatments for lowering
#' cholesterol.
#' 
#' 
#' @name quail
#' @docType data
#' @format A data frame with 39 observations on the following 2 variables.
#' \describe{ \item{list("treat")}{a factor with levels \code{1} \code{2}
#' \code{3} \code{4}} \item{list("ldl")}{a numeric vector} }
#' @source Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' @keywords datasets
#' @examples
#' 
#' data(quail)
#' boxplot(ldl~treat,data=quail)
#' 
NULL





#' Rank-Based Estimates and Inference for Linear Models
#' 
#' Package provides functions for rank-based analyses of linear models.
#' Rank-based estimation and inference offers a robust alternative to least
#' squares.
#' 
#' \tabular{ll}{ Package: \tab Rfit\cr Type: \tab Package\cr Version: \tab
#' 0.18\cr Date: \tab 2014-06-25\cr License: \tab GPL (version 2 or later)\cr
#' LazyLoad: \tab yes\cr } %~~ An overview of how to use the package, including
#' the most important ~~ %~~ functions ~~
#' 
#' @name Rfit-package
#' @aliases Rfit-package Rfit
#' @docType package
#' @author John Kloke, Joesph McKean
#' 
#' Maintainer: John Kloke <kloke@@biostat.wisc.edu>
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' 
#' Jaeckel, L. A. (1972). Estimating regression coefficients by minimizing the
#' dispersion of residuals. \emph{Annal s of Mathematical Statistics}, 43, 1449
#' - 1458.
#' 
#' Jureckova, J. (1971). Nonparametric estimate of regression coefficients.
#' \emph{Annals of Mathematical Statistics }, 42, 1328 - 1338.
#' @keywords nonparametric robust regression package
#' @examples
#' 
#' data(baseball)
#' data(wscores)
#' fit<-rfit(weight~height,data=baseball)
#' summary(fit)
#' plot(fitted(fit),rstudent(fit))
#' 
#' ### Example of the Reduction (Drop) in dispersion test ###
#' y<-rnorm(47)
#' x1<-rnorm(47)
#' x2<-rnorm(47)
#' fitF<-rfit(y~x1+x2)
#' fitR<-rfit(y~x1)
#' drop.test(fitF,fitR)
#' 
#' 
#' 
#' 
NULL





#' Class "scores"
#' 
#' A score function and it's corresponding derivative is required for
#' rank-based estimation.  This object puts them together.
#' 
#' 
#' @name scores-class
#' @aliases scores-class scores
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("scores", ...)}.
#' @author John Kloke \email{kloke@@biostat.wisc.edu}
#' @seealso \code{\linkS4class{param}} %% ~~objects to See Also as
#' \code{\link{~~fun~~}}, ~~~ %% ~~or \code{\linkS4class{CLASSNAME}} for links
#' to other classes ~~~
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' @keywords classes
#' @examples
#' 
#' showClass("scores")
#' 
NULL





#' Serum Level of luteinizing hormone (LH)
#' 
#' Hollander and Wolfe (1999) discuss a 2 by 5 factorial design for a study to
#' determine the effect of light on the release of luteinizing hormone (LH).
#' The factors in the design are: light regimes at two levels (constant light
#' and 14 hours of light followed by 10 hours of darkness) and a luteinizing
#' release factor (LRF) at 5 different dosage levels. The response is the level
#' of luteinizing hormone (LH), nanograms per ml of serum in blood samples.
#' Sixty rats were put on test under these 10 treatment combinations, six rats
#' per combination.
#' 
#' 
#' @name serumLH
#' @docType data
#' @format A data frame with 60 observations on the following 3 variables.
#' \describe{ \item{list("serum")}{a numeric vector}
#' \item{list("light.regime")}{a factor with levels \code{Constant}
#' \code{Intermittent}} \item{list("LRF.dose")}{a factor with levels \code{0}
#' \code{10} \code{1250} \code{250} \code{50}} }
#' @references Hollander, M. and Wolfe, D.A. (1999), \emph{Nonparametric
#' Statistical Methods}, New York: Wiley.
#' @source Hollander, M. and Wolfe, D.A. (1999), \emph{Nonparametric
#' Statistical Methods}, New York: Wiley.
#' @keywords datasets
#' @examples
#' 
#' data(serumLH)
#' raov(serum~light.regime + LRF.dose + light.regime*LRF.dose, data = serumLH)
#' 
NULL





#' Internal Functions for Estimating tau
#' 
#' These are internal functions used for calculating the scale parameter tau
#' necessary for estimating the standard errors of coefficients for
#' rank-regression.
#' 
#' 
#' @aliases hstarreadyscr hstar looptau pairup
#' @param ehat Full model residals
#' @param delta Window parameter (proportion) used in the Koul et al. estimator
#' of tau.  Default value is 0.80.  If the ratio of sample size to number of
#' regression parameters (n to p) is less than 5, larger values such as 0.90 to
#' 0.95 are more approporiate.
#' @param y Argument of function hstar
#' @param abdord Ordered absolute differences of residuals
#' @param wtord Standardized (by const) ordered absolute differences of
#' residuals
#' @param const Range of score function
#' @param n Sample size
#' @param x Argument for pairup
#' @param type Argument for the function pairup
#' @param asc scores
#' @param ascpr derivative of the scores
#' @author Joseph McKean, John Kloke
#' @seealso \code{\link{gettau}}, \code{\link{rfit}}
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' 
#' Koul, H.L., Sievers, G.L., and McKean, J.W. (1987) An esimator of the scale
#' parameter for the rank analysis of linear models under general score
#' functions, \emph{Scandinavian Journal of Statistics}, 14, 131-141.
NULL





#' Telephone Data
#' 
#' The number of telephone calls (in tens of millions) made in Belgium from
#' 1950-1973.
#' 
#' 
#' @name telephone
#' @docType data
#' @format A data frame with 24 observations on the following 2 variables.
#' \describe{ \item{list("year")}{years since 1950 AD}
#' \item{list("calls")}{number of telephone calls in tens of millions} }
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' @source Rousseeuw, P.J. and Leroy, A.M. (1987), \emph{Robust Regression and
#' Outlier Detection}, New York: Wiley.
#' @keywords datasets
#' @examples
#' 
#' data(telephone)
#' plot(telephone)
#' abline(rfit(calls~year,data=telephone))
#' 
NULL



