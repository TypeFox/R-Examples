##' Rmixmod a MIXture MODelling package
##'
##' Rmixmod is a package based on the existing MIXMOD software. MIXMOD is a tool for fitting a mixture model of multivariate gaussian or multinomial components to a given data set with either a clustering, a density estimation or a discriminant analysis point of view.
##'
##' \tabular{ll}{
##'   Package: \tab Rmixmod\cr 
##'   Type: \tab Package\cr 
##'   Version: \tab 2.0.3\cr
##'   Date: \tab 2015-10-07\cr 
##'   License: \tab GPL-3 + file LICENSE\cr 
##'   LazyLoad: \tab yes\cr
##' }
##'
##' The general purpose of the package is to discover, or explain, group structures in multivariate data sets with unknown (cluster analysis or clustering) or known class discriminant analysis or classification). It is an exploratory data analysis tool for solving clustering and classification problems. But it can also be regarded as a semi-parametric tool to estimate densities with Gaussian mixture distributions and multinomial distributions.
##'
##' Mathematically, mixture probability density function (pdf) \eqn{f} is a weighted sum of \eqn{K} components densities : 
##'
##' \deqn{
##'   f({\bf x}_i|\theta) = \sum_{k=1}^{K}p_kh({\bf x}_i|\lambda_k)
##' }
##' where \eqn{h(.|{\lambda}_k)} denotes a \eqn{d}-dimensional distribution parametrized by \eqn{\lambda_k}. The parameters are the mixing proportions \eqn{p_k} and the component of the distribution \eqn{\lambda_k}.\cr
##'
##' In the Gaussian case, \eqn{h} is the density of a Gaussian distribution with mean \eqn{\mu_k} and variance matrix \eqn{\Sigma_k}, and thus \eqn{\lambda_k = (\mu_k,\Sigma_k)}.
##'
##' In the qualitative case, \eqn{h} is a multinomial distribution and \eqn{\lambda_k=(a_k,\epsilon_k)} is the parameter of the distribution.
##'
##' Estimation of the mixture parameters is performed either through maximum likelihood via the EM (\emph{Expectation Maximization}, Dempster et al. 1977), the SEM (\emph{Stochastic EM}, Celeux and Diebolt 1985) algorithm or through classification maximum likelihood via the CEM  algorithm (\emph{Clustering EM}, Celeux and Govaert 1992). These three algorithms can be chained to obtain original fitting strategies (e.g. CEM then EM with results of CEM) to use advantages of each of them in the estimation process. As mixture problems usually have multiple relative maxima, the program will produce different results, depending on the initial estimates supplied by the user. If the user does not input his own initial estimates, some initial estimates procedures are proposed (random centers for instance).
##'
##' It is possible to constrain some input parameters. For example, dispersions can be equal between classes, etc.
##'
##' In the Gaussian case, fourteen models are implemented. They are based on the eigenvalue decomposition, are most generally used. They depend on constraints on the variance matrix such as same variance matrix between clusters, spherical variance matrix... and they are suitable for data sets in any dimension.
##'
##' In the qualitative case, five multinomial models are available. They are based on a reparametrization of the multinomial probabilities. 
##'
##' In both cases, the models and the number of clusters can be chosen by different criteria : BIC (Bayesian Information Criterion), ICL (Integrated Completed Likelihood, a classification version of BIC), NEC (Entropy Criterion), or Cross-Validation (CV).
##'
##'
##' @name Rmixmod-package
##' @aliases Rmixmod
##' @rdname Rmixmod-package
##' @docType package
##' @keywords package
##' @import Rcpp
##'
##' @author
##' Author: Remi Lebret and Serge Iovleff and Florent Langrognet, with contributions from C. Biernacki and G. Celeux and G. Govaert \email{contact@@mixmod.org}
##'
##' @references 
##'   R. Lebret, S. Iovleff, F. Langrognet, C. Biernacki, G. Celeux, G. Govaert (2015), "Rmixmod: The R Package of the Model-Based Unsupervised, Supervised, and Semi-Supervised Classification Mixmod Library", Journal of Statistical Software, 67(6), 1-29, doi:10.18637/jss.v067.i06
##'
##' @examples
##'   \dontrun{
##'   ## Clustering Analysis
##'   # load quantitative data set
##'   data(geyser)
##'   # Clustering in gaussian case
##'   xem1<-mixmodCluster(geyser,3)
##'   summary(xem1)
##'   plot(xem1)
##'   hist(xem1)
##'
##'   # load qualitative data set
##'   data(birds)
##'   # Clustering in multinomial case
##'   xem2<-mixmodCluster(birds, 2)
##'   summary(xem2)
##'   barplot(xem2)
##'
##'   # load heterogeneous data set
##'   data(finance)
##'   # Clustering in composite case
##'   xem3<-mixmodCluster(finance,2:6)
##'   summary(xem3)
##'
##'   ## Discriminant Analysis
##'   # start by extract 10 observations from iris data set
##'   remaining.obs<-sample(1:nrow(iris),10)
##'   # then run a mixmodLearn() analysis without those 10 observations
##'   learn<-mixmodLearn(iris[-remaining.obs,1:4], iris$Species[-remaining.obs])
##'   # create a MixmodPredict to predict those 10 observations
##'   prediction <- mixmodPredict(data=iris[remaining.obs,1:4], classificationRule=learn["bestResult"])
##'   # show results
##'   prediction
##'   # compare prediction with real results
##'   paste("accuracy= ",mean(as.integer(iris$Species[remaining.obs]) == prediction["partition"])*100
##'      	,"%",sep="")
##'   }
##'
##' @useDynLib Rmixmod
##' @exportPattern "^[[:alpha:]]+"
NULL

##' Quantitative data: Old Faithful Geyser
##' 
##' The file geyser.rda contains 272 observations from the Old Faithful Geyser in the Yellowstone National Park. Each observation consists of two measurements : the duration (in minutes) of the eruption and the waiting time (in minutes) to the next eruption.
##'
##' Old Faithful erupts more frequently than any other big geyser, although it is not the largest nor the most regular geyser in the park. Its average interval between two eruptions is about 76 minutes, varying from 45 - 110 minutes. An eruption lasts from 1.1/2 to 5 minutes, expels 3,700 - 8,400 gallons (14,000 - 32,000 liters) of boiling water, and reaches heights of 106 - 184 feet (30 - 55m). It was named for its consistent performance by members of the Washburn Expedition in 1870. Old Faithful is still as spectacular and predictable as it was a century ago.
##'
##' @format A data frame with 272 observations on the following 2 variables.
##'
##' \describe{
##'
##'   \item{\code{Duration}}{a numeric vector containing the duration (in minutes) of the eruption}
##'
##'   \item{\code{Waiting.Time}}{a numeric vector containing the waiting time (in minutes) to the next eruption}
##'
##' }
##'
##' @source \url{http://www.geyserstudy.org/geyser.aspx?pGeyserNo=OLDFAITHFUL}
##'
##' @references 
##' Hardle, W. (1991). "Smoothing Techniques with Implementation in S". Springer-Verlag, New York. 
##' Azzalini, A. and Bowman, A. W. (1990). "A look at some data on the Old Faithful geyser". Applied Statistics 39, 357-365.
##'
##' @name geyser
##' @docType data
##' @keywords datasets
##' 
##' @examples
##'   data(geyser)
NULL

##' Qualitative data: Survival of passengers on the Titanic
##' 
##' For each person on board the fatal maiden voyage of the ocean liner Titanic, this dataset records: sex, age [adult/child], economic status [first/second/third class, or crew] and whether or not that person survived.
##' Values are aligned and delimited by blanks. There are no missing values.
##'
##' The sinking of the Titanic is a famous event, and new books are still being published about it.  Many well-known facts-from the proportions of first-class passengers to the "women and children first" policy, and the fact that that policy was not entirely successful in saving the women and children in the third class-are reflected in the survival rates for various classes of passenger.
##'
##' These data were originally collected by the British Board of Trade in their investigation of the sinking.  Note that there is not complete agreement among primary sources as to the exact numbers on board, rescued, or lost.
##'
##' Due in particular to the very successful film "Titanic", the last years saw a rise in public interest in the Titanic.  Very detailed data about the passengers is now available on the Internet, at sites such as "Encyclopedia Titanica" (\url{http://www.rmplc.co.uk/eduweb/sites/phind}).
##'
##' @format A data frame with 2201 observations on the following 4 variables.
##'
##' \describe{
##'
##'   \item{\code{Class}}{0 = crew, 1 = first, 2 = second, 3 = third, which denote the economic status of the subject}
##'
##'   \item{\code{Age}}{1 = adult, 0 = child, which denote if the subject is an adult or a child}
##'
##'   \item{\code{Sex}}{1 = male, 0 = female, which denote the sex of the subject}
##'
##'   \item{\code{Survived}}{1 = yes, 0 = no, which denote if the subject lived through the fatal maiden voyage of the ocean liner Titanic} 
##'
##' }
##'
##' @source 
##' The source provides a data set recording class, sex, age, and survival status for each person on board of the Titanic, and is based on data originally collected by the British Board of Trade and reprinted in:
##' British Board of Trade (1990), "Report on the Loss of the Titanic (S.S.)".  British Board of Trade Inquiry Report (reprint).  Gloucester, UK: Allan Sutton Publishing.
##'
##' @name titanic
##' @docType data
##' @keywords datasets
##' 
##' @examples
##'   data(titanic)
NULL

##' Qualitative data : morphological description of birds
##' 
##' The dataset contains details on the morphology of birds (puffins). Each individual (bird) is described by 6 qualitative variables. One variable for the gender and 5 variables giving a morphological description of the birds. There is 69 puffins divided in 2 sub-classes: lherminieri (34) and subalaris (35).
##'
##' @format A data frame with 69 observations on the following 5 variables.
##'
##' \describe{
##'
##'   \item{\code{gender}}{a numeric vector defining the gender (2 modalities, male or female).}
##'
##'   \item{\code{eyebrow}}{a numeric vector describing the eyebrow stripe (4 modalities).}
##'
##'   \item{\code{collar}}{a numeric vector describing the collar (5 modalities).}
##'
##'   \item{\code{sub-caudal}}{a numeric vector describing the sub-caudal (5 modalities).}
##'
##'   \item{\code{border}}{a numeric vector describing the border (3 modalities).}
##'
##' }
##'
##' @source
##' Bretagnolle, V., 2007. Personal communication, source: Museum.
##'
##' @name birds
##' @docType data
##' @keywords datasets
##' 
##' @examples
##'   data(birds)
NULL


##' Qualitative data : Car Evaluation
##' 
##' Car Evaluation Database was derived from a simple hierarchical decision model originally developed for the demonstration of DEX, M. Bohanec, V. Rajkovic: Expert system for decision making.
##'
##' @format A data frame with 1728 observations on the following 6 variables.
##'
##' \describe{
##'
##'   \item{\code{buying}}{the buying price (4 modalities: vhigh, high, med, low).}
##'
##'   \item{\code{maint}}{the price of the maintenance (4 modalities: vhigh, high, med, low).}
##'
##'   \item{\code{doors}}{the number of doors (4 modalities: 2, 3, 4, 5more).}
##'
##'   \item{\code{persons}}{the capacity in terms of persons to carry (3 modalities: 2, 4, more).}
##'
##'   \item{\code{lug_boot}}{the size of luggage boot  (3 modalities: small, med, big).}
##'
##'   \item{\code{safety}}{the estimated safety of the car (3 modalities: low, med, high).}
##'
##'   \item{\code{acceptability}}{the car acceptability (4 modalities: unacc, acc, good, vgood).}
##' }
##'
##' @source
##' Creator: Marko Bohanec 
##' Donors: Marko Bohanec & Blaz Zupan
##' http://archive.ics.uci.edu/ml/datasets/Car+Evaluation
##'
##' @name car
##' @docType data
##' @keywords datasets
##' 
##' @examples
##'   data(car)
NULL


##' Composite data : Financial health of companies
##' 
##' This data set is made up of 216 healthy firms and 212 bankruptcy firms (year 2002) and also 241 healthy firms and 220 bankruptcy firms (year 2003). Companies are described by four financial ratios expected to provide some meaningful information about their health: EBITDA/Total Assets, Value Added/Total Sales, Quick Ratio, Accounts Payable/Total Sales. This data set offers the possibility to predict the company's ability to cover its financial obligations and also to study its stability over the years.
##'
##' @format A data frame with 889 companies (rows) and 6 variables (columns).
##'
##' \describe{
##'
##'   \item{\code{Year}}{categorical variable with two modalities (2002 & 2003).}
##'
##'   \item{\code{Health}}{categorical variable with two modalities (bankruptcy & healthy).}
##'
##'   \item{\code{EBITDA.Total.Assets}}{numeric variable.}
##'
##'   \item{\code{Value.Added.Total.Sales}}{numeric variable.}
##'
##'   \item{\code{Quick.Ratio}}{numeric variable.}
##'
##'   \item{\code{Accounts.Payable.Total.Sales}}{numeric variable.}
##' }
##'
##' @source
##' Lourme A, Biernacki C (2011). \emph{Simultaneous t-Model-Based Clustering for Data Differing over Time Period: Application for Understanding Companies Financial Health.} Case Studies in Business, Industry and Government Statistics, 4(2), 73-82.
##'
##' Du Jardin P, S\'everin E (2010). \emph{Dynamic analysis of the business failure process: a study of bankruptcy trajectories.} In Portuguese Finance Network. Ponte Delgada, Portugual.
##'
##' @name finance
##' @docType data
##' @keywords datasets
##' 
##' @examples
##' data(finance)
##' summary(finance)
NULL

##' Composite data with training and testing set
##' 
##' The data set is made up of 5 variables: 3 categorical variables and 2 quantitative variables.
##' The original data set contains 200 individuals. The training data set has 300 individuals while the testing data set has 100 individuals.
##'
##' @format A data frame with 200 individuals (rows) and 5 variables (columns).
##'
##' \describe{
##'
##'   \item{\code{V1}}{categorical variable with two modalities (1 & 2).}
##'
##'   \item{\code{V2}}{categorical variable with two modalities (1 & 2).}
##'
##'   \item{\code{V3}}{categorical variable with two modalities (1 & 2).}
##'
##'   \item{\code{V4}}{numeric variable.}
##'
##'   \item{\code{V5}}{numeric variable.}
##' }
##'
##' @name heterodata
##' @docType data
##' @keywords datasets
##' @seealso \code{\link{heterodatatrain}} and \code{\link{heterodatatest}}
##' 
##' @examples
##' data(heterodata)
##' summary(heterodata)
NULL


##' Composite data: A training set
##' 
##' The data set is made up of 5 variables: 3 categorical variables and 2 quantitative variables.
##' The training data set has 300 individuals.
##'
##' @format A data frame with 300 individuals (rows) and 5 variables (columns).
##'
##' \describe{
##'
##'   \item{\code{V1}}{categorical variable with two modalities (1 & 2).}
##'
##'   \item{\code{V2}}{categorical variable with two modalities (1 & 2).}
##'
##'   \item{\code{V3}}{categorical variable with two modalities (1 & 2).}
##'
##'   \item{\code{V4}}{numeric variable.}
##'
##'   \item{\code{V5}}{numeric variable.}
##' }
##'
##' @name heterodatatrain
##' @docType data
##' @keywords datasets
##' @seealso \code{\link{heterodatatest}}
##' 
##' @examples
##' data(heterodatatrain)
##' summary(heterodatatrain)
NULL

##' Composite data: A testing set
##' 
##' The data set is made up of 5 variables: 3 categorical variables and 2 quantitative variables.
##' The testing data set has 100 individuals.
##'
##' @format A data frame with 100 individuals (rows) and 5 variables (columns).
##'
##' \describe{
##'
##'   \item{\code{V1}}{categorical variable with two modalities (1 & 2).}
##'
##'   \item{\code{V2}}{categorical variable with two modalities (1 & 2).}
##'
##'   \item{\code{V3}}{categorical variable with two modalities (1 & 2).}
##'
##'   \item{\code{V4}}{numeric variable.}
##'
##'   \item{\code{V5}}{numeric variable.}
##' }
##'
##' @name heterodatatest
##' @docType data
##' @keywords datasets
##' @seealso \code{\link{heterodatatrain}}
##' 
##' @examples
##' data(heterodatatest)
##' summary(heterodatatest)
NULL
