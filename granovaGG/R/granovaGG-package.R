#' Elemental Graphics for Analysis of Variance Using ggplot2
#'
#' This collection of functions in granovaGG provides what we call
#' elemental graphics for display of anova results. The term
#' elemental derives from the fact that each function is
#' aimed at construction of graphical displays that afford
#' direct visualizations of data with respect to the
#' fundamental questions that drive the particular anova
#' methods. This package represents a modification of the
#' original granova package; the key change is to use ggplot2,
#' Hadley Wickham's package based on Grammar of Graphics
#' concepts (due to Wilkinson). The main function is granovagg.1w
#' (a graphic for one way anova); two other functions (granovagg.ds
#' and granovagg.contr) are to construct graphics for dependent
#' sample analyses and contrast-based analyses respectively. (The
#' function granova.2w, which entails dynamic displays of data, is
#' not currently part of granovaGG.) The granovaGG functions are
#' to display data for any number of groups, regardless of
#' their sizes (however, very large data sets or numbers of
#' groups can be problematic). For granovagg.1w a
#' specialized approach is used to construct data-based
#' contrast vectors for which anova data are displayed. The
#' result is that the graphics use a straight line to facilitate clear
#' interpretations while being faithful to the standard
#' effect test in anova. The graphic results are
#' complementary to standard summary tables; indeed, numerical
#' summary statistics are provided as side effects of the graphic
#' constructions. granovagg.ds and granovagg.contr provide
#' graphic displays and numerical outputs for a dependent
#' sample and contrast-based analyses. The graphics based on
#' these functions can be especially helpful for learning
#' how the respective methods work to answer the basic
#' question(s) that drive the analyses. This means they can be particularly
#' helpful for students and non-statistician analysts. But
#' these methods can be of assistance for work-a-day
#' applications of many kinds, as they can help to
#' identify outliers, clusters or patterns, as well as
#' highlight the role of non-linear transformations of data.
#' In the case of granovagg.1w and granovagg.ds
#' several arguments are provided to facilitate flexibility
#' in the construction of graphics that accommodate diverse
#' features of data, according to their corresponding
#' display requirements. See the help files for individual
#' functions.
#'
#'
#' \tabular{ll}{
#'   Package: \tab granovaGG\cr
#'   Version: \tab 1.0\cr
#'   License: \tab GPL (>= 2)\cr
#' }
#'
#' @author Brian A. Danielak \email{brian@@briandk.com}\cr
#'   Robert M. Pruzek \email{RMPruzek@@yahoo.com}
#'
#' with contributions by:\cr
#'   William E. J. Doane \email{wil@@DrDoane.com}\cr
#'   James E. Helmreich \email{James.Helmreich@@Marist.edu}\cr
#'   Jason Bryer \email{jason@@bryer.org}
#'
#' @import ggplot2
#' @name granovaGG-package
#' @aliases granovaGG-package granovaGG
#' @docType package
#' @seealso
#'
#' \code{\link{granovagg.1w}} \code{\link{granovagg.ds}}
#'   \code{\link{granovagg.contr}}
#' @keywords hplot
#' @references Wickham, H. (2009). Ggplot2: Elegant Graphics for Data Analysis. New York: Springer.
#' @references Wilkinson, L. (1999). The Grammar of Graphics. Statistics and computing. New York: Springer.

NULL

#' Virus Preparation on Tobacco Leaves
#'
#' This data is taken from Snedecor and Cochran (1980) and corresponds to a true
#' matched pairs experiment. The data originally came from Youden and Beale in
#' 1934 who "wished to find out if two preparations of a virus would produce
#' different effects on tobacco plants. Half a leaf of a tobacco plant was
#' rubbed with cheesecloth soaked in one preparation of the virus extract, and
#' the second half was rubbed similarly with the second extract." (Page 86,
#' Snedecor and Cochran, 1980) Each of the 8 points in the figure corresponds to
#' the numbers of lesions on the two halves of one leaf with sides that had been
#' treated differently.
#'
#' @name tobacco
#' @docType data
#' @format A dataframe with 8 observations on the following 2 variables, no NAs
#'
#' \describe{
#'    \item{prep1}{Virus Preparation 1}
#'    \item{prep2}{Virus Preparation 2}
#' }
#'
#' @references Snedecor, W., Cochran, W. (1980). Statistical methods. Iowa State
#'   University Press, Ames Iowa, seventh edition.
#' @source Youden, W. J., Beale, H. P. (1934). A statistical study of the local
#'   lesion method for estimating tobacco mosaic virus. In Contributions from
#'   Boyce Thompson Institute 6, page 437.
#' @keywords datasets
NULL

#' Anorexia Data on Weight Change
#'
#' The anorexia data frame has 72 rows and 3 columns. Weight change data for young female anorexia patients.
#' @name anorexia
#' @docType data
#' @format A dataframe with 72 observations of three variables:
#' 
#' \describe{
#'
#' \item{Treat}{Factor of three levels: "\code{Cont}" (control), "\code{CBT}" (Cognitive Behavioural treatment) and "\code{FT}" (family treatment).}
#' 
#' \item{Prewt}{Pretreatment weight of subject, in pounds.}
#'
#' \item{Postwt}{Postreatment weight of subject, in pounds.}
#'}
#'
#' @references Venables, W. N. and Ripley, B. D. (2002) Modern Applied
#'   Statistics with S. Fourth edition. Springer.
#' @source Hand, D. J., Daly, F., McConway, K., Lunn, D. and Ostrowski, E. eds
#'   (1993) A Handbook of Small Data Sets. Chapman & Hall, Data set 285 (p.
#'   229)
#' @keywords datasets

NULL

#' Shoe wear data of Box, Hunter and Hunter
#' 
#' A list of two vectors, giving the wear of shoes of materials A and B for one foot each of ten boys.
#' 
#' @name shoes
#' @docType data
#' @references Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S}. Fourth edition. Springer.
#' @source G. E. P. Box, W. G. Hunter and J. S. Hunter (1978) \emph{Statistics for Experimenters}. Wiley, p. 100

NULL 

#' Family Treatment Weight change data for young female anorexia patients (subset).
#'
#' The MASS package includes the dataset \code{anorexia}, containing pre and
#' post treatment weights for young female anorexia patients.  This is a subset
#' of those data, containing only those patients who received Family Treatment.
#'
#'
#' @name anorexia.sub
#' @docType data
#' @format A dataframe with 17 observations on the following 2 variables, no
#'   NAs.
#'
#' \describe{
#'
#' \item{\code{Prewt}}{Pretreatment weight of subject, in pounds.}
#'
#' \item{\code{Postwt}}{Postreatment weight of subject, in pounds.}
#'
#' }
#' @references Venables, W. N. and Ripley, B. D. (2002) Modern Applied
#'   Statistics with S. Fourth edition. Springer.
#' @source Hand, D. J., Daly, F., McConway, K., Lunn, D. and Ostrowski, E. eds
#'   (1993) A Handbook of Small Data Sets. Chapman & Hall, Data set 285 (p.
#'   229)
#' @keywords datasets
NULL

#' Arousal in Rats
#'

#'
#' 40 rats were given divided randomly into four groups and assigned to one of
#' four treatments: placebo, drug A, drug B, or both drug A and drug B.
#' Response is a standard measure of physiological arousal.
#'
#'
#' @name arousal
#' @docType data
#' @format A data frame with 40 observations, 10 in each of 4 columns the
#'   corresponding to placebo, drug A, drug B and both drug A and drug B; no
#'   NAs.
#'
#' \describe{
#'
#' \item{Placebo}{Rats receiving a placebo treatment.}
#'
#' \item{Drug.A}{Rats receiving only drug A.}
#'
#' \item{Drug.B}{Rats receiving only drug B.}
#'
#' \item{Drug.A.B}{Rats receiving both drug A and drug B.}
#'
#' }
#' @source Richard Lowry. Concepts & Applications of Inferential Statistics.
#'   Vassar College, Poughkeepsie, N.Y., 2010,
#'   http://faculty.vassar.edu/lowry/webtext.html
#' @keywords datasets
NULL

#' Blood lead levels of lead workers' children matched with similar control
#' children.
#'

#'
#' Children of parents who had worked in a factory where lead was used in
#' making batteries were matched by age, exposure to traffic, and neighborhood
#' with children whose parents did not work in lead-related industries. Whole
#' blood was assessed for lead content yielding measurements in mg/dl
#'
#'
#' @name blood_lead
#' @docType data
#' @format A dataframe with 33 observations on the following 2 variables, no
#'   NAs.
#'
#' \describe{
#'
#' \item{Exposed}{Blood lead level of exposed child, mg/dl.}
#'
#' \item{Control}{Blood lead level of exposed child, mg/dl.}
#'
#' }
#' @references See discussion in Section 2.5 of Enhancing Dependent Sample
#'   Analyses with Graphics, Journal of Statistics Education Volume 17, Number
#'   1 (March 2009).
#' @source Morton, D., Saah, A., Silberg, S., Owens, W., Roberts, M., Saah, M.
#'   (1982). Lead absorption in children of employees in a lead related
#'   industry. American Journal of Epidemiology, 115:549-555.
#' @keywords datasets
NULL

#' Poison data from Biological Experiment
#'
#' Survial times of animals in a 3 x 4 factorial experiment involving poisons
#' (3 levels) and various treatments (four levels), as described in Chapter 8
#' of Box, Hunter and Hunter.
#'
#'
#' @name poison
#' @docType data
#' @format This data frame was originally \code{poison.data} from the package
#'   \code{BHH2}, but as presented here has added columns; no NAs. \describe{
#'   \item{Poison}{Factor with three levels I, II, and III.}
#'   \item{Treatment}{Factor with four levels, A, B, C, and D.}
#'   \item{Group}{Factor with 12 levels, 1:12.}
#'   \item{SurvTime}{Numeric; survival time.}
#'   \item{RateSurvTime}{Numeric; inverse of SurvTime}
#'   \item{RankRateSurvTime}{Numeric; \code{RateSurvTime} scores have
#'   been converted to ranks, and then rescaled to have the same median as and
#'   a spread comparable to \code{RateSurvTime}} }
#' @references Box G. E. P, Hunter, J. S. and Hunter, W. C. (2005). Statistics
#'   for Experimenters II. New York: Wiley.
#' @source Box, G. E. P. and D. R. Cox, An Analysis of Transformations (with
#'   discussion), Journal of the Royal Statistical Society, Series B, Vol. 26,
#'   No. 2, pp. 211 - 254.
#' @keywords datasets
NULL

#' Weight gains of rats fed different diets
#'

#'
#' 60 rats were fed varying diets to see which produced the greatest weight
#' gain.  Two diet factors were protein type: beef, pork, chicken and protein
#' level: high and low.
#'
#'
#' @name rat
#' @docType data
#' @format A data frame with 60 observations on the following 3 variables, no
#'   NAs.
#'
#' \describe{
#'
#' \item{Weight.Gain}{Weight gain (grams) of rats fed the diets.}
#'
#' \item{Diet.Amount}{Amount of protein in diet: 1 = High, 2 = Low.}
#'
#' \item{Diet.Type}{Type of protein in diet: 1 = Beef, 2 = Pork, 3 =
#'   Cereal.}
#'
#' }
#' @source Fundamentals of Exploratory Analysis of Variance, Hoaglin D.,
#'   Mosteller F. and Tukey J. eds., Wiley, 1991, p. 100; originally from
#'   Statistical Methods, 7th ed, Snedecor G. and Cochran W. (1980), Iowa State
#'   Press.
#' @keywords datasets
NULL
