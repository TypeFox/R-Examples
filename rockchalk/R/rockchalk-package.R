##' rockchalk: regression functions
##'
##'
##' Includes an ever-growing collection of functions that assist in
##' the presentation of regression models.  The initial function was
##' \code{\link{outreg}}, which produces LaTeX tables that summarize
##' one or many fitted regression models.  It also offers plotting
##' conveniences like \code{\link{plotPlane}} and
##' \code{\link{plotSlopes}}, which illustrate some of the variables
##' from a fitted regression model. For a detailed check on
##' multicollinearity, see \code{\link{mcDiagnose}}.  The user should
##' be aware of this fact: Not all of these functions lead to models
##' or types of analysis that we endorese.  Rather, they all lead to
##' analysis that is endorsed by some scholars, and we feel it is
##' important to facilitate the comparison of competing methods.  For
##' example, the function \code{\link{standardize}} will calculate
##' standardized regression coefficients for all predictors in a
##' regression model's design matrix in order to replicate results
##' from other statistical frameworks, no matter how unwise the use of
##' such coefficients might be. The function \code{\link{meanCenter}}
##' will allow the user to more selectively choose variables for
##' centering (and possibly standardization) before they are entered
##' into the design matrix.  Because of the importance of interaction
##' variables in regression analysis, the \code{\link{residualCenter}}
##' and \code{\link{meanCenter}} functions are offered.  While mean
##' centering does not actually help with multicollinearity of
##' interactive terms, many scholars have argued that it does.  The
##' meanCenter function can be compared with the "residual centering"
##' of interaction terms.
##'
##' @name rockchalk-package
##' @aliases rockchalk-package rockchalk
##' @docType package
##' @author Paul E. Johnson \email{pauljohn@@ku.edu}
##'
##' @references http://pj.freefaculty.org/R
##' @keywords regression hplot
NULL






##' Religious beliefs and crime rates
##'
##' The data national-level summary indicators of public opinion about
##' the existence of heaven and hell as well as the national rate of
##' violent crime.
##' @name religioncrime
##' @docType data
##' @usage data(religioncrime)
##' @author Paul E. Johnson \email{pauljohn@@ku.edu} and Anonymous
##' @format data.frame: 51 obs. of 3 variables
##' @keywords datasets
##' @source  Anonymous researcher who claims the data is real.
##' @examples
##' require(rockchalk)
##' data(religioncrime)
##' mod1 <- lm(crime ~ heaven, data=religioncrime)
##' mod2 <- lm(crime ~ hell, data=religioncrime)
##' mod3 <- lm(crime ~ heaven + hell, data=religioncrime)
##' with(religioncrime,
##' mcGraph1(heaven, hell, crime)
##' )
##' with(religioncrime,
##' mcGraph2(heaven, hell, crime)
##' )
##' mod1 <- with(religioncrime,
##' mcGraph3(heaven, hell, crime)
##' )
##' summary(mod1[[1]])
##' ##TODO: Draw more with perspective matrix mod1[[2]]
NULL



##' Cheating and Looting in Japanese Electoral Politics
##'
##' Extracted from the "cheating-replication.dta" data file
##' with permission by the authors, Benjamin Nyblade and Steven
##' Reed. The Stata data file provided by the authors included
##' many constructed variables that have been omitted.  Within
##' R, these can be easily re-contructed by users.
##'
##' Special thanks to NyBlade and Reed for permission to repackage
##' this data. Also special thanks to them for creating an especially
##' transparent variable naming scheme.
##'
##' The data set includes many columns for variables that can easily
##' be re-constructed from the columns that are provided here.  While
##' Stata users might need to manually create 'dummy variables' and
##' interactions, R users generally do not do that manually.
##'
##' These variables from the original data set were omitted:
##'
##' Dummy variables for the year variable: c("yrd1", "yrd2", ...,
##' "yrd17", "yrd18")
##'
##' Dummy variables for the ku variable:
##' c("ku1", "ku2", ..., "ku141", "ku142")
##'
##' Constructed product variables: c("actualratiosq", "viabsq",
##' "viab_candcamp_divm", "viab_candothercamp_divm",
##' "viabsq_candcamp_divm", "viabsq_candothercamp_divm",
##' "absviab_candcamp", "absviab_candothercamp",
##' "absviab_candcamp_divm", "absviab_candothercamp_divm",
##' "viabsq_candcamp", "viabsq_candothercamp", "viab_candcamp",
##' "viab_candothercamp", "candothercamp_divm", "candcamp_divm",
##' "candcampminusm", "candothercampminusm", "predratiosq", "absviab")
##'
##' Mean centered variables: constr2 <- c("viab_candcampminusm",
##' "viab_candothercampminusm", "viabsq_candothercampminusm",
##' "viabsq_candcampminusm")
##'
##' In the end, we are left with these variables:
##'
##' [1] "ku"
##' [2] "prefecture"
##' [3] "dist"
##' [4] "year"
##' [5] "yr"
##' [6] "cdnr"
##' [7] "jiban"
##' [8] "cheating"
##' [9] "looting"
##' [10] "actualratio"
##' [11] "viab"
##' [12] "inc"
##' [13] "cons"
##' [14] "ur"
##' [15] "newcand"
##' [16] "jwins"
##' [17] "cons_cwins"
##' [18] "oth_cwins"
##' [19] "camp"
##' [20] "fleader"
##' [21] "incablast"
##' [22] "predratio"
##' [23] "m"
##' [24] "candcamp"
##' [25] "candothercamp"
##' [26] "kunocheat"
##' [27] "kunoloot"
##'
##' @name cheating
##' @docType data
##' @usage data(cheating)
##' @author Paul E. Johnson \email{pauljohn@@ku.edu}, on behalf of Benjamin Nyblade and Steven Reed
##' @format data.frame: 16623 obs. on 27 variables
##' @keywords datasets
##' @source \url{http://faculty.arts.ubc.ca/bnyblade/publications.html}.
##' @references
##' Benjamin Nyblade and Steven Reed, "Who Cheats? Who
##' Loots? Political Competition and Corruption in Japan, 1947-1993."
##' American Journal of Political Science 52(4): 926-41. October 2008.
##' @examples
##' require(rockchalk)
##' data(cheating)
##'
##' table1model2 <- glm(cheating ~ viab + I(viab^2) + inc + cons + ur
##' + newcand + jwins + cons_cwins + oth_cwins, family = binomial(link
##' = "logit"), data = cheating)
##'
##' predictOMatic(table1model2)
##'
##' predictOMatic(table1model2, interval = "confidence")
##'
##' ## The publication used "rare events logistic", which I'm not bothering
##' ## with here because I don't want to invoke additional imported packages.
##' ## But the ordinary logit results are proof of concept.
NULL
