#' The Stratified Analysis with Nonparametric covariable adjustment Package
#'
#' A Package for Implementation of the method in Kawaguchi, Koch, and Wang (2011)
#'
#' @name sanon-package
#' @aliases sanon-package
#' @rdname sanon-package
#' @docType package
#' @keywords documentation
#'
#' @author Atsushi Kawaguchi. \email{kawa_a24@@yahoo.co.jp}
#' @seealso \code{\link{sanon}}
#' @references 
#' Kawaguchi A., Koch, G. G. (2015). sanon: An R Package for Stratified Analysis with Nonparametric Covariable Adjustment. Journal of Statistical Software, 67(9), 1-37. doi:10.18637/jss.v067.i09
#' 
#' Kawaguchi A., Koch, G. G., Wang, X. (2011). Stratified Multivariate Mann-Whitney Estimators for the Comparison of Two Treatments with Randomization Based Covariance Adjustment. Statistics in Biopharmaceutical Research, Vol. 3, No. 2, 217-231. 
#' @importFrom stats coef get_all_vars model.extract model.frame na.action na.omit na.pass pchisq printCoefmat qnorm relevel terms vcov
NULL

#' Chronic Pain Data
#'
#' The data are from a multicenter randomized clinical trial to compare test and control treatments for the management of chronic pain, and they have had previous consideration in Stokes et al. (2000, chap. 13).
#' 
#'  \describe{
#'    \item{\code{treat}}{a factor with levels \code{active} and \code{placebo} for treatment}
#'    \item{\code{response}}{a factor with five levels \code{poor}, \code{fair}, \code{moderate}, \code{good} and \code{excel} for pain status after treatment for 4 weeks}
#'    \item{\code{center}}{a factor with two levels \code{I} and \code{II} for two centers}
#'    \item{\code{diagnosis}}{a factor with four levels \code{A}, \code{B}, \code{C}, and \code{D} for diagnoses}
#'  }
#'
#' @docType data
#' @keywords datasets
#' @name cpain
#' @usage data(cpain)
#' @format A data frame with 193 observations and 4 variables
#' @references Stokes, M. E., Davis, C. S., and Koch, G. G. (2000), Categorical Data Analysis using the SAS System, Cary: SAS Publishing.
NULL

#' Respiratory Disorder Data
#'
#' The data are from a randomized clinical trial to compare a test treatment to placebo for a respiratory disorder, and listings of the data appear in Stokes et al. (2000, chap. 15, pp. 495-496) and Koch et al. (1990).
#' The variables are as follows:
#'
#'  \describe{
#'    \item{\code{center}}{a factor vector for two centers}
#'    \item{\code{treatment}}{a factor with levels \code{A} and \code{P} for active and placebo treatments, respectively}
#'    \item{\code{sex}}{a factor with levels \code{F} and \code{M} for female and male, respectively}
#'    \item{\code{age}}{a numeric vector for age}
#'    \item{\code{baseline}}{a numeric vector for patient global ratings of symptom control according to 5 categories (4 = excellent, 3 = good, 2 = fair, 1 = poor, 0 = terrible) at baseline measurement}
#'    \item{\code{visit1}}{a numeric vector for patient global ratings of symptom control at visit 1 with same categories as \code{baseline}}
#'    \item{\code{visit2}}{a numeric vector for patient global ratings of symptom control at visit 2 with same categories as \code{baseline}}
#'    \item{\code{visit3}}{a numeric vector for patient global ratings of symptom control at visit 3 with same categories as \code{baseline}}
#'    \item{\code{visit4}}{a numeric vector for patient global ratings of symptom control at visit 4 with same categories as \code{baseline}}
#'  }
#'
#' @docType data
#' @keywords datasets
#' @name resp
#' @usage data(resp)
#' @format A data frame with 111 observations and 9 variables.
#' @references 
#' Stokes, M. E., Davis, C. S., and Koch, G. G. (2000), Categorical Data Analysis using the SAS System, Cary: SAS Publishing.
#' 
#' Koch, G. G., Carr, G. J., Amara, I. A., Stokes, M. E., and Uryniak, T. J. (1990), "Categorical Data Analysis," in Statistical Methodology in Pharmaceutical Sciences, ed. D. A. Berry, New York: Marcel Dekker, pp. 291-475.
NULL

#' Seborrheic Dermatitis Data
#'
#' The data are from a randomized clinical trial to compare a test treatment to placebo for a seborrheic dermatitis, and listings of the data appear in Ramaswamy, Koch, and Amara (1997).
#' The variables are as follows:
#'
#'  \describe{
#'    \item{\code{center}}{a factor vector for eight centers}
#'    \item{\code{treat}}{a factor with levels \code{placebo} and \code{test} for placebo and test treatments, resectively}
#'    \item{\code{score1}}{a numeric vector for patient global scores for the face according to 6 categories (0 = cleared, 1 = excellent improvement, 2 = moderate improvement, 3 = slight improvement, 4 = no change, 5 = exacerbation)}
#'    \item{\code{score2}}{a numeric vector for patient global scores for the scalp with same categories as \code{score1}}
#'    \item{\code{score3}}{a numeric vector for patient global scores for the chest with same categories as \code{score1}}
#'    \item{\code{severity1}}{a numeric vector for the baseline desease severity for the face according to 3 categories (1 = mild, 2 = moderate, 3 = severe)}
#'    \item{\code{severity2}}{a numeric vector for the baseline desease severity for the scalp with same categories as \code{severity1}}
#'    \item{\code{severity3}}{a numeric vector for the baseline desease severity for the chest with same categories as \code{severity1}}
#'  }
#'
#' @docType data
#' @keywords datasets
#' @name sebor
#' @usage data(sebor)
#' @format A data frame with 167 observations and 8 variables.
#' @references 
#' Ramaswamy R, Koch G, Amara I (1997). "Application of rank analysis of covariance methods to analysis of multiple anatomical regions with treatment for seborrheic dermatitis." Journal of Biopharmaceutical Statistics, 7(3), 403--416.
NULL

#' Skin Condition Data
#'
#' The data are from a randomized clinical trial to compare a test treatment to placebo for skin conditions, and listings of the data appear in Stanish, Gillings, Koch (1978a, b).
#' The variables are as follows:
#'
#'  \describe{
#'    \item{\code{center}}{a factor vector for two centers}
#'    \item{\code{treat}}{a factor with levels \code{A} and \code{P} for active and placebo treatments, skinectively}
#'    \item{\code{stage}}{a numeric vector for initial severity of the skin condition according to 3 categories (3 = fair, 4 = poor, 5 = exacerbation) at baseline measurement}
#'    \item{\code{res1}}{a numeric vector for extent of improvement at visit 1 according to 5 categories (1 = rapidly improving, 2 = slowly improving, 3 = stable, 4 = slowly worsening, 5 = rapidly worsening)}
#'    \item{\code{res2}}{a numeric vector for extent of improvement at visit 2 with same categories as \code{res1}}
#'    \item{\code{res3}}{a numeric vector for extent of improvement at visit 3 with same categories as \code{res1}}
#'  }
#'
#' @docType data
#' @keywords datasets
#' @name skin
#' @usage data(skin)
#' @format A data frame with 172 observations and 6 variables.
#' @references 
#' Stanish W, Gillings D, Koch G (1978a). "An application of multivariate ratio methods for the analysis of a longitudinal clinical trial with missing data." Biometrics, 34(2), pp. 305--317.
#'
#' Stanish WM, Koch GG, Landis JR (1978b). "A computer program for multivariate ratio analysis (MISCAT)." Computer Programs in Biomedicine, 8(3-4), 197--207.
NULL

#' Relief of heartburn Data
#'
#' The data are from two period cross-over design clinical trial for relief of heartburn, and listings of the data appear in Koch, Gitomer, Skalland, and Stokes (1983).
#' The variables are as follows:
#'
#'  \describe{
#'    \item{\code{center}}{a factor vector for two centers}
#'    \item{\code{sequence}}{a factor with levels \code{AP} and \code{PA} for sequence groups}
#'    \item{\code{age}}{a numeric vector for age}
#'    \item{\code{sex}}{a factor for sex with levels \code{female} and \code{male}}
#'    \item{\code{freq}}{a numeric vector for weekly frequency of condition from previous medical history}
#'    \item{\code{MD1}}{a numeric vector for time to relief from first dose during period 1}
#'    \item{\code{MD2}}{a numeric vector for time to relief from first dose during period 2}
#'    \item{\code{res1}}{a factor vector for relief status for period 1 (R = relief from first dose within 15 min, NF = no relief from first dose within 15 min)}
#'    \item{\code{ref2}}{a factor vector for relief status for period 2 with same categories as \code{res1}}
#'  }
#'
#' @docType data
#' @keywords datasets
#' @name heartburn
#' @usage data(heartburn)
#' @format A data frame with 60 observations and 9 variables.
#' @references 
#' Koch G, Gitomer S, Skalland L, Stokes M (1983). "Some non-parametric and categorical data analyses for a change-over design study and discussion of apparent carry-over effects." Statistics in Medicine, 2(3), 397--412.
NULL

#' Non-Parametric Covariable Adjustment for Stratified Rank Measures of Association
#'
#' This is a function for computing a stratified multivariate Mann-Whitney estimator that addresses the comparison between two randomized groups for a strictly ordinal response variable.
#' Response variables may have some missing completely at random (MCAR) values for some patients. 
#' Non-parametric covariable adjustment is considered through the difference estimates between mean covariable and the weighted least squares method.
#' Although such estimators can be computed directly as weighted linear combinations of within-stratum Mann-Whitney estimators, consistent estimation of their covariance matrix is done using methods for multivariate U-statistics.
#' 
#' \code{sanon} has two specifications for the input, variable and formula based.
#' In the variable based input, one can specify R objects to outcome, group, and strata variables, and covariable.
#' In the formula based input, the formula consists of variable names in a data.frame. 
#' The strata and group variables, and covariable are recognized by functions \code{\link{strt}}, \code{\link{grp}}, \code{\link{covar}}, and \code{\link{catecovar}}.
#' \code{outcome} can be contained missing values, which should be coded by \code{NA}. 
#' Five options for the management of missing values can be specifed in the argument \code{res.na.action}; 
#' \code{"default"} = the method in Kawaguchi et al. (2011), \code{"LOCF1"} and \code{"LOCF2"} = last observation carried forward with respect to kernels of U-statistics and observed velues, repsectively, \code{"replace"} = missing values are managed as tied with all other values in the same stratum, and \code{"remove"} = the complete cases analaysis.
#' For \code{res.na.action = "LOCF1"} or \code{"LOCF2"}, the order in the outcome is considered as the time order in imputing.
#' if the baseline measurement is missing, then the corresponding subject is removed.
#' \code{outcome} can be also multiple (repeatly measured).
#' If more than two strata are specified, these are taking a cross-classification. 
#' The group variable can be specifies its reference group in the argument \code{ref} in the \code{sanon} or in the function \code{grp}.
#'
#' @name sanon
#' @aliases sanon
#' @rdname sanon
#' @docType methods
#' @export
#'
#' @param formula a formula object, with the response on the left of a ~ operator, and the terms on the right.
#' @param data a data.frame in which to interpret the variables named in the formula. 
#' @param outcome vector of observations of length n, or a matrix with n rows for the response (or outcome) variables
#' @param group numeric vector of observations of length n for treatment group. The reference group can be specified in \code{ref}.
#' @param strt numeric or factor vector of observations of length n, or a matrix with n rows for strata.
#' @param covar numeric or factor vector of observations of length n, or a matrix with n rows for covariable.
#' @param catecovar numeric or factor vector of observations of length n, or a matrix with n rows for categorical covariable.
#' @param ref character for the reference group for treatment group in \code{group}.
#' @param covref character vector for the reference group for categorical covariables in \code{catecovar}.
#' @param P a matrix for weighted least squares estimation.
#' @param res.na.action character for setting NA actions. "default", "LOCF1", "LOCF2", "replace", and "remove" are available. default is "default". see the details.
#' @param x an object of class "\code{sanon}", usually, a result of a call to \code{\link{sanon}}
#' @param ... further arguments passed to or from other methods.
#' @return \item{N}{Sample size}
#' @return \item{Nna}{tne number of subjects with missing values}
#' @return \item{nhik}{Sample size in each strata, group, and response}
#' @return \item{nik}{Sample size in each group and response}
#' @return \item{xi}{(multivariate) Mann-Whitney estimate(s) that addresses the comparison between two randomized groups}
#' @return \item{g}{the difference estimates between mean covariable}
#' @return \item{f}{a vector consisting of \code{xi} and \code{g}}
#' @return \item{Vf}{estimated covariance matrix of \code{f}}
#' @return \item{b}{fully adjustmented estimators for all covariables and the strata}
#' @return \item{Vb}{covariance matrix of \code{b}}
#' @return \item{se}{standard error of \code{b}}
#' @return \item{Q}{test statistics for \code{b}}
#' @return \item{p}{p-value for \code{b}}
#' @return \item{outnames}{outcome or response names}
#' @return \item{covarnames}{covariable names}
#' @return \item{advarnames}{variable names adjusting in the weighted least squares}
#' @return \item{bnames}{variable names of adjusted in the weighted least squares}
#' @return \item{reslevels}{levels for response variables}
#' @return \item{grouplevels}{levels for the group variable}
#' @return \item{strtout}{resulting (cross-classification) strata}
#' @return \item{strtlevels}{resulting (cross-classification) strata levels}
#' @return \item{strtnames}{resulting (cross-classification) strata names}
#' @return \item{matP}{design matrix used in the weighted least squares}
#'
#' @references 
#' Kawaguchi A., Koch, G. G. (2015). sanon: An R Package for Stratified Analysis with Nonparametric Covariable Adjustment. Journal of Statistical Software, 67(9), 1-37. doi:10.18637/jss.v067.i09
#' 
#' Kawaguchi, A., Koch, G. G., Wang, X. (2011): Stratified Multivariate Mann-Whitney Estimators for the Comparison of Two Treatments with Randomization Based Covariance Adjustment. Statistics in Biopharmaceutical Research, Vol. 3, No. 2, 217-231. 
#' @examples
#' ##### Example 3.1 Randomized Clinical Trial of Chronic Pain #####
#' data(cpain)
#' out11 = sanon(response ~ grp(treat, ref="placebo") + strt(center) + strt(diagnosis), data=cpain)
#' out11
#' summary(out11)
#'
#' # R objects are also available
#' attach(cpain)
#' out12 = sanon(outcome=response, group=treat, 
#' strt=cbind(center, diagnosis), ref="placebo")
#' out12
#' summary(out12)
#'
#' ##### Example 3.2 Randomized Clinical Trial of Respiratory Disorder #####
#' data(resp)
#' out21 = sanon(cbind(baseline, visit1, visit2, visit3, visit4) 
#' ~ grp(treatment, ref="P") + strt(center) + strt(sex) + covar(age), data=resp)
#' out21
#' summary(out21)
#'
#' # the matrix P can be specified
#' P = rbind(rep(0, 4), diag(4), rep(0, 4))
#' out22 = sanon(cbind(baseline, visit1, visit2, visit3, visit4) 
#' ~ grp(treatment, ref="P") + strt(center) + strt(sex) + covar(age), data=resp, P=P)
#' out22
#' summary(out22)
#'

sanon = function(outcome, ...) UseMethod("sanon")

#' @rdname sanon
#' @export
sanon.formula = function(formula, data=list(), ...)
{
mf = model.frame(formula=formula, data=data, na.action = na.pass)
outcome = model.extract(mf, "response")
#if(is.null(ncol(outcome))){outcome = data.frame(outcome); colnames(outcome) = names(mf)[1]}
if(is.null(ncol(outcome))){outnames = names(mf)[1]}else{outnames = colnames(model.extract(mf, "response"))}
outcome = get_all_vars(formula=formula, data=data)[outnames]

special = c("catecovar", "covar", "strt", "grp")
formula2 = terms(formula, special, data=data)
mf2 = model.frame(formula2, data=data, na.action = na.pass)

x = lapply(attr(formula2, "specials"), function(x){if(is.null(x)){x}else{mf2[x]}})

group = x$grp; if(!is.null(group)) names(group) = strsplit(substr(names(group), 5, nchar(names(group))-1), ",")[[1]][1]
ref = levels(group[,1])[1]
strt=x$strt; if(!is.null(strt)) names(strt) = substr(names(strt), 6, nchar(names(strt))-1)
covar=x$covar; if(!is.null(covar)) names(covar) = substr(names(covar), 7, nchar(names(covar))-1)
catecovar=x$catecovar; if(!is.null(catecovar)) names(catecovar) = unlist(lapply(strsplit(substr(names(catecovar), 11, nchar(names(catecovar))-1), ","), function(z) z[1]))
if(!is.null(ncol(catecovar))) covref = unlist(sapply(1:ncol(catecovar), function(x) levels(catecovar[,x])[1]))

est = sanon.default(outcome=outcome, group=group, strt=strt, covar=covar, catecovar = catecovar, ref=ref, covref = covref,...)
est$call = match.call()
est$formula = formula
est
}

#' @rdname sanon
#' @method sanon default
#' @export
sanon.default = function(outcome, group, strt=NULL, covar=NULL, catecovar=NULL, ref=NULL, covref=NULL, P=NULL, res.na.action = "default", ...)
{
#######################
##### Requirement #####
#######################
if(missing(outcome)) stop("outcome should be specified")
if(missing(group)) stop("group should be specified")
if(class(P) != "matrix" & !is.null(P)) stop("P should be the matrix class")
if(!(res.na.action %in% c("default", "LOCF1", "LOCF2", "replace", "remove"))) stop(paste(res.na.action, "not a option for res.na.action"))

############################
##### Set up variables #####
############################

##### response variable #####
Y = as.data.frame(outcome); outnames = colnames(Y)
if(res.na.action == "remove") Y = na.omit(Y)
naY = na.action(Y)
if(!is.null(ncol(Y)))
{
reslevels = lapply(1:ncol(Y), function(i) levels(as.factor(Y[,i])))
for(i in 1:ncol(Y)){if(!(class(Y[,i]) %in% c("numeric", "integer"))) Y[,i] = as.numeric(Y[,i])}
}else
{
Y = as.matrix(Y, 1)
reslevels = levels(as.factor(Y))
if(class(Y) != "numeric") Y = as.numeric(Y)
}

names(reslevels) = outnames

### Last observation carried forward (LOCF) for observed values ###
if(res.na.action == "LOCF2"){
if(ncol(Y) == 1) stop("LOCF2 can not be applied to a single response")
# naY = na.action(na.omit(Y[,1])); Y = Y[-naY,]
Y = t(apply(Y, 1, function(x){ for(i in 2:length(x)){ x[i] = ifelse(is.na(x[i]), x[i-1], x[i])};x }))
}

## check ##
#resuniq = apply(Y, 2, function(x) length(unique(x)) == 1)
#if(any(resuniq)) stop("Responses should not have a unique value")

##### Numbers of subjects and responses #####
N = ifelse(!is.null(ncol(Y)), nrow(Y), length(Y))
Nna = length(naY)
r = ncol(Y)

##### Strata variable #####
if(is.null(strt)){strtnames = NULL; strt = rep(1, N+Nna)}else{strtnames = paste(colnames(strt), collapse="*")}
S = as.data.frame(strt)
if(!is.null(na.action(na.omit(S)))) stop("strata should not have missing values (NA)")
if(!is.null(ncol(S))) S = apply(S, 1, function(x) paste(x, collapse="*"))

if(!is.null(naY)) S = S[-naY]

strtout = as.factor(S) # for output
strtlevels = levels(strtout)
S = as.numeric(as.factor(S))

##### group variable #####
if(class(group) == "data.frame"){ t = group[,1]}else{t = group}
if(!is.null(na.action(na.omit(t)))) stop("group should not have missing values (NA)")
if(!is.null(ncol(t))){if(ncol(t) > 1) stop("duplicated group variable")}
if(length(unique(t)) != 2) stop("group should have two categories")

if(is.null(ref)){grouplevels = levels(factor(t))}else{grouplevels = levels(relevel(factor(t), ref=ref))}
t = ifelse(t == grouplevels[1], -1, 1)
if(class(t) != "numeric") t = as.numeric(t)
if(!is.null(naY)) t = t[-naY]

## check ##
#resuniq2 = apply(Y, 2, function(x) aggregate(x, by=list(t), function(y) length(unique(y)) == 1)[,2])
#if(any(resuniq2)) stop("Responses should not have a unique value in each group")

##### covariable #####
X = NULL
if(!is.null(covar)){ X = as.data.frame(covar)}
## dummy variables for categorical cavariables ##
if(!is.null(catecovar)){
catecovar = as.data.frame(catecovar)
catecovar = do.call(cbind, lapply(1:ncol(catecovar), function(x){
tmp = catecovar[,x]
tmpcate = unique(tmp)
tmpref = tmpcate[length(tmpcate)]
if(!is.null(covref[x])){
if(!any(tmpcate %in% covref[x])) stop(paste(covref[x], "not in categories"))
tmpref = covref[x]
}
tmpcate2 = tmpcate[!(tmpcate %in% tmpref)]
out = sapply(tmpcate2, function(y) as.numeric(tmp == y))
colnames(out) = paste(colnames(catecovar)[x], "[", tmpcate2, "/", tmpref, "]", sep="")
out
}))
X = cbind(X, catecovar)
}

## ##
if(!is.null(X)){
if(!is.null(na.action(na.omit(X)))) stop("covariables should not have missing values (NA)")
if(!is.null(ncol(X)))
{
for(i in 1:ncol(X)){if(!(class(X[,i]) %in% c("numeric", "integer"))) X[,i] = as.numeric(X[,i])}
}else
{
X = as.matrix(X, 1)
if(class(X) != "numeric") X = as.numeric(X)
}

covarnames = colnames(X)
if(!is.null(naY)) X = X[-naY,]
if(class(X) != "matrix") X = as.matrix(X)
}else{covarnames = NULL}

## Number of covariables ##
M = ifelse(!is.null(X), ncol(X), 0)

##### Missing Indicator #####
if(res.na.action == "LOCF1"){
if(ncol(Y) == 1) stop("LOCF1 can not be applied to a single response")
Z = matrix(1, nrow(Y), ncol(Y))
Z[,1] = as.numeric(!is.na(Y[,1]))
Y[,1] = ifelse(is.na(Y[,1]), 0, Y[,1])
}else{
Z = apply(Y, 2, function(x) as.numeric(!is.na(x)))
Y = apply(Y, 2, function(x) ifelse(is.na(x), 0, x))
}
if(class(Y) != "matrix") Y = as.matrix(Y)

##### sample size within treatment group and strata #####
n0 = lapply(1:ncol(Z), function(x) table(factor(t)[Z[,x] == 1], factor(S)[Z[,x] == 1]))
n = sapply(1:r, function(i) sapply(1:N, function(x) n0[[i]][as.character(t[x]), as.character(S[x])]))
n20 = sapply(1:N, function(x) table(t, S)[as.character(t[x]), as.character(S[x])])
n2 = n20 %o% rep(1, M)
n2r = n20 %o% rep(1, r)

# for output
names(n0) = outnames
for(i in 1:length(n0)) dimnames(n0[[i]]) = list(grouplevels, strtlevels)

# for check
n1 = lapply(n0, function(x) rowSums(x))
njudge1 = do.call(cbind, n1) < 50
njudge2 = lapply(n0, function(x) x < 4)

if( any(njudge1) )
{
njudge1wh = which(njudge1, arr.ind = TRUE)
njudge1names = apply(njudge1wh, 1, function(x) c(rownames(njudge1)[x[1]], colnames(njudge1)[x[2]]))
#warning("Sample size is not greater than 50 in:\n", paste(apply(njudge1names, 2, function(x) paste(x[1], "group for", x[2])), collapse=", "))
}

if( any(unlist(njudge2)) )
{
njudge2wh = lapply(njudge2, function(x) which(x, arr.ind = TRUE))
njudge2names = lapply(njudge2wh, function(y) apply(y, 1, function(x) c(rownames(njudge2[[1]])[x[1]], colnames(njudge2[[1]])[x[2]])))
njudge2names = njudge2names[unlist(lapply(njudge2names, class)) == "matrix"]
njudge2names2 = lapply(njudge2names, function(y) paste(apply(y, 2, function(x) paste(x[1], "group in strata", x[2])), collapse=","))
#warning("Sample size is not greater than 4 in:\n", paste(sapply(1:length(njudge2names2), function(x) paste(njudge2names2[[x]], "for", names(njudge2names2)[x])), collapse=", "))
}

###############################################
#####========== Computing Start ==========#####
###############################################

##### Function for computing kernels of U-statistics #####
U = function(j1)
{
j2 = (1:N)[-j1]

## ##
tmpS = ifelse(S[j1] - S[j2] == 0, 1, 0)
tmpY = ifelse((Y[rep(j1,N-1),] - Y[j2,]) == 0, 1, 0)
tmptY = ifelse((t[j1] - t[j2]) * (Y[rep(j1,N-1),] - Y[j2,]) * Z[rep(j1,N-1),] * Z[j2,] > 0, 1, 0)
tmpt = ifelse((t[j1] - t[j2])^2 * Z[rep(j1,N-1),] * Z[j2,]  > 0, 1, 0)
if(res.na.action %in% c("replace", "LOCF1", "LOCF2")) tmpt2 = ifelse((t[j1] - t[j2])^2 * (1 - Z[rep(j1,N-1),] * Z[j2,]) > 0, 1, 0)

## ##
if(res.na.action %in% c("replace", "LOCF1", "LOCF2")){
tmpU1 = (tmpS * (tmptY + 0.5 * tmpt * tmpY + 0.5 * tmpt2)) / (n2r[rep(j1,N-1),] + n2r[j2,] + 1)
tmpU2 = (tmpS * (tmpt + tmpt2)) / (n2r[rep(j1,N-1),] + n2r[j2,] + 1)
}else
{
tmpU1 = (tmpS * (tmptY + 0.5 * tmpt * tmpY)) / (n[rep(j1,N-1),] + n[j2,] + 1)
tmpU2 = (tmpS * tmpt) / (n[rep(j1,N-1),] + n[j2,] + 1)
}

## carrying forward ##
dim(tmpU1) = dim(tmpU2) = c(N-1, r)
if(res.na.action == "LOCF1"){
misidx1 = apply(tmpU1, 1, function(x) any(is.na(x)))
misidx2 = apply(tmpU2, 1, function(x) any(is.na(x)))
tmpU1[misidx1,] = t(apply(tmpU1[misidx1,], 1, function(x){ for(i in 2:length(x)){ x[i] = ifelse(is.na(x[i]), x[i-1], x[i])};x }))
tmpU2[misidx2,] = t(apply(tmpU2[misidx2,], 1, function(x){ for(i in 2:length(x)){ x[i] = ifelse(is.na(x[i]), x[i-1], x[i])};x }))
}

## ##
c(apply(tmpU1, 2, mean), apply(tmpU2, 2, mean))
}

##### Implementation of function #####
tmpU = do.call(rbind, lapply(1:N, function(x) U(x)))

U1j = tmpU[,1:r]
U2j = tmpU[,r+(1:r)]

thetas = apply(tmpU, 2, mean)
theta1 = thetas[1:r]
theta2 = thetas[r+(1:r)]

##### #####
F = cbind(U1j, U2j)
meanF = apply(F, 2, mean)

##### #####
tmpF = sweep(F, 2, meanF)
VF = 4 * t(tmpF) %*% (tmpF) / (N*(N-1))

##### #####
if(r > 1){
Dtheta1 = diag(c(theta1))
Dtheta2 = diag(c(theta2))
xi = solve(Dtheta2) %*% theta1
}else
{
xi = theta1 / theta2
}

##### #####
if(r > 1){
Dxi = diag(c(xi))
#Dthetas = cbind(solve(Dtheta1), -solve(Dtheta2))
#Vxi = Dxi %*% Dthetas %*% VF %*% t(Dthetas) %*% Dxi
Dthetas = cbind(diag(r), -solve(Dtheta2)^2)
Vxi = Dthetas %*% VF %*% t(Dthetas)
}else
{
Dthetas = cbind(1/theta1, -1/theta2)
Vxi = xi * Dthetas %*% VF %*% t(Dthetas) * xi
}

##### the estimator from two-way analysis of variance for the difference between the stratification adjusted means of the mth covariable for the two groups. #####
if(!is.null(X)){
tU = function(j1)
{
j2 = (1:N)[-j1]

tmpS = ifelse(S[j1] - S[j2] == 0, 1, 0) 
tmpt1 = 0.5 * (t[j1] - t[j2])
tmpt2 = ifelse(t[j1] - t[j2] != 0, 1, 0)
tmpX = X[rep(j1,N-1),] - X[j2,]

tmpU1 = (tmpS * tmpt1 * tmpX) / (n2[rep(j1,N-1),] + n2[j2,])
tmpU2 = (tmpS * tmpt2) / (n2[rep(j1,N-1),] + n2[j2,])

dim(tmpU1) = dim(tmpU2) = c(N-1, M)
c(apply(tmpU1, 2, mean), apply(tmpU2, 2, mean))
}

tmptU = do.call(rbind, lapply(1:N, function(x) tU(x)))

tU1j = tmptU[,1:M]
tU2j = tmptU[,M+(1:M)]

phis = apply(tmptU, 2, mean)
phi1 = phis[1:M]
phi2 = phis[M+(1:M)]

g = matrix(phi1 / phi2, ncol=1)
#if(any(g==0)) stop(paste(paste(covarnames[which(g==0)], collapse=", "), "completely balanced"))

}else{
g = NULL
}

##### #####
if(!is.null(X)){
G = cbind(U1j, tU1j, U2j, tU2j)
}else{
G = cbind(U1j, U2j)
}
meanG = apply(G, 2, mean)

##### #####
tmpG = sweep(G, 2, meanG)
VG = 4 * t(tmpG) %*% (tmpG) / (N*(N-1))

##### #####
if(!is.null(X)){
f = rbind(xi, g)
tG = c(theta1, phi1, theta2, phi2)
}else{
f = xi
tG = c(theta1, theta2)
}

##### Variance Covariance Matrix for f #####
if(length(f) > 1){
#H0 = rbind(cbind(diag(r), matrix(0, r, M), -diag(r), matrix(0, r, M)), cbind(matrix(0, M, r), diag(M), matrix(0, M, r), -diag(M))); #H = diag(c(f)) %*% H0 %*% diag(1/tG)
if(r > 1){tmptheta = diag(theta2^(-2))%*%diag(theta1); invDtheta2 = solve(Dtheta2)}else{tmptheta = theta1/theta2^2; invDtheta2 = 1/theta2}
if(M > 1){tmpphi = diag(phi2^(-2))%*%diag(phi1); invphi2 = diag(phi2^(-1))}else if(M > 0){tmpphi = phi1/phi2^2; invphi2 = 1/phi2}else{tmpphi = invphi2 = matrix(0, M, M)}
H = rbind(cbind(invDtheta2, matrix(0, r, M), -tmptheta, matrix(0, r, M)), cbind(matrix(0, M, r), invphi2, matrix(0, M, r), -tmpphi))
}else{
H0 = cbind(1, -1)
H = f %*% H0 %*% diag(1/tG)
}
Vf = H %*% VG %*% t(H)

##### Covariable Adjusted Estimators #####
if(!is.null(X)){
if(is.null(P)) P = rbind(diag(r), matrix(0, M, r))
f = rbind(xi-0.5, g)
}else{
if(is.null(P)) P = diag(r)
f = xi-0.5
}

allnames = c(outnames, covarnames)
if(nrow(P) != length(allnames)){stop("The number of row of P should be r+M")}else{rownames(P) = allnames}
advarnames = allnames[apply(P, 1, function(x) all(x == 0))]
bnames = apply(P, 2, function(x) paste(rownames(P)[which(x == 1)], collapse=" + "))#allnames[!(allnames %in% advarnames)]

invVf = try(solve(Vf), silent = TRUE)
if(class(invVf) == "try-error"){
warning("Vf is computationally singular.")
e = 0.0000001
while(class(invVf) == "try-error"){
invVf = try(solve(Vf + e*diag(ncol(Vf))), silent = TRUE)
e = 2*e
}}

b = solve(t(P) %*% invVf %*% P) %*% (t(P) %*% invVf %*% f)
Vb = solve(t(P) %*% invVf %*% P)

if(r > 1){se = sqrt(diag(Vb))}else{se = sqrt(Vb)}

##### Test statistic and p-value #####
Q = (b / se)^2
pQ = 1 - pchisq(Q, 1)

##### Output #####
out = list(N=N, Nna=Nna, nhik=n0, nik=n1,
xi=xi, g=g, f=f, Vf=Vf, b=b, Vb=Vb, se=se, Q=Q, p=pQ, call = match.call(), 
outnames = outnames, covarnames = covarnames, advarnames = advarnames, bnames = bnames,
reslevels = reslevels, grouplevels = grouplevels,
strtout = strtout, strtlevels = strtlevels, strtnames = strtnames,
matP = P)
class(out) = "sanon"
out
}

#' @rdname sanon
#' @method print sanon
#' @export
print.sanon = function(x, ...)
{
tmpxi = c(x$xi); names(tmpxi) = x$outnames
tmpb = c(x$b); names(tmpb) = x$bnames; tmpb[!(x$bnames %in% x$covarnames)] = tmpb[!(x$bnames %in% x$covarnames)] + 0.5
cat("Call:\n")
print(x$call)
cat("\n")
cat(paste("Sample size:", x$N))
if(x$Nna > 0) cat(paste(" (", x$Nna, "samples removed)"))
cat("\n\n")
if(!is.null(x$strtnames)) cat(paste("Strata (", x$strtnames, "):", paste(x$strtlevels, collapse=", "), "\n\n"))
cat("Response levels:\n")
reslevels = lapply(x$reslevels, function(z){if(length(z)>7){c(z[1:3], "...", z[length(z)+(-2:0)])}else{z}})
cat( paste(sapply(1:length(reslevels), function(z) paste("[", names(reslevels)[z], "; ", length(x$reslevels[[z]]), " levels] (lower) ", paste(reslevels[[z]], collapse=", "), " (higher)", sep="")), collapse="\n"))
cat("\n\n")
cat("Design Matrix:\n")
print(x$matP)
cat("\n")
if(!is.null(x$strtnames)) cat("Stratification Adjusted ")
cat(paste("Mann-Whitney Estimate \n for comparison [",  x$grouplevels[2], "/", x$grouplevels[1], "] :\n"))
if(length(x$advarnames) > 0){ 
tmpadvarnames = x$advarnames
lim = max(which(cumsum(nchar(tmpadvarnames)) < 50))
if(lim < length(tmpadvarnames)){cat(paste("(adjusted by",  paste(c(tmpadvarnames[1:lim], ""), collapse=", "), "\n", paste(tmpadvarnames[(lim+1):length(tmpadvarnames)], collapse=", "), ")\n")) }else
{cat(paste("(adjusted by",  paste(tmpadvarnames, collapse=", "), ")\n"))}
}
print(round(tmpb, 4))
cat("\n")
}


#' Summarizing Weighted Least Squares Fits
#'
#' summary method for class "sanon". 
#'
#' This function provide the p value for the hypothesis test of coefficient in the model of weighted least squares method.
#' Note that the estimates in the output are for the (xi_k - 0.5).
#'
#' @name summary.sanon
#' @aliases summary.sanon
#' @rdname summary.sanon
#' @method summary sanon
#' @docType methods
#' @export
#'
#' @param object,x an object of class "\code{sanon}", usually, a result of a call to \code{\link{sanon}}
#' @param ... further arguments passed to or from other methods.
#' @return \item{coefficients}{a p x 4 matrix with columns for the estimated coefficient, its standard error, chi-squared statistic and corresponding (two-sided) p-value.}
#' @return \item{advarnames}{adjust variable names in weighted least squares method}
#'
#' @examples
#' ##### Example 3.1 Randomized Clinical Trial of Chronic Pain #####
#' data(cpain)
#' sum1 = summary(sanon(response ~ grp(treat, ref="placebo") + strt(center) + strt(diagnosis)
#' , data=cpain))
#' sum1
#'
#' ##### Example 3.2 Randomized Clinical Trial of Respiratory Disorder #####
#' data(resp)
#' sum22 = summary(sanon(cbind(baseline, visit1, visit2, visit3, visit4) 
#' ~ grp(treatment, ref="P") + strt(center) + strt(sex) + covar(age), data=resp))
#' sum22
#'

summary.sanon = function(object, ...)
{
tmpest = c(object$b); #tmpest[!(object$bnames %in% covarnames)] = tmpest[!(object$bnames %in% covarnames)] + 0.5
TAB = cbind(tmpest, object$se, object$Q, object$p)
colnames(TAB) = c("Estimate", "Std.Err", "Chisq", "Pr(>Chisq)")
rownames(TAB) = object$bnames
res = list(call=object$call, coefficients=TAB, matP = object$matP, advarnames = object$advarnames)
class(res) = "summary.sanon"
res
}

#' @rdname summary.sanon
#' @method print summary.sanon
#' @family print
#' @export
print.summary.sanon = function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\n")
if(length(x$advarnames) > 0)
{ 
tmpadvarnames = x$advarnames
lim = max(which(cumsum(nchar(tmpadvarnames)) < 50))
if(lim < length(tmpadvarnames)){cat(paste("Randomization-Based Covariance Adjusted Analysis \n (adjusted by",  paste(c(tmpadvarnames[1:lim], ""), collapse=", "), "\n", paste(tmpadvarnames[(lim+1):length(tmpadvarnames)], collapse=", "), "):"))}else
{cat(paste("Randomization-Based Covariance Adjusted Analysis \n (adjusted by",  paste(tmpadvarnames, collapse=", "), "):"))}


}
cat("\n")
printCoefmat(x$coefficients, digits=3, dig.tst = 2, P.values = TRUE, has.Pvalue = TRUE)
cat("Note that the estimates of responses are for the (MW estimate - 0.5).\n")
}

#' Identify Group Variables
#'
#' This is a special function used in the context of \code{\link{sanon}}. It identifies group variables when they appear on the right hand side of a formula. 
#'
#' @name grp
#' @aliases grp
#' @rdname grp
#' @docType methods
#' @export
#'
#' @param x variable name
#' @param ref character for the reference group for treatment group.
#'

grp = function(x, ref=NULL){
if(!is.null(ref)){
x = as.factor(x)
x = relevel(x, ref=ref)
}
x
}
#' Identify Stratification Variables
#'
#' This is a special function used in the context of \code{\link{sanon}}. It identifies stratification variables when they appear on the right hand side of a formula. 
#'
#' @name strt
#' @aliases strt
#' @rdname strt
#' @docType methods
#' @export
#'
#' @param x variable name
#'
#' @usage strt(x)
#'

strt = function(x) x

#' Identify Covariables
#'
#' This is a special function used in the context of \code{\link{sanon}}. It identifies covariables when they appear on the right hand side of a formula. 
#'
#' @name covar
#' @aliases covar
#' @rdname covar
#' @docType methods
#' @export
#'
#' @param x variable name
#'
#' @usage covar(x)
#'

covar = function(x) x

#' Identify Categorical Covariables
#'
#' This is a special function used in the context of \code{\link{sanon}}. It identifies categorical covariables when they appear on the right hand side of a formula. 
#'
#' In the \code{sanon}, the categorical covariable is converted into a dummy variable. The reference group is specified in the \code{ref} argument.
#'
#' @name catecovar
#' @aliases catecovar
#' @rdname catecovar
#' @docType methods
#' @export
#'
#' @param x variable name
#' @param ref character for the reference group for the categorical covariable.
#'

catecovar = function(x, ref=NULL){
if(!is.null(ref)){
x = as.factor(x)
x = relevel(x, ref=ref)
} 
x
}

#' Extract Model Coefficients
#' 
#' coef is a generic function which extracts model coefficients from objects returned by modeling functions. coefficients is an alias for it. 
#' 
#' All object classes which are returned by model fitting functions should provide a coef method or use the default one. 
#' 
#' @name coef.sanon
#' @aliases coef.sanon
#' @rdname coef.sanon
#' @method coef sanon
#' @docType methods
#' @export
#'
#' @param object an object of class "\code{sanon}", usually, a result of a call to \code{\link{sanon}}
#' @param ... further arguments passed to or from other methods.
#' @return Coefficients extracted from the model object object. 
#'
#' @examples
#' ##### Example 3.1 Randomized Clinical Trial of Chronic Pain #####
#' data(cpain)
#' out1 = sanon(response ~ grp(treat, ref="placebo") + strt(center) + strt(diagnosis), data=cpain)
#' coef(out1)
#' coefficients(out1)
#' 
#' ##### Example 3.2 Randomized Clinical Trial of Respiratory Disorder #####
#' data(resp)
#' P = rbind(rep(0, 4), diag(4), rep(0, 4))
#' out23 = sanon(cbind(baseline, visit1, visit2, visit3, visit4) ~ grp(treatment, ref="P")
#'  + strt(center) + strt(sex) + covar(age), data=resp, P=P)
#' # each four visits
#' coef(out23)
#' coefficients(out23)
#'

coef.sanon = function(object, ...){tmpb = c(object$b); names(tmpb) = object$bnames; tmpb}

#' Calculate Variance-Covariance Matrix for a Fitted Model Object
#' 
#' Returns the variance-covariance matrix of the main parameters of a fitted model object. 
#' 
#' This is a generic function. 
#' 
#' @name vcov.sanon
#' @aliases vcov.sanon
#' @rdname vcov.sanon
#' @method vcov sanon
#' @docType methods
#' @export
#'
#' @param object an object of class "\code{sanon}", usually, a result of a call to \code{\link{sanon}}
#' @param ... further arguments passed to or from other methods.
#' @return Coefficients extracted from the model object object. 
#'
#' @examples
#' ##### Example 3.1 Randomized Clinical Trial of Chronic Pain #####
#' data(cpain)
#' out1 = sanon(response ~ grp(treat, ref="placebo") + strt(center) + strt(diagnosis), data=cpain)
#' vcov(out1)
#' 
#' ##### Example 3.2 Randomized Clinical Trial of Respiratory Disorder #####
#' data(resp)
#' P = rbind(rep(0, 4), diag(4), rep(0, 4))
#' out23 = sanon(cbind(baseline, visit1, visit2, visit3, visit4) ~ grp(treatment, ref="P")
#'  + strt(center) + strt(sex) + covar(age), data=resp, P=P)
#' # each four visits
#' vcov(out23)
#'
vcov.sanon = function(object, ...){tmpVb = object$Vb; rownames(tmpVb) = colnames(tmpVb) = object$bnames; tmpVb}

#' Confidence Intervals for Model Parameters
#' 
#' Computes confidence intervals for one or more parameters in a fitted model.
#' 
#' Confidence intervals for adjusted parameters in the weighted least squares are computed based on an asymptotic normal.
#' 
#' @name confint.sanon
#' @aliases confint.sanon
#' @rdname confint.sanon
#' @method confint sanon
#' @docType methods
#' @export
#'
#' @param object,x an object of class "\code{sanon}", usually, a result of a call to \code{\link{sanon}}
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param ... further arguments passed to or from other methods.
#' @return \item{ci}{A matrix (or vector) with columns giving Mann-Whiteney estimates and their lower and upper confidence limits for each parameter with estimates. The interval will be labelled as Lower for (1 - level)/2 limit and Upper for 1 - (1 - level)/2 limit (by default 0.025 and 0.975). }
#' @return \item{level}{Confidence level}
#' @return \item{advarnames}{Adjust variable names in the weighted least squares method}
#'
#' @examples
#' ##### Example 3.1 Randomized Clinical Trial of Chronic Pain #####
#' data(cpain)
#' out1 = sanon(response ~ grp(treat, ref="placebo") + strt(center) + strt(diagnosis), data=cpain)
#' confint(out1)
#'
#' ##### Example 3.2 Randomized Clinical Trial of Respiratory Disorder #####
#' data(resp)
#' P = rbind(rep(0, 4), diag(4), rep(0, 4))
#' out23 = sanon(cbind(baseline, visit1, visit2, visit3, visit4) ~ grp(treatment, ref="P")
#'  + strt(center) + strt(sex) + covar(age), data=resp, P=P)
#' # each four visits
#' confint(out23)
#'

confint.sanon = function(object, parm = NULL, level = 0.95, ...)
{
b2 = coef(object)
se = sqrt(diag(vcov(object)))
za = qnorm(1-(1-level)/2)
TAB = cbind(b2 + 0.5, b2 - za * se + 0.5, b2 + za * se + 0.5)
colnames(TAB) = c("Estimate", "Lower", "Upper")
rownames(TAB) = object$bnames
if(!is.null(parm)) TAB = TAB[parm,]
res = list(call=object$call, ci=TAB, level=level, advarnames = object$advarnames)
class(res) = "confint.sanon"
res
}
#' @rdname confint.sanon
#' @method print confint.sanon
#' @export
print.confint.sanon = function(x, ...)
{
cat(paste("M-W Estimate and ", 100*x$level, "% Confidence Intervals \n", sep=""))
if(length(x$advarnames) > 0){ 
tmpadvarnames = x$advarnames
lim = max(which(cumsum(nchar(tmpadvarnames)) < 50))
if(lim < length(tmpadvarnames)){cat(paste("(adjusted by",  paste(c(tmpadvarnames[1:lim], ""), collapse=", "), "\n", paste(tmpadvarnames[(lim+1):length(tmpadvarnames)], collapse=", "), ")")) }else
{cat(paste("(adjusted by",  paste(tmpadvarnames, collapse=", "), ")"))}
}
cat(":\n")
print(round(x$ci, 4))
cat("\n")
}

#' Contrast for Model Parameters
#' 
#' Inference by contrast of parameters in a fitted model.
#' 
#' This function provide the inference based on contrast after applying the function \code{\link{sanon}}.
#' The contrast matrix C should be defined by the user. If the the number of row of C = 1, the confidence interval for the estimator is produced.
#' 
#' @name contrast
#' @aliases contrast
#' @rdname contrast
#' @docType methods
#' @export
#'
#' @param object,x an object of class "\code{sanon}", usually, a result of a call to \code{\link{sanon}}.
#' @param C contrast matrix. The number of column should be same as the length of \code{b} in outputs of \code{sanon}.
#' @param confint logical value for whether the confidence interval is computed (only if C has one row).
#' @param level the confidence level required (only if C has one row).
#' @param ... further arguments passed to or from other methods.
#' @return \item{C}{contrast matrix}
#' @return \item{Cb}{contrast estimates}
#' @return \item{VCb}{variance and covariance matrix of \code{Cb}}
#' @return \item{se}{standard error of \code{Cb}}
#' @return \item{level}{confidence level}
#' @return \item{UL}{upper confidence limit (only if the number of row of C = 1, otherwise \code{NULL})}
#' @return \item{LL}{lower confidence limit (only if the number of row of C = 1, otherwise \code{NULL})}
#' @return \item{Q}{test statistic}
#' @return \item{df}{degree of freedom}
#' @return \item{p}{p-value}
#'
#' @examples
#' ##### Example 3.2 Randomized Clinical Trial of Respiratory Disorder #####
#' data(resp)
#' P = rbind(rep(0, 4), diag(4), rep(0, 4))
#' out23 = sanon(cbind(baseline, visit1, visit2, visit3, visit4) ~ grp(treatment, ref="P")
#'  + strt(center) + strt(sex) + covar(age), data=resp, P=P)
#'
#' # Homogeneity of the xi_k across the four visits
#' contrast(out23, C=cbind(diag(3), rep(-1, 3)))
#'
#' # Comparison between treatments for the average of the xi_k across the 4 visits
#' contrast(out23, C=matrix(rep(1, 4)/4, ncol=4))
#'
contrast = function(object, C = diag(length(object$b)), confint = FALSE, level = 0.95, ...)
{
if(class(C) != "matrix") stop("C should be matrix")
if(ncol(C) != length(object$b)) stop("column length of C should be same as length of b")
if(nrow(C) != 1 & confint == TRUE) stop("Confidence interval is computed only if C has one row")
b2 = C %*% object$b
Vb2 = C %*% object$Vb %*% t(C)
Q2 = t(b2) %*% solve(Vb2) %*% b2
df2 = nrow(C)
za = qnorm(1-(1-level)/2)
Cb = C %*% (object$b+0.5)
UL = LL = NULL
if(df2 > 1){ se = sqrt(diag(Vb2))}else{ 
se = sqrt(Vb2)
if(confint){
UL = Cb + za * se
LL = Cb - za * se
}}
pQ2 = 1 - pchisq(Q2, df2)
colnames(C) = object$bnames
res = list(C=C, Cb=Cb, VCb=Vb2, se=se, level = level, UL=UL, LL=LL, Q=Q2, df=df2, p=pQ2)
class(res) = "contrast"
res
}

#' @rdname contrast
#' @method print contrast
#' @export
print.contrast = function(x, ...)
{
TAB = cbind(x$Q, x$df, x$p)
colnames(TAB) = c("Chisq", "df", "Pr(>Chisq)")
#rownames(TAB) = x$bnames
cat("Contrast Matrix:\n")
print(x$C)
cat("\n")
cat("Contrast Inference:\n")
printCoefmat(TAB, digits=3, dig.tst = 2, P.values = TRUE, has.Pvalue = TRUE)
cat("\n")
if(!is.null(x$UL))
{
cat(paste("Contrast M-W Estimate (", 100*x$level, "% Confidence Interval)\n", sep=""))
cat(paste(sprintf("%.4f", x$Cb), "(", sprintf("%.4f", x$LL), "-", sprintf("%.4f", x$UL), ")"))
cat("\n")
}
}

