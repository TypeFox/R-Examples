#' Graphical model stability and model selection procedures
#'
#' @name mplot-package
#' @docType package
#' @title Graphical model stability and model selection procedures
#' @keywords package
NULL


#' Body fat data set
#'
#' A data frame with 128 observations on 15 variables.
#'
#' @name bodyfat
#' @format A data frame with 128 observations on 15 variables.
#' \describe{
#' \item{Id}{Identifier}
#' \item{Bodyfat}{Bodyfat percentage}
#' \item{Age}{Age (years)}
#' \item{Weight}{Weight (kg)}
#' \item{Height}{Height (inches)}
#' \item{Neck}{Neck circumference (cm)}
#' \item{Chest}{Chest circumference (cm)}
#' \item{Abdo}{Abdomen circumference (cm) "at the umbilicus
#'            and level with the iliac crest"}
#' \item{Hip}{Hip circumference (cm)}
#' \item{Thigh}{Thigh circumference (cm)}
#' \item{Knee}{Knee circumference (cm)}
#' \item{Ankle}{Ankle circumference (cm)}
#' \item{Bic}{Extended biceps circumference (cm)}
#' \item{Fore}{Forearm circumference (cm)}
#' \item{Wrist}{Wrist circumference (cm) "distal to the
#'             styloid processes"}
#' }
#' @details A subset of the 252 observations available in the \code{mfp} package.
#'   The selected observations avoid known high leverage points and
#'   outliers.  The unused points from the data set could be used to validate
#'   selected models.
#' @docType data
#' @keywords datasets
#' @usage data(bodyfat)
#' @references Johnson W (1996, Vol 4). Fitting percentage of
#'   body fat to simple body measurements. Journal of Statistics
#'   Education. Bodyfat data retrieved from
#'   http://www.amstat.org/publications/jse/v4n1/datasets.johnson.html
#'   An expanded version is included in the \code{mfp} R package.
#' @examples
#' data(bodyfat)
#' full.mod = lm(Bodyfat~.,data=subset(bodyfat,select=-Id))
NULL

#' Rock-wallabies data set
#'
#' On Chalkers Top in the Warrumbungles (NSW, Australia) 200 evenly distributed
#' one metre squared plots were surveyed. Plots were placed at a density
#' of 7-13 per hectare. The presence or absence of fresh
#' (<1 month old) scats of rock-wallabies was recorded for each plot
#' along with location and a selection of predictor variables.
#'
#' @name wallabies
#' @format A data frame with 200 observations on 9 variables.
#' \describe{
#' \item{rw}{Presence of rock-wallaby scat}
#' \item{edible}{Percentage cover of edible vegetation}
#' \item{inedible}{Percentage cover of inedible vegetation}
#' \item{canopy}{Percentage canopy cover}
#' \item{distance}{Distance from diurnal refuge}
#' \item{shelter}{Whether or not a plot occurred within a shelter point (large
#'                rock or boulder pile)}
#' \item{lat}{Latitude of the plot location}
#' \item{long}{Longitude of the plot location}
#' }
#' @details Macropods defaecate randomly as they forage and scat 
#'   (faecal pellet) surveys are a reliable method for detecting the
#'   presence of rock-wallabies and other macropods. 
#'   Scats are used as an indication of spatial foraging patterns 
#'   of rock-wallabies and sympatric macropods. Scats deposited while
#'   foraging were not confused with scats deposited while
#'   resting because the daytime refuge areas of rock-wallabies
#'   were known in detail for each colony and no samples were
#'   taken from those areas. Each of the 200 sites were 
#'   examined separately to
#'   account for the different levels of predation risk and the
#'   abundance of rock-wallabies.
#' @docType data
#' @keywords datasets
#' @usage data(wallabies)
#' @references 
#'    Tuft KD, Crowther MS, Connell K, Mueller S and McArthur C (2011), 
#'    Predation risk and competitive interactions affect foraging of 
#'    an endangered refuge-dependent herbivore. Animal Conservation, 
#'    14: 447-457. doi: 10.1111/j.1469-1795.2011.00446.x
#' @examples
#' data(wallabies)
#' wdat = data.frame(subset(wallabies,select=-c(lat,long)), 
#'   EaD = wallabies$edible*wallabies$distance,
#'   EaS = wallabies$edible*wallabies$shelter,
#'   DaS = wallabies$distance*wallabies$shelter)
#' M1 = glm(rw~., family = binomial(link = "logit"), data = wdat)
NULL


#' Blood and other measurements in diabetics
#'
#' The diabetes data frame has 442 rows and 11 columns.
#' These are the data used in Efron et al. (2004).
#'
#' @name diabetes
#' @format A data frame with 442 observations on 11 variables.
#' \describe{
#' \item{age}{Age}
#' \item{sex}{Gender}
#' \item{bmi}{Body mass index}
#' \item{map}{Mean arterial pressure (average blood pressure)}
#' \item{tc}{Total cholesterol (mg/dL)? Desirable range: below 200 mg/dL}
#' \item{ldl}{Low-density lipoprotein ("bad" cholesterol)? 
#'            Desirable range: below 130 mg/dL }
#' \item{hdl}{High-density lipoprotein ("good" cholesterol)? 
#'            Desirable range: above 40 mg/dL}
#' \item{tch}{Blood serum measurement}
#' \item{ltg}{Blood serum measurement}
#' \item{glu}{Blood serum measurement (glucose?)}
#' \item{y}{A quantitative measure of disease progression 
#'          one year after baseline}
#' }
#' @details Data sourced from http://web.stanford.edu/~hastie/Papers/LARS
#' @docType data
#' @keywords datasets
#' @usage data(diabetes)
#' @references Efron, B., Hastie, T., Johnstone, I., Tibshirani, R., (2004).
#'   Least angle regression. The Annals of Statistics 32(2) 407-499.
#'   DOI: 10.1214/009053604000000067
#' @examples
#' data(diabetes)
#' full.mod = lm(y~.,data=diabetes)
NULL




#' Artificial example
#'
#' An artificial data set which causes stepwise regression
#' procedures to select a non-parsimonious model.
#' The true model is a simple linear regression of
#' y against x8.
#'
#' @name artificialeg
#' @format A data frame with 50 observations on 10 variables.
#' @details Inspired by the pathoeg data set in the MPV pacakge.
#' @docType data
#' @keywords datasets
#' @usage data(artificialeg)
#' @examples
#' data(artificialeg)
#' full.mod = lm(y~.,data=artificialeg)
#' step(full.mod)
#' \dontrun{
#' # generating model
#' n=50
#' set.seed(8) # a seed of 2 also works
#' x1 = rnorm(n,0.22,2)
#' x7 = 0.5*x1 + rnorm(n,0,sd=2)
#' x6 = -0.75*x1 + rnorm(n,0,3)
#' x3 = -0.5-0.5*x6 + rnorm(n,0,2)
#' x9 = rnorm(n,0.6,3.5)
#' x4 = 0.5*x9 + rnorm(n,0,sd=3)
#' x2 = -0.5 + 0.5*x9 + rnorm(n,0,sd=2)
#' x5 = -0.5*x2+0.5*x3+0.5*x6-0.5*x9+rnorm(n,0,1.5)
#' x8 = x1 + x2 -2*x3 - 0.3*x4 + x5 - 1.6*x6 - 1*x7 + x9 +rnorm(n,0,0.5)
#' y = 0.6*x8 + rnorm(n,0,2)
#' artificialeg = round(data.frame(x1,x2,x3,x4,x5,x6,x7,x8,x9,y),1)
#' }
NULL




#' Forced Expiratory Volume
#'
#' This data set consists of 654 observations on youths aged 3 to 19 from 
#' East Boston recorded duing the middle to late 1970's. 
#' Forced expiratory volume (FEV), a measure of lung capacity, is the 
#' variable of interest. Age and height are two continuous predictors. 
#' Sex and smoke are two categorical predictors.
#'
#' @name fev
#' @format A data frame with 654 observations on 5 variables.
#' \describe{
#' \item{age}{Age (years)}
#' \item{fev}{Forced expiratory volume (liters).  Roughly the amount 
#'            of air an individual can exhale in the first second of 
#'            a forceful breath.}
#' \item{height}{Height (inches).}
#' \item{sex}{Female is 0. Male is 1.}
#' \item{smoke}{A binary variable indicating whether or not the 
#'              youth smokes. Nonsmoker is 0. Smoker is 1.}
#' }
#' @details Copies of this data set can also be found in the 
#'  \code{coneproj} and \code{tmle} packages.
#' @references 
#'  Tager, I. B., Weiss, S. T., Rosner, B., and Speizer, F. E. (1979). 
#'  Effect of parental cigarette smoking on pulmonary function in children. 
#'  \emph{American Journal of Epidemiology}, \bold{110}, 15-26.
#'  
#'  Rosner, B. (1999).
#'  \emph{Fundamentals of Biostatistics}, 5th Ed., Pacific Grove, CA: Duxbury.
#'   
#'   Kahn, M.J. (2005). An Exhalent Problem for Teaching Statistics.
#'   \emph{Journal of Statistics Education},  \bold{13}(2). 
#'    http://www.amstat.org/publications/jse/v13n2/datasets.kahn.html
#' @docType data
#' @keywords datasets
#' @usage data(fev)
#' @examples
#' data(fev)
#' full.mod = lm(fev~.,data=fev)
#' step(full.mod)
NULL




#' Extract model elements
#'
#' This function extracts things like the formula,
#' data matrix, etc. from a lm or glm object
#'
#' @param model a fitted 'full' model, the result of a call
#'   to lm or glm (and in the future lme or lmer).
#' @param screen logical, whether or not to perform an initial
#'   screen for outliers.  Highly experimental, use at own risk.
#'   Default = FALSE.
#' @param redundant logical, whether or not to add a redundant
#'   variable.  Default = TRUE.
#' @noRd
mextract = function(model, screen = FALSE, redundant = TRUE){
  # what's the name of the dependent variable?
  yname = deparse(stats::formula(model)[[2]])
  # Set up the data frames for use
  data = stats::model.frame(model)
  X = stats::model.matrix(model)
  n = nrow(X)
  # full model plus redundant variable
  exp.vars = names(model$coefficients)[names(model$coefficients) != "(Intercept)"]

  if (redundant) {
    REDUNDANT.VARIABLE = stats::runif(n, min = 0, max = 1)
    X = cbind(X,REDUNDANT.VARIABLE)
    data = cbind(data,REDUNDANT.VARIABLE)
    exp.vars = c(exp.vars,"REDUNDANT.VARIABLE")
  }
  if (colnames(X)[1] == "(Intercept)") {
    # overwrite intercept with y-variable
    X[,1] = stats::model.frame(model)[,yname]
  } else {
    X = cbind(stats::model.frame(model)[,yname],X)
  }
  colnames(X)[1] = yname
  X = data.frame(X)
  fixed = stats::as.formula(c(paste(yname, "~"),
                       paste(colnames(X)[-1], collapse = "+")))
  Xy = X[c(2:ncol(X),1)]

  k = length(exp.vars) + 1 # +1 for intercept
  if (screen) {
    if (!requireNamespace("mvoutlier", quietly = TRUE)) {
      stop("mvoutlier package needed when screen=TRUE. Please install it.",
           call. = FALSE)
    }
    x.mad = apply(Xy, 2, stats::mad)
    Xy.sub = Xy[,which(x.mad != 0)]
    Xy = Xy[mvoutlier::pcout(Xy.sub)$wfinal01 == 1,]
    n = dim(Xy)[1]
    if (k >= n) {
      warning("Screening deleted too many observations.")
      return()
    }
  }
  wts = model$weights
  if (is.element("glm",class(model))) {
    wts = model$prior.weights
    Xy[,yname] = model$y
  }
  if (is.null(wts)) {
    wts = rep(1,n)
  }

  return(list(yname = yname, fixed = fixed,
              wts = wts, X = Xy, k = k,
              n = n, exp.vars = exp.vars,
              data = data, family = stats::family(model)))

  # MIXED MODELS NOT IMPLEMENTED IN FIRST RELEASE
  #   if(class(model)=="lmerMod"){ # lme4 package
  #     X = stats::model.matrix(model)
  #     data = data.frame(stats::model.frame(model)) # or slot(model, "frame")
  #     if(isGLMM(model)){
  #       family = stats::family(model) #
  #       # pass to glm ??
  #       # i.e. wild bootstrap appropriate for glmm's too?
  #       # though should really call glmer for the full model
  #       # so it may be the case that class(model) won't be lmerMod
  #       # if the family argument was passed to it as lmer()
  #       # calls glmer() if there is a family argument present
  #     }
  #     yname = deparse(stats::formula(model)[[2]])
  #     X = data.frame(data[,yname],X)
  #     colnames(X)[1] = yname
  #     fixed.formula = paste(yname,"~",
  #                           paste(names(fixef(model))[-1],collapse="+"))
  #     model = lm(fixed.formula,data = X)
  #   } else if(class(model)=="lme"){ # nlme package
  #     data = model$data
  #     X = stats::model.matrix(model,data = data)
  #     # nlme package assumes gaussian errors
  #     yname = deparse(stats::formula(model)[[2]])
  #     X = data.frame(data[,yname],X)
  #     colnames(X)[1] = yname
  #     fixed.formula = paste(yname,"~",
  #                           paste(names(stats::coef(model))[-1],collapse="+"))
  #     model = lm(fixed.formula, data = X)
  #   }
}

#' Safe deparse
#'
#' Supports long formula construction
#'
#' @param expr expression to be safely deparsed
#'
#' @noRd
safeDeparse <- function(expr){
  ret <- paste(deparse(expr), collapse = "")
  #rm whitespace
  gsub("[[:space:]][[:space:]]+", " ", ret)
}

#' Print text for fence methods
#'
#' This function provides the text for the case when trace=TRUE
#' when using lmfence and glmfence functions.
#'
#' @param score realised value
#' @param UB upper bound
#' @param obj fitted model object
#' @keywords internal
txt.fn = function(score,UB,obj){
  cat("\n")
  cat(paste("hatQm:", round(score,2),"; Upper bound:", round(UB,2)),"\n")
  cat(paste("hatQm <= UB:",score<=UB,"\n"))
  cat(deparse(stats::formula(obj)))
  cat("\n")
}

#' Process results within af function
#'
#' This function is used by the af function to process
#' the results when iterating over different boundary values
#'
#' @param fence.mod set of fence models
#' @param fence.rank set of fence model ranks
#' @noRd
process.fn = function(fence.mod,fence.rank){
  del2 = function(x) x[-c(1:2)]
  splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))

  # best only (true fence)
  fence.mod.bo = fence.mod[fence.rank==1]
  #  temp.bo = sort(table(sapply(fence.mod.bo,deparse,
  #                              width.cutoff=500)), decreasing=TRUE)
  temp.bo = sort(table(unlist(lapply(lapply(fence.mod.bo,as.character),del2))),
                 decreasing=TRUE)
  pstarj.bo = as.numeric(temp.bo[1]/length(fence.mod.bo))
  pstarnamej.bo = names(temp.bo)[1]

  # all that pass the fence
  #temp.all = sort(table(sapply(fence.mod,deparse,
  #                             width.cutoff=500)),decreasing=TRUE)
  temp.all = sort(table(unlist(lapply(lapply(fence.mod,as.character),del2))),
                  decreasing=TRUE)
  #old version
  #pstarj.all = as.numeric(temp.all[1]/length(fence.mod))
  #pstarnamej.all = names(temp.all)[1]
  # new version
  unlist.fence.rank = unlist(fence.rank)
  fence.rank.split = splitAt(unlist.fence.rank,which(unlist.fence.rank==1))
  custom.p = 1/unlist(lapply(fence.rank.split,length))
  custom.names = unlist(lapply(lapply(fence.mod.bo,as.character),del2))
  agg = stats::aggregate(custom.p,by=list(custom.names),sum)
  agg = agg[order(agg$x,decreasing = TRUE),]
  pstarj.all = agg[1,2]/sum(unlist.fence.rank==1)
  pstarnamej.all = agg[1,1]

  return(c(pstarj.bo,pstarnamej.bo,pstarj.all,pstarnamej.all))
}
