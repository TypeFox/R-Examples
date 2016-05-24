##' Given covariates, a treatment variable, and a stratifying factor,
##' calculates standardized mean differences along each covariate,
##' with and without the stratification and tests for conditional
##' independence of the treatment variable and the covariates within
##' strata.
##'
##' In the unstratified case, the standardized difference of covariate
##' means is the mean in the treatment group minus the mean in the
##' control group, divided by the S.D. (standard deviation) in the
##' same variable estimated by pooling treatment and control group
##' S.D.s on the same variable.  In the stratified case, the
##' denominator of the standardized difference remains the same but
##' the numerator is a weighted average of within-stratum differences
##' in means on the covariate.  By default, each stratum is weighted
##' in proportion to the harmonic mean \eqn{1/[(1/a +
##' 1/b)/2]=2*a*b/(a+b)} of the number of treated units (a) and
##' control units (b) in the stratum; this weighting is optimal under
##' certain modeling assumptions (discussed in Kalton 1968, Hansen and
##' Bowers 2008).  This weighting can be modified using the
##' \code{stratum.weights} argument; see below.
##'
##' When the treatment variable, the variable specified by the
##' left-hand side of \code{fmla}, is not binary, \code{xBalance}
##' calculates the covariates' regressions on the treatment variable,
##' in the stratified case pooling these regressions across strata
##' using weights that default to the stratum-wise sum of squared
##' deviations of the treatment variable from its stratum mean.
##' (Applied to binary treatment variables, this recipe gives the same
##' result as the one given above.)  In the numerator of the
##' standardized difference, we get a ``pooled S.D.'' from separating
##' units into two groups, one in which the treatment variable is 0 or
##' less and another in which it is positive.  If \code{report}
##' includes "adj.means", covariate means for the former of these
##' groups are reported, along with the sums of these means and the
##' covariates' regressions on either the treatment variable, in the
##' unstratified (``pre'') case, or the treatment variable and the
##' strata, in the stratified (``post'') case.
##'
##' \code{stratum.weights} can be either a function or a numeric
##' vector of weights.  If it is a numeric vector, it should be
##' non-negative and it should have stratum names as its names. (i.e.,
##' its names should be equal to the levels of the factor specified by
##' \code{strata}.) If it is a function, it should accept one
##' argument, a data frame containing the variables in \code{data} and
##' additionally \code{Tx.grp} and \code{stratum.code}, and return a
##' vector of non-negative weights with stratum codes as names; for an
##' example, do \code{getFromNamespace("harmonic", "RItools")}.
##'
##' If \code{covariate.scaling} is not \code{NULL}, no scaling is
##' applied. This behavior is likely to change in future versions.
##' (If you want no scaling, set \code{covariate.scaling=1}, as this
##' is likely to retain this meaning in the future.)
##'
##' \code{adj.mean.diffs.null.sd} returns the standard deviation of
##' the Normal approximated randomization distribution of the
##' strata-adjusted difference of means under the strict null of no
##' effect.
##' @title Standardized Differences for Stratified Comparisons
##' @param fmla A formula containing an indicator of treatment
##'   assignment on the left hand side and covariates at right.
##' @param strata A list of right-hand-side-only formulas containing
##'   the factor(s) identifying the strata, with \code{NULL} entries
##'   interpreted as no stratification; or a factor with length equal
##'   to the number of rows in data; or a data frame of such
##'   factors. See below for examples.
##' @param data A data frame in which \code{fmla} and \code{strata}
##'   are to be evaluated.
##' @param report Character vector listing measures to report for each
##'   stratification; a subset of \code{c("adj.means",
##'   "adj.mean.diffs", "adj.mean.diffs.null.sd", "chisquare.test",
##'   "std.diffs", "z.scores", "p.values", "all")}. P-values reported
##'   are two-sided for the null-hypothesis of no effect. The option
##'   "all" requests all measures.
##' @param stratum.weights Weights to be applied when aggregating
##'   across strata specified by \code{strata}, defaulting to weights
##'   proportional to the harmonic mean of treatment and control group
##'   sizes within strata.  This can be either a function used to
##'   calculate the weights or the weights themselves; if
##'   \code{strata} is a data frame, then it can be such a function, a
##'   list of such functions, or a data frame of stratum weighting
##'   schemes corresponding to the different stratifying factors of
##'   \code{strata}.  See details.
##' @param na.rm Whether to remove rows with NAs on any variables
##'   mentioned on the RHS of \code{fmla} (i.e. listwise deletion).
##'   Defaults to \code{FALSE}, wherein rows aren't deleted but for
##'   each variable with \code{NA}s a missing-data indicator variable
##'   is added to the variables on which balance is calculated and
##'   medians are imputed for the variable with missing data (in
##'   RItools versions 0.1-9 and before the default imputation was the
##'   mean, in RItools versions 0.1-11 and henceforth the default is
##'   the median). See the example below.
##' @param covariate.scaling A scale factor to apply to covariates in
##'   calculating \code{std.diffs}.  If \code{NULL}, \code{xBalance}
##'   pools standard deviations of each variable in the treatment and
##'   control group (defining these groups according to whether the
##'   LHS of \code{formula} is greater than or equal to 0).  Also, see
##'   details.
##' @param normalize.weights If \code{TRUE}, then stratum weights are
##'   normalized so as to sum to 1.  Defaults to \code{TRUE}.
##' @param impfn A function to impute missing values when
##'   \code{na.rm=FALSE}. Currently \code{\link{median}}. To impute
##'   means use \code{\link{mean.default}}.
##' @param post.alignment.transform Optional transformation applied to
##'   covariates just after their stratum means are subtracted off.
##' @return An object of class \code{c("xbal", "list")}.  There are
##'   \code{plot}, \code{print}, and \code{xtable} methods for class
##'   \code{"xbal"}; the \code{print} method is demonstrated in the
##'   examples.
##' @note Evidence pertaining to the hypothesis that a treatment
##'   variable is not associated with differences in covariate values
##'   is assessed by comparing the differences of means (or regression
##'   coefficients), without standardization, to their distributions
##'   under hypothetical shuffles of the treatment variable, a
##'   permutation or randomization distribution.  For the unstratified
##'   comparison, this reference distribution consists of differences
##'   (more generally, regression coefficients) when the treatment
##'   variable is permuted without regard to strata.  For the
##'   stratified comparison, the reference distribution is determined
##'   by randomly permuting the treatment variable within strata, then
##'   re-calculating the treatment-control differences (regressions of
##'   each covariate on the permuted treatment variable). Significance
##'   assessments are based on the large-sample Normal approximation
##'   to these reference distributions.
##' @export
##' @references Hansen, B.B. and Bowers, J. (2008), ``Covariate
##'   Balance in Simple, Stratified and Clustered Comparative
##'   Studies,'' \emph{Statistical Science} \bold{23}.
##'
##'   Kalton, G. (1968), ``Standardization: A technique to control for
##'   extraneous variables,'' \emph{Applied Statistics} \bold{17},
##'   118--136.
##' @author Ben Hansen and Jake Bowers and Mark Fredrickson
##' @keywords design nonparametric
##' @import SparseM svd
##' @examples
##' data(nuclearplants)
##' ##No strata, default output
##' xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
##'          data=nuclearplants)
##'
##' ##No strata, all output
##' xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
##'          data=nuclearplants,
##'          report=c("all"))
##'
##' ##Stratified, all output
##' xBalance(pr~.-cost-pt, strata=factor(nuclearplants$pt),
##'          data=nuclearplants,
##'          report=c("adj.means", "adj.mean.diffs",
##'                   "adj.mean.diffs.null.sd",
##'                   "chisquare.test", "std.diffs",
##'                   "z.scores", "p.values"))
##'
##' ##Comparing unstratified to stratified, just adjusted means and
##' #omnibus test
##' xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
##'          strata=list(unstrat=NULL, pt=~pt),
##'          data=nuclearplants,
##'          report=c("adj.means", "chisquare.test"))
##'
##' ##Comparing unstratified to stratified, just adjusted means and
##' #omnibus test
##' xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
##'          strata=data.frame(unstrat=factor('none'),
##'            pt=factor(nuclearplants$pt)),
##'          data=nuclearplants,
##'          report=c("adj.means", "chisquare.test"))
##'
##' ##Missing data handling.
##' testdata<-nuclearplants
##' testdata$date[testdata$date<68]<-NA
##'
##' ##na.rm=FALSE by default
##' xBalance(pr ~ date, data = testdata, report="all")
##' xBalance(pr ~ date, data = testdata, na.rm = TRUE,report="all")
##'
##' ##To match versions of RItools 0.1-9 and older, impute means
##' #rather than medians.
##' ##Not run, impfn option is not implemented in the most recent version
##' \dontrun{xBalance(pr ~ date, data = testdata, na.rm = FALSE,
##'            report="all", impfn=mean.default)}
##'
##' ##Comparing unstratified to stratified, just one-by-one wilcoxon
##' #rank sum tests and omnibus test of multivariate differences on
##' #rank scale.
##' xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
##'          strata=data.frame(unstrat=factor('none'),
##'            pt=factor(nuclearplants$pt)),
##'          data=nuclearplants,
##'          report=c("adj.means", "chisquare.test"),
##' 	 post.alignment.transform=rank)
xBalance <- function(fmla, strata=list(unstrat=NULL),
                     data,
                     report=c("std.diffs","z.scores","adj.means","adj.mean.diffs","adj.mean.diffs.null.sd",
                         "chisquare.test","p.values", "all")[1:2],
                     #                     include.means=FALSE, chisquare.test=FALSE,
                     stratum.weights=harmonic, na.rm=FALSE,
                     covariate.scaling=NULL, normalize.weights=TRUE,impfn=median,
                     post.alignment.transform=NULL) {
  stopifnot(class(fmla)=="formula",
            is.null(strata) || is.factor(strata) || is.list(strata),
            !is.data.frame(strata) || !any(is.na(names(strata))),
            !is.data.frame(strata) || all(names(strata)!=""),
            !is.data.frame(strata) || all(sapply(strata, is.factor)),
            is.null(data) || is.data.frame(data),
            is.null(post.alignment.transform) || is.function(post.alignment.transform)
            )

  if (any(grepl("strata", fmla))) {
    splitstrat <- findStrata(fmla, data)

    if (!is.null(splitstrat$strata)) {
      fmla <- splitstrat$newx

      # apply was giving trouble here; not ideal but we shouldn't
      # be having more than a few strata, so shouldn't be a
      # performance hit
      strata <- list()
      for (i in paste("~", splitstrat$strata)) {
        strata <- c(strata, list(formula(i)))
      }

      names(strata) <- splitstrat$strata

      # Automatically add the unadjusted version. Maybe make this an optional argument later.
      strata <- c(list("Unadj" = NULL), strata)
    }
  }

  # Using charmatch instead of pmatch to distinguish between no match and ambiguous match. It reports
  # -1 for no match, and 0 for ambiguous (multiple) matches.
  valid.for.report <- c("adj.means","adj.mean.diffs","adj.mean.diffs.null.sd","chisquare.test",
                                     "std.diffs","z.scores","p.values","all")
  report.good <- charmatch(report, valid.for.report, -1)
  if (any(report.good == -1)) {
    stop(paste("Invalid option(s) for report:", paste(report[report.good == -1], collapse=", ")))
  }
  if (any(report.good == 0)) {
    stop(paste("Option(s) for report match multiple possible values:", paste(report[report.good == 0], collapse=", ")))
  }

  # Now that we've found the partial matches, get their proper names
  report <- valid.for.report[report.good]

  if (is.null(strata))
    warning("Passing NULL as a 'strata=' argument is depracated;\n for balance w/o stratification pass 'list(nostrat=NULL)' instead.\n (Or did you mean to pass a non-NULL 'strata=' argument? Then check for typos.)")

  if (is.list(strata) && !is.data.frame(strata) && !all(sapply(strata, function(x) (is.null(x) | inherits(x,"formula")))))
    stop("For balance against multiple alternative stratifications,\n please make 'strata' either a data frame or a list containing formulas or NULL entries.")

  if("all" %in% report)
    report<-c("adj.means","adj.mean.diffs","adj.mean.diffs.null.sd","chisquare.test", "std.diffs","z.scores","p.values")

  ### NA Handling ##
  if (na.rm==TRUE) {
    tfmla <- terms.formula(fmla,data=data, keep.order=TRUE)
  } else {
    data <- naImpute(fmla,data,impfn)
    tfmla <- attr(data, 'terms')
  }
  ### End NA handling ###

  ###Extract the treatment var
  if (!attr(tfmla, "response")>0)
    stop("fmla must specify a treatment group variable")

  zz <- eval(tfmla[[2]], data, parent.frame()) # changed for v.93, see comment in log
  zzname<- deparse(tfmla[[2]])
  if (!is.numeric(zz) & !is.logical(zz))
    stop("LHS of fmla should be logical or numeric")
  if (any(is.na(zz)))
    stop('NAs on LHS of fmla not allowed.')
  ### End extract treatment var

  mm1 <- xBalance.makeMM(tfmla,data)

  ### Prepare ss.df, data frame of strata
  if (is.null(strata))
    ss.df <- data.frame(unstrat=factor(numeric(length(zz))))

  if (is.factor(strata) & length(strata)!=length(zz))
    stop("length of strata doesn\'t match dim of data")

  if (is.factor(strata))
    ss.df <- data.frame(strat=factor(strata))
  if (is.data.frame(strata))
    ss.df <- as.data.frame(lapply(strata,factor))
  if (is.list(strata) & !is.data.frame(strata)) {
    ### In this case strata should be a list of formulas

    pfr <- parent.frame()
    ss.df <-
      lapply(strata,
             function(fmla) {
               if (is.null(fmla)) factor(numeric(length(zz))) else {
                 ss <- eval(attr(terms(fmla), "variables"), data,
                            pfr)
                 if (length(ss)-1) interaction(ss, drop=TRUE) else factor(ss[[1]])
               }
             })
    ss.df <- as.data.frame(ss.df)
  }
  ### End prepare ss.df, data frame of strata

  ### Remove stratification variables without levels (e.g., all NAs), with a warning
  if (any(ss.rm <- !sapply(ss.df, nlevels))) {
    if (length(ss.df)==1)
      stop("'strata=' variable contains no strata.  Perhaps it evaluates to NAs?")
    if (all(ss.rm))
      stop("'strata=' variables contain no strata.  Perhaps they all evaluate to NAs?")
    ss.rm.nms <- if (is.null(names(ss.df))) which(ss.rm) else names(ss.df)[ss.rm]
    ss.rm.nms <- paste(ss.rm.nms, collapse=" ,")
    warning(paste("Removing the following strata entries, which contained no strata.\n(Perhaps they evaluate to NAs?)\n",
                  ss.rm.nms))
    ss.df <- ss.df[!ss.rm]
  }
  ### End remove stratification variables without levels
  gs.df <- xBalance.find.goodstrats(ss.df,zz,mm1)

  swt.ls <- xBalance.make.stratwts(stratum.weights,ss.df, gs.df, zz, data, normalize.weights)

  s.p <- if (is.null(covariate.scaling)) {
    xBalance.makepooledsd(zz,mm1,dim(mm1)[1])
  } else 1

  ### Call xBalanceEngine here.

  RES <- lapply(names(ss.df),
                function(nm) {
                  ###                  workingswt.ls<-swt.ls[[nm]]  # shouldn't be neccessary after r216 change to xBalance.make.stratwts
                  ###                  workingswt.ls[["wtratio"]]<-swt.ls[[nm]][["wtratio"]][gs.df[[nm]]]
                  xBalanceEngine(factor(ss.df[gs.df[[nm]],nm]),
                                 zz[gs.df[[nm]]],
                                 mm1[gs.df[[nm]],,drop=FALSE],
                                 report, swt.ls[[nm]],
                                 s.p, normalize.weights,zzname,
                                 post.alignment.transform)
                })
  names(RES) <- names(ss.df)
  ##nms <- paste(rep(names(ss.df), rep(length(RES[[1]]$dfr),length(ss.df))),
  ##            names(RES[[1]]$dfr), sep=".")
  ans <- list() ##the overall function still returns a list because of the overall test info.
  ##results is an array of variables by balance statistics by stratification.
  ##here assuming that the variables and statistics are the same across stratifications (including unstratified).
  ans$results<-array(dim=c(vars=nrow(RES[[1]][["dfr"]]),stat=ncol(RES[[1]][["dfr"]]),strata=length(RES)),
                     dimnames=list(vars=rownames(RES[[1]][["dfr"]]),stat=colnames(RES[[1]][["dfr"]]),strata=names(RES)))

  attr(ans$results, "originals") <- attr(mm1, "originals")

  for (i in names(RES)) {
    ##print(i);print(RES[[i]][["dfr"]])
    ans$results[,,i]<-as.matrix(RES[[i]][["dfr"]])
  }
  ##dimnames(ans)[["stat"]][grep("Tx",dimnames(ans)[["stat"]])]<-c("adj.mean.strata=0","adj.mean.strata=1")
  ##ans$by.variable <- do.call(cbind, lapply(RES, function(x) x[['dfr']]) )
  ##colnames(ans$by.variable) <- nms
  attr(ans, "fmla") <- formula(tfmla)

  if ("chisquare.test" %in% report) {
    ans$overall <- data.frame(chisquare = numeric(length(RES)),
                              df        = numeric(length(RES)),
                              p.value   = numeric(length(RES)),
                              row.names = names(RES))
    for (nn in names(RES)) {
      ans$overall[nn,'chisquare'] <- RES[[nn]]$chisq['chisquare']
      ans$overall[nn,'df']        <- RES[[nn]]$chisq['df']
      ans$overall[nn,'p.value']   <- pchisq(RES[[nn]]$chisq['chisquare'],
                                            df = RES[[nn]]$chisq['df'],
                                            lower.tail = FALSE)

    }

    attr(ans$overall, "tcov") <- lapply(RES, function(r) {
      r$tcov
    })
  }
  class(ans) <- c("xbal", "list")
  ans
}

xBalance.make.stratum.mean.matrix <- function(ss, mm) {

  post.nobs <- dim(mm)[1]
  nlev <- nlevels(ss)

  # for this matrix, finding the indices of the rows is easy, as there is only one
  # item per row, and there post.nobs number of rows.
  tR <- new("matrix.csr",
            ja = as.integer(as.integer(ss)),
            ia = as.integer(1:(post.nobs+ 1)),
            ra = unsplit(1/tapply(ss,ss,length),ss),
            dimension = c(post.nobs,nlev))

  # With many items per row, we need to break the list of strata
  # down into indices of where each row starts and ends
  # e.g. Say ss = 0 0 0 0 1 1 1 1 1 1 1 2 2 2 2 the row indices would be
  # 1 5 12 16 (where 16 is the start of the non existant 4th row)

  L <- new("matrix.csr", #ifelse(oldver,"tripletMatrix","dgTMatrix"),
           ia = as.integer(1:(post.nobs + 1)),
           ja = as.integer(as.integer(ss)),
           ra = rep(1,length(ss)),
           dimension = c(post.nobs,nlev))


  msmn <- t(tR) %*% mm
  msmn <- L %*% msmn

  msmn <- as.matrix(msmn)

  return(msmn)
}


# Extract `strata(...)` arguments from a formula.
findStrata <- function(x, data) {

  t <- terms(x, specials = "strata", data = data)

  strata <- rownames(attr(t, "factors"))[attr(t, "specials")$strata]
  if (length(strata) > 0) {
    # Trying to update(x) directly was causing errors about having a "."
    # and no data. Updating the terms returns a fmla and bypasses the bug.
    x <- update(terms(x, data=data),
                as.formula(paste("~ . - ", paste(strata, collapse="-"))))

    # The gsubs return only the `...` inside `strata(...)`
    return(list(newx = x,
                strata = gsub("\\)", "", gsub("strata\\(", "", strata))))
  }

  return(list(newx = x, strata = NULL))
}
