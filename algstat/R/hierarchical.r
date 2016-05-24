#' Fitting Hierarchical Log-linear Models with Algebraic Methods
#'
#' Run the Metropolis-Hastings algorithm using a Markov basis computed with 4ti2 to sample from the conditional distribution of the data given the sufficient statistics of a hierarchical model.
#'
#' Hierarchical fits and tests a hierarchical log-linear model on a contingency table.  In many ways, hierarchical is like stats::loglin or MASS::loglm; however, there are a few key differences in the functionality of hierarchical.
#'
#' The first difference is methodological.  The tests conducted with hierarchical are exact tests based on the conditional distribution of the data given the sufficient statistics for the model.  In other words, they are Fisher's exact test analogues for log-linear models.  These tests are made possible by advances in algebraic statistics; see the first and second references below.  In particular, hierarchical leverages Markov bases through the software 4ti2 to construct a Metropolis-Hastings algorithm to sample from the conditional distribution of the table given the sufficient statistics.
#'
#' A second way that hierarchical differs from stats::loglin or MASS::loglm is in generalizing the kinds of tests performed.  While those allow for the asymptotic unconditional testing using Pearson's X^2 test and the likelihood ratio test (MASS::loglm is simply a wrapper for stats::loglin), hierarchical gives several test statistics: Pearson's X^2, the likelihood ratio G^2, Freeman-Tukey, Cressie-Read (lambda = 2/3), and Neyman's modified X^2., see the last reference.  In other words, to compute the exact p-value, iter = 1e4 samples are sampled from the conditional distribution of the table given the sufficient statistics, and then the proportion of tables that have X^2, G^2, etc. values greater than or equal to that of the observed table is the p value for the (conditional) exact test.  A similar, but perhaps preferable approach, simply adds up the probabilities of the tables that have probabilities less than or equal to that of the observed table; this is the first line output in hierarchical and does not use a test statistic.
#'
#' Some authors (see the third reference) suggest that for discrete problems, a "mid p value" is preferable to the traditional p value, and when presented should be interepreted in the same way.  If the p value is defined to be, say, P(samps >= obs), the mid p value is defined to be P(samps > obs) + P(samps == obs)/2.  The mid p value is computed for each test.
#'
#' Since the tests make use of Monte Carlo sampling, standard errors (SE) are reported for each statistic.  For the test statistics, this is just the standard deviation of the samples divided by the square root of the sample size, iter; they are computed and returned by the print method.  The standard errors of the p values use the CLT asymptotic approximation and, therefore, warrant greater consideration when the p value is close to 0 or 1.  
#' 
#' @param formula formula for the hierarchical log-linear model
#' @param data data, typically as a table but can be in different formats.  see \code{\link{teshape}}
#' @param iter number of chain iterations
#' @param burn burn-in
#' @param thin thinning
#' @param engine C++ or R? (C++ yields roughly a 20-25x speedup)
#' @param method should the expected value (exp) be fit using iterative proportional fitting (via loglin) or the MCMC as the average of the steps?
#' @param moves the markov moves for the mcmc
#' @return a list containing named elements
#' \itemize{
#'   \item \code{steps}: an integer matrix whose columns represent individual samples from the mcmc.
#'   \item \code{moves}: the moves used for the proposal distribution in the mcmc, computed with 4ti2 (note that only the positive moves are given).
#'   \item \code{acceptProb}: the average acceptance probability of the moves, including the thinned moves.
#'   \item \code{param}: the fitted parameters of the log linear model.
#'   \item \code{df}: parameters per term in the model
#'   \item \code{quality}: model selection statistics AIC, AICc, and BIC.
#'   \item \code{residuals}: the (unstandardized) pearson residuals (O - E) / sqrt(E)
#'   \item \code{call}: the call.
#'   \item \code{obs}: the contingency table given.
#'   \item \code{exp}: the fit contingency table as an integer array.
#'   \item \code{A}: the sufficient statistics computing matrix (from Tmaker).
#'   \item \code{p.value}: the exact p-values of individual tests, accurate to Monte-Carlo error.  these are computed as the proportion of samples with statistics equal to or larger than the oberved statistic.
#'   \item \code{mid.p.value}: the mid p.values, see Agresti pp.20--21.
#'   \item \code{statistic}: the pearson's chi-squared (X2), likelihood ratio (G2), Freeman-Tukey (FT), Cressie-Read (CR), and Neyman modified chi-squared (NM) statistics computed for the table given.
#'   \item \code{sampsStats}: the statistics computed for each mcmc sample.
#'   \item \code{cells}: the number of cells in the table.
#'   \item \code{method}: the method used to estimate the table.
#' }	
#' @export hierarchical
#' @author David Kahle
#' @seealso \code{\link{loglin}}, \code{\link{loglm}}, \code{\link{metropolis}}
#' @references Diaconis, P. and B. Sturmfels (1998). Algebraic Algorithms for Sampling from Conditional Distributions. \emph{The Annals of Statistics} 26(1), pp.363-397.
#' @references Drton, M., B. Sturmfels, and S. Sullivant (2009). \emph{Lectures on Algebraic Statistics}, Basel: Birkhauser Verlag AG.
#' @references Agresti, A. (2002). \emph{Categorical Data Analysis}, Basel: John Wiley & Sons, 2ed.
#' @references Agresti, A. (1992). A Survey of Exact Inference for Contingency Tables \emph{Statistical Science} 7(1), pp.131-153.
#' @references Read, T. and Cressie, N. (1998). \emph{Goodness-of-Fit Statistics for Discrete Multivariate Data}, Springer-Verlag.
#' @examples
#'
#' \dontrun{
#'
#' 
#' ## handedness introductory example
#' ############################################################
#' 
#' data(handy)
#' 
#' (out <- hierarchical(~ Gender + Handedness, data = handy))
#' 
#' # hierarchical performs the same tasks as loglin and loglm,
#' # but hierarchical gives the exact test p values and more statistics
#' statsFit <- stats::loglin(handy, list(c(1),c(2)), fit = TRUE, param = TRUE)
#' massFit <- MASS::loglm(~ Gender + Handedness, data = handy)
#' # loglm is just a wrapper of loglin  
#'   
#'
#'
#' 
#' 
#' 
#' 
#' 
#' 
#'   
#'
#' 
#'   
#'
#' 
#' # comparisons between hierarchical and loglin
#' ############################################################
#'
#' # the expected table given the sufficient statistics can be computed
#' # via two methods, iterative proportional fitting, and the mcmc itself:
#' out$exp # ipf
#' hierarchical(~ Gender + Handedness, data = handy, method = "mcmc")$exp
#' statsFit$fit # the equivalent in loglin; this is used by default in hierarchical
#'
#'
#'
#' 
#' # the parameter values of the loglinear model can be accessed
#' out$param
#' statsFit$param
#'
#'
#'
#'
#' # the p-value for the overall model is available as well
#' # hierarchical gives the exact conditional p-value
#' # (conditional on the sufficient statistics)
#' # the five numbers correspond the probability of tables that are
#' # "more weird" than the observed table, where "more weird" is determined
#' # by having a larger X2 value (or G2, FT, CR, or NM)
#' out$p.value
#' fisher.test(handy)$p.value # out$p.value["X2"] is accurate to monte carlo error
#' 
#' 
#' # loglin gives the p-values using the unconditional asymptotic distributions
#' c(
#'   "X2" = pchisq(statsFit$pearson, df = statsFit$df, lower.tail = FALSE),
#'   "G2" = pchisq(statsFit$lrt, df = statsFit$df, lower.tail = FALSE)
#' ) 
#'
#' out$mid.p.value # the mid (exact conditional) p-value is also available
#'
#'
#'
#'
#' # the test statistics based on the observed table and the expected
#' # table under the model are available
#' out$statistic 
#' c(statsFit$pearson, statsFit$lrt) # loglin only gives X2 and G2
#' 
#' 
#'
#'
#' # the markov basis used for the proposal distribution of the metropolis-hastings
#' # algorithm are returned. the proposal distribution is uniform on +/- 
#' # the moves added to the current table
#' out$moves
#' # they are easier understood as tables
#' vec2tab(out$moves, dim(handy))
#' # notice that the marginals stay fixed:
#' handy + vec2tab(out$moves, dim(handy))
#' 
#' 
#' 
#' 
#' # these were computed as the markov basis of the integer matrix
#' out$A
#' markov(out$A) 
#' out$moves
#' 
#' 
#' 
#' 
#' # the moves are also sometimes written in tableau form (LAS p.13)
#' tableau(out$moves, dim(handy))
#' # that's +1 the the table in elements [1,1] and [2,2]
#' # and -1 in the table in elements [1,2] and [2,1]
#' 
#'
#'
#'
#' # the acceptance probability of the MCMC is retained
#' out$acceptProb
#' 
#' 
#' 
#' 
#' # various model assessment measures are also available
#' out$quality
#' 
#' 
#' 
#'
#' # the number of independent parameters per term are in df
#' out$df
#' 
#' 
#' 
#' 
#' # as an added help, you may find the visuals in vcd useful:
#' # library(vcd)
#' # mosaic(~ Gender + Handedness, data = handy, shade = TRUE, legend = TRUE)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' ## politics example - with computing the exact p value by hand
#' ############################################################
#' 
#' data(politics)
#' 
#' (out <- hierarchical(~ Personality + Party, data = politics))
#' statsFit <- stats::loglin(politics, as.list(1:2), fit = TRUE, param = TRUE)
#' 
#' out$p.value
#' # exact without monte-carlo error
#' sum(dhyper(c(0:3,6:9), 10, 10, 9))
#' fisher.test(politics)$p.value
#' round(dhyper(0:9, 10, 10, 9), 4)
#'
#' 
#' # comparisons :
#' out$exp
#' statsFit$fit
#' 
#' out$param
#' statsFit$param
#'
#' out$p.value # exact
#' c(
#'   "X2" = pchisq(statsFit$pearson, df = statsFit$df, lower.tail = FALSE),
#'   "G2" = pchisq(statsFit$lrt, df = statsFit$df, lower.tail = FALSE)
#' ) # asymptotic approximation
#' fisher.test(politics)$p.value # accurate to monte carlo error
#'
#' out$statistic # accurate to monte carlo error
#' c(statsFit$pearson, statsFit$lrt)
#' 
#' # mosaic(~ Personality + Party, data = politics, shade = TRUE, legend = TRUE)
#' 
#' 
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' ## eyeHairColor from the Diaconis and Sturmfels reference
#' ############################################################
#' 
#' data(HairEyeColor)
#' eyeHairColor <- margin.table(HairEyeColor, 2:1)
#' 
#' outC <- hierarchical(~ Eye + Hair, data = eyeHairColor)
#' outR <- hierarchical(~ Eye + Hair, data = eyeHairColor, engine = "R")
#' 
#' # doesn't work even with workspace = 2E9 (with over 4.5Gb in memory)
#' #fisher.test(eyeHairColor, hybrid = TRUE, workspace = 2E9)
#' 
#' tableau(outC$moves, dim(eyeHairColor))
#' 
#' 
#' # library(microbenchmark)
#' # microbenchmark(
#' #   hierarchical(~ Eye + Hair, data = eyeHairColor),
#' #   hierarchical(~ Eye + Hair, data = eyeHairColor, engine = "R")
#' # )
#' # 5-10 times faster; much faster with increased iter
#' 
#' 
#' # mosaic(~ Eye + Hair, data = HairEyeColor, shade = TRUE, legend = TRUE)
#' 
#' 
#' 
#'
#'
#' 
#' ## abortion preference example from the 
#' ## Diaconis and Sturmfels reference pp. 379--381
#' ## a no 3-way interaction model
#' ############################################################
#' 
#' data(abortion)
#' 
#' out <- hierarchical(
#'   ~ Education*Abortion + Abortion*Denomination + Education*Denomination, 
#'   data = abortion,
#'   iter = 10000, burn = 50000, thin = 50
#' )
#' out$p.value
#' 
#' 
#' vec2tab(rowMeans(out$steps), dim(abortion)) # cf. p. 380
#' loglin(abortion, subsets(1:3, 2), fit = TRUE)$fit
#' 
#'
#'
#' out$param
#' loglin(abortion, subsets(1:3, 2), param = TRUE)$param
#'
#'
#'
#' qqplot(rchisq(1055, df = 8), out$sampsStats$X2s)
#' curve(1*x, from = 0, to = 30, add = TRUE, col = "red")
#' 
#' ( nMoves <- 2*ncol(out$moves) ) # DS uses 110
#' # the markov basis is larger than it needs to be
#'
#' 
#' 
#'
#'
#' 
#' 
#'
#'
#' 
#' ## loglin no three-way interaction model example
#' ############################################################
#' 
#' # the help for fits the no three-way interaction model on HairEyeColor,
#' # finds a .66196 p-value using the asymptotic distribution, and concludes
#' # a good fit:
#' data(HairEyeColor)
#' 
#' fit <- loglin(HairEyeColor, subsets(1:3, 2), fit = TRUE, param = TRUE)
#' mod <- hierarchical(~ Eye*Hair + Hair*Sex + Eye*Sex, data = HairEyeColor)
#' 
#'
#' 
#' 
#' # p values
#' pchisq(fit$lrt, fit$df, lower.tail = FALSE) # see ?loglin
#' mod$p.value
#'
#' # test statistics
#' c(fit$pearson, fit$lrt)
#' mod$statistic 
#'
#' # fits (estimated tables)
#' fit$fit
#' mod$exp
#' mod$obs
#'
#'
#' # checking the autocorrelation
#' acf(mod$sampsStats$PRs)
#' mod <- hierarchical(~ Eye*Hair + Hair*Sex + Eye*Sex, data = HairEyeColor, thin = 100)
#' acf(mod$sampsStats$PRs) # got it!
#'
#' 
#' # the slight differences in fit$fit and mod$exp (both done with ipf from loglin)
#' # are due to differences in variable order:
#' loglin(HairEyeColor, subsets(1:3, 2), fit = TRUE)$fit
#' loglin(HairEyeColor, subsets(1:3, 2)[c(1,3,2)], fit = TRUE)$fit
#'
#' # a few model moves
#' vec2tab(mod$moves[,1], dim(HairEyeColor))
#' vec2tab(mod$moves[,50], dim(HairEyeColor))
#' -vec2tab(mod$moves[,50], dim(HairEyeColor))
#' 
#' # they contribute 0 to the marginals of the table
#' vec2tab(mod$moves[,50], dim(HairEyeColor))
#' mod$A %*% mod$move[,50]
#' vec2tab(mod$A %*% mod$move[,50], dim(HairEyeColor))
#'
#' HairEyeColor 
#' HairEyeColor + vec2tab(mod$moves[,50], dim(HairEyeColor))
#' 
#'
#'
#'
#' 
#' 
#'
#'
#' ## a table with positive marginals but no MLE for
#' ## the no-three way interaction model
#' ############################################################
#'
#'
#' data(haberman)
#'
#' mod <- hierarchical(~ X1*X2 + X2*X3 + X1*X3, data = haberman)
#' 
#' statsFit <- loglin(haberman, subsets(1:3, 2), param = TRUE, fit = TRUE) 
#' statsFit$fit
#' statsFit$param
#' c(statsFit$pearson, statsFit$lrt)
#' 
#' algstatFit <- hierarchical(~ X1*X2 + X2*X3 + X1*X3, data = haberman, method = "mcmc")
#' algstatFit$exp
#' algstatFit$param
#' algstatFit$statistic
#' 
#'
#'
#'
#' 
#' 
#'
#'
#'
#' 
#' 
#'
#' ## an example from agresti, p.322
#' ############################################################
#' 
#' data(drugs) 
#' ftable(aperm(drugs, c(3, 1, 2))) # = table 8.3
#'
#' out <- hierarchical(~Alcohol + Cigarette + Marijuana, data = drugs)
#' matrix(round(aperm(out$exp, c(2,1,3)), 1), byrow = FALSE)
#' 
#' loglin(drugs, as.list(1:3), fit = TRUE)$fit
#' loglin(drugs, as.list(1:3), param = TRUE)$param
#' 
#' # # the saturated model issues a warning from markov, but works :
#' # out <- hierarchical(~Alcohol * Cigarette * Marijuana, data = drugs)
#' # matrix(round(aperm(out$exp, c(2,1,3)), 1), byrow = FALSE)
#' 
#' 
#' ftable(aperm(out$exp, c(3,1,2)))
#' 
#' stats <- loglin(drugs, as.list(1:3), fit = TRUE, param = TRUE)
#' 
#'
#' # considered via glm
#'
#' df <- as.data.frame(drugs)
#' mod <- glm(Freq ~ Alcohol + Cigarette + Marijuana, data = df, family = poisson)
#' summary(mod)
#' mod$fitted.values
#' 
#' 
#' # the same can be done with glm :
#' 
#' mod <- glm(
#'   Freq ~ Alcohol + Cigarette + Marijuana, 
#'   data = as.data.frame(drugs), family = poisson
#' )
#' summary(mod)
#' matrix(round(mod$fitted.values[c(1,3,2,4,5,7,6,8)],1))
#' 
#' 
#' 
#' mod <- glm(
#'   Freq ~ Alcohol * Cigarette + Marijuana, 
#'   data = as.data.frame(drugs), family = poisson
#' )
#' summary(mod)
#' matrix(round(mod$fitted.values[c(1,3,2,4,5,7,6,8)],1))
#' 
#' 
#' mod <- glm(
#'   Freq ~ Alcohol * Cigarette * Marijuana, 
#'   data = as.data.frame(drugs), family = poisson
#' )
#' summary(mod)
#' matrix(round(mod$fitted.values[c(1,3,2,4,5,7,6,8)],1))
#' 
#' 
#'
#'
#'
#' 
#' 
#'
#' 
#' 
#'
#'
#' 
#' }
#'
#' 
hierarchical <- function(formula, data, iter = 1E4, burn = 1000, thin = 10,
  engine = c("Cpp","R"), method = c("ipf", "mcmc"), moves){

  engine <- match.arg(engine)
  method <- match.arg(method)  

  
  ## reshape data
  ##################################################
  
  data <- suppressMessages(teshape(data, "tab"))
  varsNlevels <- dimnames(data)  
  p <- length(varsNlevels)
  
  
  ## check for sampling zeros
  ##################################################
  if(any(data == 0L)) message("care ought be taken with tables with sampling zeros to ensure the MLE exists.")


  ## parse formula for vector of r_k's
  ##################################################
  
  fString <- as.character(formula)
  predString <- fString[2]
  varsInFormula <- attr(attr(terms(formula), "factors"), "dimnames")[[1]]
  varsNndcs <- 1:length(varsInFormula)
  names(varsNndcs) <- varsInFormula
  # varsNndcs = c("Survived" = 1, ...)


  ## make list of facets
  ##################################################
  # attr(data, "dimnames")
  facets <- strsplit(predString, " \\+ ")[[1]]
  facets <- strsplit(facets, " \\* ")
  facets <- lapply(facets, function(x) unname(varsNndcs[x]) )


  ## construct A matrix and compute moves
  ##################################################  
  A <- hmat(dim(data), facets)

  if(missing(moves)){
    message("Computing moves... ", appendLF = FALSE)  	
    moves <- markov(A)
    message("done.", appendLF = TRUE)      
  }


  ## run metropolis-hastings
  ##################################################  
  current <- unname(tab2vec(data)) # init
  out <- metropolis(current, moves, 
    iter = iter, burn = burn, thin = thin,
    engine = engine)  



  ## compute data chi square
  ##################################################  
  if(method == "ipf"){
    exp <- loglin(data, facets, fit = TRUE, print = FALSE)$fit
  } else if(method == "mcmc"){
    exp <- vec2tab(rowMeans(out$steps), dim(data))
    dimnames(exp) <- dimnames(data)
  }
  e <- unname(tab2vec(exp))
  u <- t(t(unname(tab2vec(data))))
  PR <- computeUProbsCpp(matrix(u))  # unnormd prob; numers LAS 1.1.10
  X2 <- computeX2sCpp(u, e)  
  G2 <- computeG2sCpp(u, e)    
  FT <- computeCRsCpp(u, e, -.5)      
  CR <- computeCRsCpp(u, e, 2/3)
  NM <- computeNMsCpp(u, e)  

  
  ## compute MCMC chi squares
  ##################################################  
  PRs <- computeUProbsCpp(out$steps) # unnormd probs; numers LAS 1.1.10
  X2s <- computeX2sCpp(out$steps, e)  
  G2s <- computeG2sCpp(out$steps, e) 
  FTs <- computeCRsCpp(out$steps, e, -.5)
  CRs <- computeCRsCpp(out$steps, e, 2/3)  
  NMs <- computeNMsCpp(out$steps, e)   


  ## compute parameters
  ##################################################      
  # in principle, there should be one parameter for every cell.
  # there are prod(dim(data)) cells.
  # a good reference is BFH, p. 35 (and to a lesser extent 43)
  # the prod(dim(data)[terms[[j]]] - 1) line below is like
  # (I - 1) (J - 1) (K - 1)
  # CDA p.79 also helpful
  nCells <- length(data)
  dimSatModel <- nCells - 1
  degFreedom <- rep.int(0, 2^p) # there are 2^p possible subsets of vars, and
                                # therefore there are 2^p possible terms
                                
  # possibleTerms are more "types of terms" as opposed to individual terms
  # for example, an entry c(1,3) would refer to all combinations of levels
  # of variables 1 and 3; ie (# var 1 levels - 1) * (# var 3 levels - 1)
  # individual terms (parameters)
  possibleTerms <- subsets(p, include_null = TRUE)
  names(possibleTerms) <- sapply(possibleTerms, paste, collapse = " ")
  names(possibleTerms)[which(names(possibleTerms) == "")] <- "(Intercept)"    
  nVarLvls <- dim(data)
  # paramsPerTerm <- lapply(possibleTerms, function(x){
  #   if(length(x) == 0) return(1L)
  #   prod(nVarLvls[x] - 1)
  # })
  
  
  # similarly, there are the terms in the model
  termsInModel <- unique(unlist(lapply(
    lapply(facets, as.character), # to avoid subsets(2)
    subsets, include_null = TRUE), 
    recursive = FALSE
  ))
  termsInModel <- lapply(termsInModel, as.integer)
  names(termsInModel) <- sapply(termsInModel, paste, collapse = " ")  
  names(termsInModel)[which(names(termsInModel) == "")] <- "(Intercept)"
  paramsPerTermInModel <- lapply(termsInModel, function(x){
    if(length(x) == 0) return(1L) 
    prod(nVarLvls[x] - 1)
  })
  names(paramsPerTermInModel) <- unname(sapply(termsInModel, function(x){
    if(length(x) == 0) return("(Intercept)")
    paste(names(dimnames(data))[x], collapse = ".")
  }))
  nParamsInModel <- sum(unlist(paramsPerTermInModel))
  dimModel <- nParamsInModel - 1 # the - 1 accounts for the overall mean
  overallAsymptoticDegFreedom <- (dimSatModel - dimModel)
  

  # compute the parameters  
  log_fit <- exp
  log_fit[exp > 0] <- log(exp[exp > 0])  
  param <- as.list(rep(NA, length(termsInModel)))
  names(param) <- names(paramsPerTermInModel) 
  for(k in seq_along(param)){
    if(length(termsInModel[[k]]) == 0){
      param[[k]] <- mean(log_fit)
      log_fit <- log_fit - param[[k]]
    } else {
      param[[k]] <- apply(log_fit, termsInModel[[k]], mean)
      log_fit <- sweep(log_fit, termsInModel[[k]], param[[k]])
    }
  }
  # for every step, fit mle
  # then decompose mle
  # problem : they all have the same marginals, so the same
  # mles!
  # idea 1 : sample from the multinomial with the same sample
  # size (so different marginals), estimate, then decompose
  # idea 2 : bootstrap sample from the table, estimate, decompose
  # i think i like idea 2 better.
  

  # reorder the param estimates in the order of subsets
  # so you have the intercept, then all first order terms, and so on
  goodOrder <- sapply(
    c("(Intercept)", subsets(names(dimnames(data)))),
    paste, collapse = "."
  )
  param <- param[goodOrder[goodOrder %in% names(param)]]
  out$param <- param
  
  
  
  ## compute residuals and model selection, agresti p.81, 216, 324
  ##################################################  
  out$residuals <- exp
  out$residuals[exp > 0] <- 
    (data[exp > 0] - exp[exp > 0]) / sqrt(exp[exp > 0])
  
  k <- nParamsInModel  # = number of params 
  n <- sum(data)       # = sample size
  L <- dmultinom(u, sum(u), e, TRUE) # maximized log-likelihood
  BIC  <- log(n)*k - 2*L
  AIC  <-      2*k - 2*L
  AICc <- AIC + 2*k*(k+1)/(n-k-1)
  out$df <- paramsPerTermInModel
  out$quality <- c(AIC = AIC, AICc = AICc, BIC = BIC)
  

  ## add A matrix, p.value and return
  ##################################################  
  out$call <- match.call()   
  out$obs <- data  
  out$exp <- exp
  out$A <- A
  out$p.value <- c(
    PR = mean(PRs <= PR),   
    X2 = mean(X2s >= X2), 
    G2 = mean(G2s >= G2),   
    FT = mean(FTs >= FT),
    CR = mean(CRs >= CR),
    NM = mean(NMs >= NM)
  )
  out$p.value.std.err <- c(
    PR = sqrt(mean(PRs <= PR)*(1-mean(PRs <= PR))/iter), 
    X2 = sqrt(mean(X2s >= X2)*(1-mean(X2s >= X2))/iter), 
    G2 = sqrt(mean(G2s >= G2)*(1-mean(G2s >= G2))/iter),   
    FT = sqrt(mean(FTs >= FT)*(1-mean(FTs >= FT))/iter),
    CR = sqrt(mean(CRs >= CR)*(1-mean(CRs >= CR))/iter), 
    NM = sqrt(mean(NMs >= NM)*(1-mean(NMs >= NM))/iter)     
  )  
  out$mid.p.value <- c(
    PR = mean(PRs < PR) + mean(PRs == PR)/2,
    X2 = mean(X2s > X2) + mean(X2s == X2)/2, 
    G2 = mean(G2s > G2) + mean(G2s == G2)/2,
    FT = mean(FTs > FT) + mean(FTs == FT)/2,
    CR = mean(CRs > CR) + mean(CRs == CR)/2,
    NM = mean(NMs > NM) + mean(NMs == NM)/2    
  )  
  out$iter <- iter
  out$burn <- burn
  out$thin <- thin
  out$statistic = c(PR = PR, X2 = X2, G2 = G2, FT = FT, CR = CR, NM = NM)
  out$sampsStats = list(PRs = PRs, X2s = X2s, G2s = G2s, FTs = FTs, CRs = CRs, NMs = NMs)
  out$cells <- nCells
  out$method <- method

  class(out)   <- "hierarchical"
  out
}



