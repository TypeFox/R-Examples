#' @include terms.R
{}

#' Generate design and model information for \code{spikeSlabGAM}
#'
#' This function generates the design matrix for a generalized additive (mixed)
#' model, based on smoothing spline ANOVA-like orthogonal decomposition of the
#' model terms and their interactions. It parses the formula given to
#' \code{\link{spikeSlabGAM}} to provide all the arguments necessary for the
#' MCMC sampler.
#'
#' Setting \code{lowRankInteractions} to FALSE can result in very large models,
#' especially if higher-order interactions or interactions between terms with
#' lots of parameters are involved. Note that numeric covariates with fewer
#' unique values than \code{minUniqueValues} are treated as factors unless
#' wrapped in a special argument.
#'
#' This function is not meant to be called directly by the user,
#' \code{\link{spikeSlabGAM}} is the user interface.
#' @param formula the model formula. Follows the usual R syntax described in
#'   \code{\link[stats]{formula}}. Terms in the formula that are not in the list
#'   of specials are automatically assigned an appropriate special, i.e.
#'   numerical covariates \code{x} are treated as \code{lin(x) + sm(x)}, factors
#'   \code{f} are treated as \code{fct(f)} (see \code{specials} argument). See
#'   \code{\link{spikeSlabGAM}} for more details.
#' @param data \code{data.frame} containing all the variables in the function
#' @param specials a vector of the types of possible model terms. These must be
#'   implemented as functions generating a design matrix with a label attribute.
#'   See also \code{\link{sm}()} for an example. The documentation for
#'   \code{\link{spikeSlabGAM}} contains a list of implemented model term types
#'   and usage examples.
#' @param minUniqueValues the minimal number of unique values a covariate has to
#'   have in order to not be treated as a factor. Defaults to 6.
#' @param lowRankInteractions should a low-rank approximation of the design
#'   matrix for interaction terms based on a (truncated) spectral decomposition
#'   of the implied covariance matrix be used? defaults to TRUE.
#' @param orthogonalizeInteractions should the design matrices for interaction
#'   terms be projected into the complement of the column space of the
#'   respective main effects?  Can help separate marginal and interaction
#'   effects. Defaults to TRUE.
#' @param decomposition which decomposition to use, see \code{\link{sm}}.
#'   Defaults to the default of \code{\link{sm}}.
#' @return a list with components: \item{response}{the left hand side of the
#'   model formula} \item{Design}{the design matrix of the model specified in
#'   formula} \item{groupIndicators}{a factor that maps the columns of
#'   \code{Design} to the different model terms} \item{H}{a matrix containing
#'   the hierarchy of the penalization groups (not used, included for backwards
#'   compatibility)}
#' @export
#' @author Fabian Scheipl
ssGAMDesign <-  function(formula, data,
  specials = c('u', 'lin', 'fct', 'sm', 'rnd', 'mrf', 'srf'),
  minUniqueValues = 6,
  lowRankInteractions = TRUE,
  orthogonalizeInteractions = TRUE,
  decomposition = NULL) {

  replaceRawTerms <- function(formula, specials, decomposition) {
    # replace all non-specials with specials
    # (i.e. numeric:  x -> lin(x) + sm(x) (or lin(x) if GAM), factor x -> fct(x))
    # and rule out lin(x):sm(x) interactions in higher order terms

    rawTrms <- terms.formula(formula, specials)
    # find out which terms do not involve special terms
    rawLabels <- if(length(unlist(attr(rawTrms, "specials")))) {
      whereSpecials <- unique(unlist(apply(
        attr(rawTrms, "factors")[unlist(attr(rawTrms, "specials")), , drop = F] > 0,
        1, which)))
      if(length(whereSpecials)) {
        attr(rawTrms, "term.labels")[-whereSpecials]
      } else {
        attr(rawTrms, "term.labels")
      }
    } else {
      attr(rawTrms, "term.labels")
    }
    # split up interactions:
    rawLabels <- unique(unlist(sapply(rawLabels, strsplit, split = c(":","*"))))

    replace1Term <- function(x) {
      decompArg <- ifelse(!is.null(decomposition),
        paste(", decomposition ='", decomposition,"'", sep =""), "")
      ifelse(is.numeric(data[, x]),
        ifelse(length(unique(data[, x])) > minUniqueValues,
          paste("(lin(", x,") + sm(", x, decompArg,"))", sep =""),
          paste("fct(", x,")", sep ="")),
        paste("fct(", x,")", sep =""))
    }

    rawLabelReplacements <-  sapply(rawLabels, replace1Term)

    replaceRawLabels <- function(x) {
      #step through parse tree unless at terminal node (which are always 'name's):
      if (!is.name(x)) {
        #leave special terms (which are function 'call's, not 'name's) unchanged
        if(as.character(x[[1]]) %in% specials) {
          return(x)
        } else{
          #go further "left"
          if(length(x)>1) x[[2]] <- replaceRawLabels(x[[2]])
          #go further "right"
          if(length(x)== 3) x[[3]] <- replaceRawLabels(x[[3]])
        }
      } else {
        #find replacement and change entry
        which <- grep(paste("^", as.character(x),"$", sep =""),
          names(rawLabelReplacements))
        if(length(which)) x <- as.name(rawLabelReplacements[which])
      }
      return(x)
    }
    newFormulaString <- {
      tmp <- paste(deparse(replaceRawLabels(formula)), collapse ="")
      gsub("  ","", gsub("`","", tmp))
    }

    rmNonsenseInteractions <- function(newFormulaString)
    {
      # rm all interactions that contain interactions between linear and
      # smooth effects of the same covariate
      trms <- terms(formula(newFormulaString), specials)
      mainEfs <- attr(trms, "term.labels")[attr(trms, "order")== 1]
      intActs <- attr(trms, "term.labels")[attr(trms, "order")>1]
      tmp <- c("")
      if(length(intActs)) {
        # for each covariate with a linear term...
        for(l in grep("lin(", mainEfs, fixed = T)) {
          x <- substr(mainEfs[l], 5, nchar(mainEfs[l])-1)
          # ...remove interactions in which the covariate occurs twice
          if(length(x)) {
            nonsense <- sapply(gregexpr(x, intActs, fixed = T), length) > 1
            if(any(nonsense)) tmp <- c(tmp, paste("-", intActs[nonsense],
              collapse =""))
          }
        }
      }
      tmp
    }
    noNonsense <- rmNonsenseInteractions(newFormulaString)

    # TODO: make all interactions with u() into u()-terms
    newFormulaString <- paste(newFormulaString, paste(noNonsense, collapse =""))
    return(formula(newFormulaString))
  }

  formula <- replaceRawTerms(formula, specials, decomposition)

  trms <- terms.formula(formula, specials = specials )
  labels <- attr(trms,"term.labels")

  # make designs for all main effects
  mainEffects <- trms[attr(trms, "order")== 1]
  mainEffectDesigns <- sapply(1:sum(attr(trms, "order")== 1), function(i) {
    with(data, eval(mainEffects[i][[3]]))
  }, simplify = FALSE)
  names(mainEffectDesigns) <- labels[attr(trms, "order")== 1]

  # make designs for all interactions
  interactions <- labels[attr(trms, "order")!= 1]

  makeInteractions <- function(term, Designs,
    center = orthogonalizeInteractions,
    orthogonalize = lowRankInteractions) {
    terms <- strsplit(term, ":")[[1]]
    Xs <- Designs[terms]
    label <- paste(sapply(Xs, attr, which ="label"), collapse =":")
    # design for interaction contains column-wise product
    # of main effect designs:
    B <- matrix(apply(Xs[[1]], 2, function(x) {
      x * Xs[[2]]
    }),
      nrow = nrow(Xs[[1]]), ncol = NCOL(Xs[[1]])* NCOL(Xs[[2]]))
    colnames(B) <- namesB <- paste(rep(colnames(Xs[[1]]), e = NCOL(Xs[[2]])),
      rep(colnames(Xs[[2]]), t = NCOL(Xs[[1]])), sep =":")
    if(length(Xs)>2) {
      for(i in 3:length(Xs)) {
        B <- matrix(apply(B, 2, function(x) {
          x * Xs[[i]]
        }),
          nrow = nrow(B), ncol = NCOL(B)* NCOL(Xs[[i]]))
        colnames(B) <- namesB <- paste(rep(namesB, e = NCOL(Xs[[i]])),
          rep(colnames(Xs[[i]]), t = length(namesB)), sep =":")
      }
    }
    qrB <- qr(B)
    interActRank <- qrB$rank

    if(center) {
      proj <- qr(cbind(1, do.call(cbind, Xs)))
      B <- qr.resid(proj, B)
      interActRank <- interActRank - proj$rank
    }
    if(orthogonalize && NCOL(B)>1) {
      # use truncated svd of B instead of eigen/svd of BB' as basis
      # i.e. orthoDesign(C = tcrossprod(B), rank = interActRank)
      B <- {
        #irlba works well only for larger matrices
        if(ncol(B)>20 & interActRank > 10) {
          eC <- try(irlba(B, max(interActRank, sapply(Xs, NCOL)),
            max(interActRank, sapply(Xs, NCOL)), adjust = 0))
        } else {
          sv <- svd(B, nv = 0)
          eC <- list(vectors = sv$u, values = sv$d)
        }
        if(class(eC)=="try-error") {
          sv <- svd(B, nv = 0)
          eC <- list(vectors = sv$u, values = sv$d)
        }
        eC$values <- eC$values^2

        nullvals <- eC$values < 10e-10
        colsZ <- max(3,
          min(which( cumsum(eC$values[!nullvals])/sum(eC$values[!nullvals]) >
              .999)))

        colsZ <- min(colsZ, sum(!nullvals))
        use <- 1:colsZ
        t(eC$values[use]* t(eC$vectors[, use]))
      }
    }
    colnames(B) <- paste(label,".", 1:NCOL(B), sep ="")

    B <- scaleMat(B)

    return(structure(B, label = label, predvars = list(qrB = qrB)))
  }

  interactionDesigns <- sapply(interactions, makeInteractions,
    Designs = mainEffectDesigns, simplify = FALSE)

  #initialize Design
  y  <- with(data, eval(attr(trms, "variables")))[[attr(trms,"response")]]
  n <- length(y)
  if(attr(trms,"intercept")== 1) {
    X <-matrix(1, nrow = n, ncol = 1)
    colnames(X) <- "u(Int)"
    grpInds <- c("u")
  }else {
    X <- matrix(NA, nrow = n, ncol = 0)
    grpInds <- c()
  }

  X <- cbind(X, do.call(cbind, mainEffectDesigns),
    do.call(cbind, interactionDesigns))
  grpInds <- unlist(c(grpInds,
    rep(sapply(mainEffectDesigns, attr, which ="label"),
      sapply(mainEffectDesigns, NCOL)),
    rep(sapply(interactionDesigns, attr, which ="label"),
      sapply(interactionDesigns, NCOL))))
  names(grpInds) <- grpInds
  grpInds[grep("u(", grpInds, fixed = T)] <- "u"

  makeHierarchy <- function(formula, grpInds) {
    trms <- terms(formula)
    H <- attr(trms, "factors")
    H <- H[-1,, drop = FALSE]
    order <- attr(trms, "order")
    labels <- c(sapply(mainEffectDesigns, attr, which ="label"),
      sapply(interactionDesigns, attr, which ="label"))
    colnames(H) <- labels
    if(any(order > 1)) {
      for(ord in 2:max(order)) {
        newH <- matrix(0, nrow = sum(order == ord), ncol(H))
        H <- rbind(H, newH)
      }
    }
    H <- (H + t(H))
    diag(H) <- 1
    rownames(H) <- labels
    #TODO: hierarchy for higher order interactions (i.e. which 2-interactions
    #are involved in which 3-interaction missing)
    return(H)
  }
  H <- makeHierarchy(formula, grpInds)
  predvars <- c(lapply(mainEffectDesigns, attr, "predvars"),
    lapply(interactionDesigns, attr, "predvars"))
  names(predvars) <- c(sapply(mainEffectDesigns, attr, which ="label"),
    sapply(interactionDesigns, attr, which ="label"))

  return(list(response = y, Design = X, groupIndicators = grpInds, H = H,
    formula = formula, predvars = predvars))
}


