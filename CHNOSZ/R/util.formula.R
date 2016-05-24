# CHNOSZ/util.formula.R
# functions to compute some properties of chemical formulas

get.formula <- function(formula) {
  # return the argument if it's a matrix
  if(is.matrix(formula)) return(formula)
  # return the argument as matrix if it's a data frame
  if(is.data.frame(formula)) return(as.matrix(formula))
  # return the values in the argument, or chemical formula(s) 
  # for values that are species indices
  # for numeric values, get the formulas from those rownumbers of thermo$obigt
  i <- suppressWarnings(as.numeric(formula))
  # we can't have more than the number of rows in thermo$obigt
  thermo <- get("thermo")
  iover <- i > nrow(thermo$obigt)
  iover[is.na(iover)] <- FALSE
  if(any(iover)) stop(paste("species number(s)",paste(i[iover],collapse=" "),
    "not available in thermo$obigt"))
  # we let negative numbers pass as formulas
  i[i < 0] <- NA
  # replace any species indices with formulas from thermo$obigt
  formula[!is.na(i)] <- thermo$obigt$formula[i[!is.na(i)]]
  return(formula)
}

i2A <- function(formula) {
  ## assemble the stoichiometric matrix (A)
  ## for the given formulas  20120108 jmd
  if(is.matrix(formula)) {
    # do nothing if the argument is already a matrix
    A <- formula
  } else {
    # get the elemental makeup of each formula, counting
    # zero for elements that appear only in other formulas
    msz <- makeup(formula, count.zero=TRUE)
    # convert formulas into a stoichiometric matrix with elements on the columns
    A <- t(sapply(msz, c))
    # add names from character argument
    # or from thermo$obigt for numeric argument
    if(is.numeric(formula[1])) rownames(A) <- get("thermo")$obigt$name[formula]
    else rownames(A) <- formula
  }
  return(A)
}

as.chemical.formula <- function(makeup, drop.zero=TRUE) {
  # make a formula character string from the output of makeup()
  # or from a stoichiometric matrix (output of i2A() or protein.formula())
  # first define a function to work with a single makeup object
  cffun <- function(makeup) {
    # first strip zeroes if needed
    if(drop.zero) makeup <- makeup[makeup!=0]
    # Z always goes at end
    makeup <- c(makeup[names(makeup)!="Z"], makeup[names(makeup)=="Z"])
    # the elements and coefficients
    elements <- names(makeup)
    coefficients <- as.character(makeup)
    # any 1's get zapped
    coefficients[makeup==1] <- ""
    # any Z's get zapped (if they're followed by a negative number)
    # or turned into a plus sign (to indicate a positive charge)
    elements[elements=="Z" & makeup < 0] <- ""
    elements[elements=="Z" & makeup >= 0] <- "+"
    # put the elements and coefficients together
    formula <- paste(elements, coefficients, sep="", collapse="")
    # if the formula is uncharged, and the last element has a negative
    # coefficient, add an explicit +0 at the end
    if(!"Z" %in% names(makeup) & tail(makeup,1) < 0) 
      formula <- paste(formula, "+0", sep="")
    return(formula)
  }
  # call cffun() once for a single makeup, or loop for a matrix
  if(is.matrix(makeup)) out <- sapply(1:nrow(makeup), function(i) {
    mkp <- makeup[i, ]
    return(cffun(mkp))
  }) else out <- cffun(makeup)
  return(out)
}

mass <- function(formula) {
  # calculate the mass of elements in chemical formulas
  thermo <- get("thermo")
  formula <- i2A(get.formula(formula))
  ielem <- match(colnames(formula), thermo$element$element)
  if(any(is.na(ielem))) stop(paste("element(s)",
    colnames(formula)[is.na(ielem)], "not available in thermo$element"))
  mass <- as.numeric(formula %*% thermo$element$mass[ielem])
  return(mass)
}

entropy <- function(formula) {
  # calculate the standard molal entropy at Tref of elements in chemical formulas
  thermo <- get("thermo")
  formula <- i2A(get.formula(formula))
  ielem <- match(colnames(formula), thermo$element$element)
  if(any(is.na(ielem))) warning(paste("element(s)",
    paste(colnames(formula)[is.na(ielem)], collapse=" "), "not available in thermo$element"))
  entropy <- as.numeric( formula %*% (thermo$element$s[ielem] / thermo$element$n[ielem]) )
  return(entropy)
}

GHS <- function(formula, G=NA, H=NA, S=NA, T=get("thermo")$opt$Tr) {
  # for all NA in G, H and S, do nothing
  # for no  NA in G, H and S, do nothing
  # for one NA in G, H and S, calculate its value from the other two:
  # G - standard molal Gibbs energy of formation from the elements
  # H - standard molal enthalpy of formation from the elements
  # S - standard molal entropy
  # argument checking
  if(!all(diff(sapply(list(formula, G, H, S), length)) == 0))
    stop("formula, G, H and S arguments are not same length")
  # calculate Se (entropy of elements)
  Se <- entropy(formula)
  # calculate one of G, H, or S if the other two are given
  GHS <- lapply(seq_along(G), function(i) {
    G <- G[i]
    H <- H[i]
    S <- S[i]
    Se <- Se[i]
    if(is.na(G)) G <- H - T * (S - Se)
    else if(is.na(H)) H <- G + T * (S - Se)
    else if(is.na(S)) S <- (H - G) / T + Se
    return(list(G, H, S))
  })
  # turn the list into a matrix and add labels
  GHS <- t(sapply(GHS, c))
  colnames(GHS) <- c("G", "H", "S")
  rownames(GHS) <- formula
  return(GHS)
}

ZC <- function(formula) {
  # calculate average oxidation state of carbon 
  # from chemical formulas of species
  # if we haven't been supplied with a stoichiometric matrix, first get the formulas
  formula <- i2A(get.formula(formula))
  # is there carbon there?
  iC <- match("C", colnames(formula))
  if(is.na(iC)) stop("carbon not found in the stoichiometric matrix")
  # the nominal charges of elements other than carbon
  # FIXME: add more elements, warn about missing ones
  knownelement <- c("H", "N", "O", "S", "Z")
  charge <- c(-1, 3, 2, 2, 1)
  # where are these elements in the formulas?
  iknown <- match(knownelement, colnames(formula))
  # any unknown elements in formula get dropped with a warning
  iunk <- !colnames(formula) %in% c(knownelement, "C")
  if(any(iunk)) warning(paste("element(s)",paste(colnames(formula)[iunk], collapse=" "),
    "not in", paste(knownelement, collapse=" "), "so not included in this calculation"))
  # contribution to charge only from known elements that are in the formula
  formulacharges <- t(formula[,iknown[!is.na(iknown)]]) * charge[!is.na(iknown)]
  # sum of the charges; the arrangement depends on the number of formulas
  if(nrow(formula)==1) formulacharge <- rowSums(formulacharges) 
  else formulacharge <- colSums(formulacharges)
  # numbers of carbons
  nC <- formula[,iC]
  # average oxidation state
  ZC <- as.numeric(formulacharge/nC)
  return(ZC)
}
