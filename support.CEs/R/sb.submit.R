sb.submit <- function (operation, nalternatives, nlevels, nattributes,
  design, generators, effect, interactions, determinant) {
# Arguments: see sb.design()


# Convert inputs
    nattributes.chr   <- paste(nattributes, collapse = " ")
    nlevels.chr       <- paste(nlevels, collapse = " ")
    nalternatives.chr <- paste(nalternatives)

  if (operation == "construct") {
    startingdesign.chr <- paste(apply(design, 1, paste, collapse = " "), collapse = "\n")
    choicesets.chr     <- c("")
    if (is.matrix(generators)) {
      generators.chr <- paste(apply(generators, 1, paste, collapse = " "), collapse = "\n")
    } else {
      generators.chr <- paste(generators, collapse = " ")
    }
  } else {  # operation == "check"
    choicesets.chr     <- paste(apply(design, 1, paste, collapse = " "), collapse = "\n")
    startingdesign.chr <- c("")
    generators.chr     <- c("")
  }

  if (effect == "mplussome") {
    if (is.list(interactions)) {
      interactions.chr <- sapply(interactions, paste, collapse = ",")
    } else {
      interactions.chr <- paste(interactions, collapse = ",")
    }
    interactions.chr <- paste(interactions.chr, collapse = " ")
  } else {
    interactions.chr <- c("")
  }

  if (!is.null(determinant)) {
    determinant.chr <- determinant
  } else {
    determinant.chr <- c("")
  }


# Extract software version
  process.page <- htmlTreeParse("http://130.56.248.113/choice/")
  version      <- xmlValue(xmlChildren(process.page$children$html)$body[9][[1]])


# Request to construct/check choice sets
  rtn <- postForm(
    uri     = c("http://130.56.248.113/choice/process/"),
    chsets  = choicesets.chr, 
    corc    = operation, 
    det     = determinant.chr,
    effect  = effect, 
    factors = nattributes.chr, 
    levels  = nlevels.chr, 
    gens    = generators.chr,
    msize   = nalternatives.chr, 
    tmts    = startingdesign.chr, 
    twofis  = interactions.chr,
    style   = "post")


# Return
  return(rtn <- list(html = rtn, version = version))
}

