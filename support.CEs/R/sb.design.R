sb.design <- function(operation = "construct", 
  nattributes, nalternatives, nlevels, attribute.names,
  design = NULL, generators = NULL, effect = "main",
  interactions = NULL, determinant = NULL,
  nblocks = 1, seed = NULL, ...) {

# Generate an OMED using oa.design()
  if (operation == "construct" & is.null(design)) {
    # Generate an OMED
    design <- data.matrix(oa.design(nlevels = nlevels, seed = seed, ...)) - 1

    # Validate nblocks and number of rows of design generated
    if (!isTRUE(all.equal(nrow(design) %% nblocks, 0))) {
      cat("Number of rows of design generated using oa.design() is ", nrow(design), "\n")
      cat("'nblocks' is", nblocks, "\n") 
      stop("'nblocks' must be divisors of number of rows of design")
    }
  }


# Validate arguments
  sb.check.args(operation       = operation,
               nalternatives   = nalternatives, 
               nlevels         = nlevels,
               nattributes     = nattributes,
               attribute.names = attribute.names,
               design          = design, 
               generators      = generators,
               effect          = effect,
               interactions    = interactions,
               determinant     = determinant,
               nblocks         = nblocks,
               seed            = seed)


  if (!is.null(design) & !is.matrix(design)) {
    design <- data.matrix(design)
  }
  storage.mode(design) <- "integer"


# Submit HTML form to the WEB site "Discrete Choice Experiments"
  html.content <- sb.submit(operation     = operation,
                            nalternatives = nalternatives, 
                            nlevels       = nlevels,
                            nattributes   = nattributes,
                            design        = design,
                            generators    = generators,
                            effect        = effect,
                            interactions  = interactions,
                            determinant   = determinant)


# Extract and restructuring output from the WEB site  
  sb.output <- sb.extract(x         = html.content,
                          effect    = effect, 
                          operation = operation)


# Return when construction/check on the WEB site is stopped
  if (sb.output$calculation == FALSE) {
    print(sb.output$messages)
    cat("See html content returned for detail.\n")
    return(sb.output)
  }


# Format output
  if (operation == "construct") {
    # choice sets
    ALTS <- vector("list", nalternatives)
    choicesets <- sb.output$output$choice.sets
    nquestions <- nrow(choicesets)
    nquestions_nblocks <- nquestions / nblocks
    block <- rep(1:nblocks, each = nquestions_nblocks)
    if (nblocks >= 2) { 
      set.seed(seed)
      block.rnd <- sample(block, nquestions, replace = FALSE)
      choicesets <- cbind(choicesets, block.rnd)
      choicesets <- choicesets[order(choicesets[, "block.rnd"]), ]
    }
    ques  <- rep(1:nquestions_nblocks, nblocks)
    for (i in 1:nalternatives) {
      ALTS[[i]] <- data.frame(cbind(block,
                                    ques,
                                    rep(i, nquestions),
                                    choicesets[, c(c(1:nattributes) + nattributes * (i - 1))]))
      colnames(ALTS[[i]]) <- c("BLOCK", "QES", "ALT", names(attribute.names))
      for (j in 1:nattributes) {
        ALTS[[i]][, (j + 3)] <- factor(ALTS[[i]][, (j + 3)])
        levels(ALTS[[i]][, (j + 3)]) <- attribute.names[[j]]
      }
    }
    names(ALTS) <- paste("alt.", 1:nalternatives, sep = "")

    # design information
    desinf <- list(nalternatives = nalternatives,
                   nblocks       = nblocks, 
                   nquestions    = nquestions_nblocks,
                   nattributes   = nattributes)
    # output
    rtn <- list(alternatives       = ALTS, 
                candidate          = sb.output$input$treatment,
                design.information = desinf,
                sb                 = sb.output)
  } else {
    # output
    rtn <- list(alternatives       = NULL,
                candidate          = NULL,
                design.information = NULL,
                sb                 = sb.output)
  }
  class(rtn) <- c("sb", "cedes")


# Return output
  return(rtn)
}