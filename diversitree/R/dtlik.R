get.info <- function(x)
  get.cache(x)$info
get.cache <- function(x) {
  if ( inherits(x, "big.brother") || is.constrained(x) )
    get.cache(attr(x, "func"))
  else
    environment(x)$cache
}
set.info <- function(x, value)
  environment(x)$cache$info <- value

## TODO: The formatting here could be abstracted (using a list of
## lists to represent levels of bulleting)
print.dtlik <- function(x, with.args=TRUE, with.phy=TRUE,
                        with.ref=TRUE, ...) {
  info <- get.info(x)
  argnames <- argnames(x)
  str <- 
    c(sprintf("%s likelihood function:", info$name.pretty),
      sprintf("  * Parameter vector takes %d elements:",
              length(argnames)),
      strwrap(paste(argnames, collapse=", "), initial="     - ",
              exdent=7))
  if ( is.constrained(x) )
    str <- c(str, .print.constraints(x))
  if ( with.args )
    str <- c(str, .print.func.args(x))
  if ( !is.null(info$phy) && with.phy )
    str <- c(str, .print.func.phy(info))
  if ( with.ref )
    str <- c(str, .print.func.reference(info))
  cat(paste(str, collapse="\n"), "\n", sep="")
  cat("R definition:\n")
  str(args(x))
}

.print.constraints <- function(f) {
  info <- get.info(f)
  str <-
    sprintf("  * Function constrained (original took %d elements):",
            length(info$argnames))
  if ( inherits(f, "constrained.i") ) {
    e <- environment(f)
    if ( !is.null(e$argnames) )
      str <- c(str,
               sprintf("     - %s ~ %s",
                       e$argnames[-e$i.free],
                       prettyNum(e$p[-e$i.free])))
  } else {
    rels <- environment(f)$rels
    str <- c(str,
             sprintf("     - %s ~ %s", names(rels),
                     unlist(lapply(rels, deparse))))
  }
  str
}
  
.print.func.args <- function(f) {
  args <- formals(f)
  common <-
    list(pars="Parameter vector",
         condition.surv="Condition likelihood on survial?",
         root="Type of root treatment",
         root.p="Vector of root state probabilities",
         intermediates="Also return intermediate values?",
         pars.only="Return full parameter vector?")
  if ( is.constrained(f) )
    common[['...']] <- "Additional arguments to underlying function"
  
  ## It is possible that when an argument is truly "" as default, this
  ## will miss it.
  f1 <- function(v) {
    def <- deparse(args[[v]])
    str <- if (def=="") v else sprintf("%s [%s]", v, def)
    if ( v %in% names(common) )
      str <- sprintf("%s: %s", str, common[[v]])
    str
  }
  c(sprintf("  * Function takes arguments (with defaults)"),
    strwrap(sapply(names(args), f1), prefix="     - ", exdent=12))
}

.print.func.phy <- function(info) {
  phy <- info$phy
  str <- sprintf("  * Phylogeny with %s tips and %s nodes",
                 length(phy$tip.label), phy$Nnode)
  spp <- paste(phy$tip.label, collapse=", ")
  spp <- strwrap(spp, initial="     - Taxa: ", .9*getOption("width")-5)
  if ( length(spp) > 1 )
    spp <- sprintf("%s ...", spp[1])
  c(str, spp)
}

.print.func.reference <- function(info) {
  reference <- info$reference
  if ( length(reference) == 0 )
    return(character(0))
  str <- if ( length(reference) == 1 )
    "  * Reference:" else "  * References:"
  c(str, sprintf("     - %s", reference))
}

argnames.dtlik <- function(x, ...) {
  info <- get.info(x)
  info$argnames
}
`argnames<-.dtlik` <- function(x, value) {
  ## Here, np will always be the *base* np, not the actual number of
  ## parameters.  Care will be needed when processing...
  info <- get.info(x)
  if ( is.character(value) ) {
    if ( length(value) != info$np )
      stop(sprintf("Wrong argnames length: expected %d, got %d",
                   info$np, length(value)))
    if ( any(is.na(value)) )
      stop("NA values in argnames replacement")
    if ( any(duplicated(value)) )
      stop("Duplicated names in argnames replacement")
    set.argnames(x, value)
    x
  } else {
    stop("Still to write...")
  }
}

set.argnames <- function(x, value)
  environment(x)$cache$info$argnames <- value

find.mle.dtlik <- function(func, x.init, method,
                           fail.value=NA, ...) {
  info <- get.info(func)
  name <- info$name
  if ( is.constrained(func) )
    x.init <- guess.constrained.start(func, x.init)
  if ( missing(method) ) method <- info$ml.default
  ans <- do.mle.search(func, x.init, method, fail.value, ...)
  class(ans) <- c(sprintf("fit.mle.%s", name), class(ans))
  ans
}
