checklaw <- function(law.index, sample.size = 50000, law.pars = NULL, density = NULL, trunc = c(-Inf, Inf), center = FALSE, scale = FALSE) {

  if ((length(law.index) != 1) &&  !(is.numeric(law.index))) stop("law.index argument should contain a single integer value.")
  if (is.character(density)) density <- eval(parse(text = density))
  if (!is.null(density)) {if (! (class(density) == "function")) stop("density argument should contain a function")}
  
  tab <- getindex()$mat.laws  # Table of all the laws.
  index <- which(tab[, 1] == law.index)
  name <- tab[index, 2]
  nbparams <- tab[index, 3]

  out <- law.cstr(law.index, law.pars)

  if (any(names(formals(density)) %in% "pars")) {
    pars.arg <- TRUE
    density.nbpars <- length(strsplit(as.character(body(density))[2], "pars", fixed =  TRUE)[[1]]) - 1
  } else {
    pars.arg <- FALSE  
    density.nbpars <- length(formals(density)) - 1
  }

  if (!is.null(density) && (density.nbpars != out$nbparams)) {
    warning("The number of parameters (", density.nbpars, ") of the density function you have provided is different from the default number of parameters (", out$nbparams, ") of the law specified. Have you considered using dlaw", law.index, "?")
  }

  # We generate a sample of values:
  x <- gensample(law.index = law.index, n = sample.size, law.pars = out$law.pars, center = center, scale = scale)
  values <- x$sample[(x$sample >= trunc[1]) & (x$sample <= trunc[2])] # We keep only those between trunc[1] and trunc[2]
  if (center) values <- values - mean(values)
  if (scale) values <- values / sd(values)
  
  # We draw a histogram:
  hist(values, prob = TRUE, breaks = 100, xlim = range(values), main = out$name, xlab = paste("Sample size =", sample.size))

  # We add a curve of the density:
  if (!is.null(density)) {
    if (pars.arg) {
      eval(parse(text = paste("curve(density(x,pars=", paste("c(", paste(out$law.pars, collapse = ",") , ")", sep = "") , "),add=TRUE,col='blue')", sep = "")))
    } else {
      eval(parse(text = paste("curve(density(x,", paste(out$law.pars, sep = "", collapse = ","), "),add=TRUE,col='blue')", sep = "")))
    }
  }

  return(invisible(x))

}
