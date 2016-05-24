#######################################
### doubledecker plot

doubledecker <- function(x, ...)
  UseMethod("doubledecker")

doubledecker.formula <-
function(formula, data = NULL, ..., main = NULL)
{
    if (is.logical(main) && main)
      main <- deparse(substitute(data))

    if (is.structable(data))
      data <- as.table(data)

    m <- match.call(expand.dots = FALSE)
    edata <- eval(m$data, parent.frame())

    fstr <- strsplit(paste(deparse(formula), collapse = ""), "~")
    vars <- strsplit(strsplit(gsub(" ", "", fstr[[1]][2]), "\\|")[[1]], "\\+")
    dep <- gsub(" ", "", fstr[[1]][1])
    varnames <- vars[[1]]
    if (dep == "")
      stop("Need a dependent variable!")
    varnames <- c(varnames, dep)

    if(inherits(edata, "ftable")
       || inherits(edata, "table")
       || length(dim(edata)) > 2) {
        dat <- as.table(data)
        if(all(varnames != ".")) {

          ind <- match(varnames, names(dimnames(dat)))
          if (any(is.na(ind)))
            stop(paste("Can't find", paste(varnames[is.na(ind)], collapse=" / "), "in", deparse(substitute(data))))

          dat <- margin.table(dat, ind)
        } else {
          ind <- match(dep, names(dimnames(dat)))
          if (is.na(ind))
            stop(paste("Can't find", dep, "in", deparse(substitute(data))))
          dat <- aperm(dat, c(seq_along(dim(dat))[-ind], ind))
        }
        doubledecker.default(dat, main = main, ...)
      } else {
        tab <- if ("Freq" %in% colnames(data))
          xtabs(formula(paste("Freq~", varnames, collapse = "+")),
                data = data)
        else
          xtabs(formula(paste("~", varnames, collapse = "+")),
                data = data)

        doubledecker.default(tab, main = main, ...)
      }
  }

doubledecker.default <- function(x,
                         depvar = length(dim(x)),
                         margins = c(1, 4, length(dim(x)) + 1, 1),
                         gp = gpar(fill = rev(gray.colors(tail(dim(x), 1)))),
                         labeling = labeling_doubledecker,
                         spacing = spacing_highlighting,
                         main = NULL,
                         keep_aspect_ratio = FALSE,
                         ...) {
  x <- as.table(x)
  d <- dim(x)
  l <- length(d)
  if (is.character(depvar))
    depvar <- match(depvar, names(dimnames(x)))
  condvars <- (1:l)[-depvar]
  x <- aperm(x, c(condvars, depvar))
  strucplot(x,
            core = struc_mosaic(zero_split = FALSE, zero_shade = FALSE),
            condvars = l - 1,
            spacing = spacing,
            split_vertical = c(rep.int(TRUE, l - 1), FALSE),
            gp = gp,
            shade = TRUE,
            labeling = labeling,
            main = main,
            margins = margins,
            legend = NULL,
            keep_aspect_ratio = keep_aspect_ratio,
            ...
            )
}
