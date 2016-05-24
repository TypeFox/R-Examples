wtd.boxplot.formula <-
function(formula, weights=NULL, data = NULL, ..., subset, na.action = NULL)
{
    if(missing(formula) || (length(formula) != 3))
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m$... <- NULL
    m$na.action <- na.action # force use of default for this method
    m[[1]] <- as.name("model.frame")
    if(!is.null(weights)) m$weights<-NULL
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")
    myfactors<-mf[-response]
    if (is.null(weights)){
    wtd.boxplot(split(mf[[response]], myfactors),
         weights=NULL, ...)
    }
    else  {
    if (missing(subset))
    wtd.boxplot(split(mf[[response]], myfactors),
        weights=split(weights,myfactors),
     ...)
    else {
    if (is.null(subset))
    wtd.boxplot(split(mf[[response]], myfactors),
        weights=split(weights,myfactors),
     ...)
    else
    wtd.boxplot(split(mf[[response]], myfactors),
        weights=split(weights[subset],myfactors),
     ...)
      }
    }
}

