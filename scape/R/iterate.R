iterate <- function(model, ceiling=Inf, p=1, digits.n=0, digits.sigma=2)
{
  ## 1  Define functions
  summaryN <- function(what, series=NULL)
  {
    data.frame(SS=getN(model, what=what, series=series, digits=digits.n),
               nhat=estN(model, what=what, series=series, init=FALSE, ceiling=ceiling, digits=digits.n),
               candbar=estN(model, what=what, series=series, ceiling=ceiling, digits=digits.n),
               candmed=estN(model, what=what, series=series, FUN=median, ceiling=ceiling, digits=digits.n),
               candbar1=estN(model, what=what, series=series, init=1, ceiling=ceiling, digits=digits.n),
               candmed1=estN(model, what=what, series=series, init=1, FUN=median, ceiling=ceiling, digits=digits.n))
  }
  summarySigmaI <- function(what, series=NULL)
  {
    data.frame(sigma=getSigmaI(model, what=what, series=series, digits=digits.sigma),
               sigmahat=estSigmaI(model, what=what, series=series, init=1, p=p, digits=digits.sigma),
               candbar=estSigmaI(model, what=what, series=series, p=p, digits=digits.sigma),
               candmed=estSigmaI(model, what=what, series=series, FUN=median, p=p, digits=digits.sigma))
  }

  output <- list()

  ## 2  Dev
  if(!is.null(model$Dev))
    output$Dev <- data.frame(sigmaR=getSigmaR(model,digits=digits.sigma),
                             sigmahat=estSigmaR(model,digits=digits.sigma))

  ## 3  CPUE
  if(!is.null(model$CPUE))
  {
    series <- unique(model$CPUE$Series)
    if(length(series) == 1)
    {
      output$CPUE <- summarySigmaI("c")
    }
    else
    {
      output$CPUE <- lapply(series, summarySigmaI, what="c")
      names(output$CPUE) <- series
    }
  }

  ## 4  Survey
  if(!is.null(model$Survey))
  {
    series <- unique(model$Survey$Series)
    if(length(series) == 1)
    {
      output$Survey <- summarySigmaI("s")
    }
    else
    {
      output$Survey <- lapply(series, summarySigmaI, what="s")
      names(output$Survey) <- series
    }
  }

  ## 5  CAc
  if(!is.null(model$CAc))
  {
    series <- unique(model$CAc$Series)
    if(length(series) == 1)
    {
      output$CAc <- summaryN("CAc")
    }
    else
    {
      output$CAc <- lapply(series, summaryN, what="CAc")
      names(output$CAc) <- series
    }
  }

  ## 6  CAs
  if(!is.null(model$CAs))
  {
    series <- unique(model$CAs$Series)
    if(length(series) == 1)
    {
      output$CAs <- summaryN("CAs")
    }
    else
    {
      output$CAs <- lapply(series, summaryN, what="CAs")
      names(output$CAs) <- series
    }
  }

  ## 7  CLc
  if(!is.null(model$CLc))
  {
    series <- unique(model$CLc$Series)
    if(length(series) == 1)
    {
      output$CLc <- summaryN("CLc")
    }
    else
    {
      output$CLc <- lapply(series, summaryN, what="CLc")
      names(output$CLc) <- series
    }
  }

  ## 8  CLs
  if(!is.null(model$CLs))
  {
    series <- unique(model$CLs$Series)
    if(length(series) == 1)
    {
      output$CLs <- summaryN("CLs")
    }
    else
    {
      output$CLs <- lapply(series, summaryN, what="CLs")
      names(output$CLs) <- series
    }
  }

  return(output)
}
