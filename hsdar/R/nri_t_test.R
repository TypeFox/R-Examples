setMethod("t.test", signature(x = "Nri"),
          function(x, ...)
{
  ttest_apply <- function(indices)
  {
    x <- get("x", envir= environment_apply)@nri[indices[1],indices[2],]
    y <- get("y", envir= environment_apply)@nri[indices[1],indices[2],]
    t_res <- eval(parse(text = get("eval_text", envir= environment_apply)))
    return(c(t_res$statistic, t_res$p.value))
  }
  
  stopifnot(class(x) == "Nri")
  
  dots <- list(...)
  goty <- 0
  if (length(dots) > 0)
  {
    if (any(names(dots)=="y"))
    {
      goty <- which(names(dots)=="y")
      goty <- goty[length(goty)]
    } else {
      i = length(dots) + 1
      while (i > 1 & goty == 0)
      {
        i <- i - 1
        if (class(dots[[i]]) == "nri")
          goty <- i      
      }
    }
  }
  
  if (goty > 0)
  {
    y <- dots[[goty]]
    if (length(dots) > 1)
    {
      dots <- dots[-1*goty]
      eval_text <- ""
      for (i in 1:length(dots))
        eval_text <- paste(eval_text, ", ", names(dots)[i], " = dots$", names(dots)[i], sep = "")
      eval_text <- paste("t.test(",paste("x = x, y = y", eval_text, sep=""),")",sep="")
    } else {
      eval_text <- paste("t.test(",paste("x = x, y = y", sep=""),")",sep="")
    }
    mc <- paste("x =", deparse(substitute(x)), ", y =", deparse(substitute(list(...)[[goty]])))
  } else {
    if (length(dots) > 0)
    {
      eval_text <- ""
      for (i in 1:length(dots))
        eval_text <- paste(eval_text, ", ", names(dots)[i], " = dots$", names(dots)[i], sep = "")
      eval_text <- paste("t.test(",paste("x = x", eval_text, sep=""),")",sep="")
    } else {
      eval_text <- "t.test(x = x)"
    }
    mc <- paste("x =", deparse(substitute(x)))
  }
  if (class(y) == "Nri")
  {
    environment_apply <- new.env(parent = .GlobalEnv)
    assign("eval_text", eval_text, environment_apply)
    assign("x", x, environment_apply)
    assign("y", y, environment_apply)
    assign("dots", dots, environment_apply)
    
    dummy_matrix <- matrix(data = NaN, nrow = dim(x$nri)[1], ncol = dim(x$nri)[2])
    
    indices <- matrix(c(rep.int(1:dim(x$nri)[1],
                                dim(x$nri)[1])[row(dummy_matrix) > col(dummy_matrix)],
                        rep.int(1:(dim(x$nri)[1]-1),(dim(x$nri)[1]-1):1)),
                      ncol = 2
                    )
    ttest_data <- apply(indices, MARGIN = 1, FUN = ttest_apply)
    statistic <- distMat3D(ttest_data[1,], dim(x$nri)[1], 1)
    p.value <- distMat3D(ttest_data[2,], dim(x$nri)[1], 1)
    
    final <- list(statistic = statistic,
                  p.value = p.value)

    attr(final,"is.predictor.nri") <- FALSE
    attr(final,"function") <- "t.test"
    attr(final, "call") <- mc
    x@multivariate <- final
    return(x)
  } else {
    y <- try(as.factor(y), silent = TRUE)
    if (inherits(y, "try-error"))
      stop("Error while converting y to factor")
    if (nlevels(y)!=2)
      stop("Grouping factor must have exactly two levels")
  
    nri1 <- x
    nri2 <- x

    nri1@nri <- distMat3D(nri1@nri[,,y == levels(y)[1]])
    nri2@nri <- distMat3D(nri2@nri[,,y == levels(y)[2]])
    x <- nri1
    y <- nri2
    final <- eval(parse(text = eval_text))
    attr(final@multivariate, "call") <- mc
    return(final)
  }
}
)
