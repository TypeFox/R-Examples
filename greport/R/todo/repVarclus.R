#' Variable Clustering Diagrams
#'
#' Generate variable clustering diagrams for each time period.
#'
#' Variables are grouped according to how they are
#' correlated with one another, as measured by the square of the
#' Spearman \eqn{\rho} rank correlation coefficient computed on all
#' pairs of variables.  Variables connected on
#' lower branches are more highly correlated with one another.
#' Variables missing in more than 75% of the observations or
#' categorical variables having more than 20 levels are ignored.
#' Categories less than 0.1 prevalent are pooled with other rare categories.
#'
#' @param data data.frame. Data used for report.
#' @param time numeric vector. Time for each record.
#' @param times numeric vector. Subset of times to use.
#' @param nmin numeric. Variables must have four responses or at least
# 'two responses with more than \code{nmin} counts to be included.
#' @export
#' @examples
#' \dontrun{
#'   dat <- data.frame(age=sample(41:80, replace=TRUE, 1000), height=round(rnorm(1000, 68, 4)))
#'   dat$weight <- sample(140:200, replace=TRUE, 1000)+(dat$height-60)*3
#'   repVarclus(dat, numeric(1000))
#' }

repVarclus <- function(data, time, times=sort(unique(time)), nmin=10) {
  g <- function(y, nmin) {
    y <- table(oldUnclass(y))
    length(y) > 3 || (length(y) > 1 && sort(y)[length(y)-1] > nmin)
  }

  if(any(times==0)) {
      cat('A variable clustering diagram is shown in the figure',
          'which follows.  Variables are grouped according to how they',
          'are correlated with one another, as measured by the square',
          'of the Spearman $\\rho$ rank correlation coefficient',
          'computed on all pairs of variables.  Variables connected on',
          'lower branches are more highly correlated with one another.',
          'Variables missing in more than 0.75 of the observations or',
          'categorical variables having more than 20 levels are ignored.',
          'Categories less than 0.1 prevalent are pooled with other rare categories.\n',
          file=file.path(TexDirName(), 'Ovarclus.tex'))
      startPlot('Ovarclus', h=4)
      d <- data[time==0,]
      r <- sapply(d, g, nmin=nmin)
      d <- d[r]
      v <- varclus(~., data=d, fracmiss=.75, maxlevels=20, minprev=.1)
      plot(v)
      putFig('Ovarclus', 'Ovarclus', 'Clustering of variables at baseline', append=TRUE)
      endPlot()
  }
  cat('Variable clustering diagrams are shown in the',
      if(length(times) > 1)'figures that follow.' else 'figure that follows.', 
      'Variables are grouped according to how they are',
      'correlated with one another, as measured by the square of the',
      'Spearman $\\rho$ rank correlation coefficient computed on all',
      'pairs of variables.  Variables connected on',
      'lower branches are more highly correlated with one another.',
      'Variables missing in more than 0.75 of the observations or',
      'categorical variables having more than 20 levels are ignored.',
      'Categories less than 0.1 prevalent are pooled with other rare categories.\n',
      file=file.path(TexDirName(), 'varclus.tex'))
  startPlot('varclus%d', h=4)
  i <- 0
  for(x in times) {
      i <- i + 1
      d <- data[time==x,]
      r <- sapply(d, g, nmin=nmin)
      d <- d[r]
      v <- varclus(~., data=d, fracmiss=.75, maxlevels=20, minprev=.1)
      plot(v)
      title(paste(label(time), x))
      putFig('varclus', paste('varclus',i,sep=''),
             paste('Clustering of variables at',label(time),x),
             append=TRUE)
  }
  endPlot()
  invisible()
}
