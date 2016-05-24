#' @import methods stats

harmonic.mean <- function(x) 1 / mean(1 / x)

fround <- function(x, digits=4) sprintf(paste("%.",digits,"f",sep=""), x)

get.p.value <- function(statistic, alternative="two.sided", distribution="norm", df=NULL) {
  if(alternative == "two.sided") statistic <- abs(statistic)

  dist.fun <- paste("p", tolower(distribution), sep="")

  params <- list(statistic)
  if(!is.null(df)) params <- c(params, df)
  x <- do.call(dist.fun, params)

  switch(alternative,
    two.sided=2 * (1 - x),
    greater=1 - x,
    less=x
  )
}

get.conf.int <- function(alpha, n, items, conf.level=.95) {
  conf.int <- matrix(nrow=length(alpha), ncol=2, dimnames=list(paste("CI", 1:length(alpha), sep=""),c("lower.bound","upper.bound")))
  for(i in 1:length(alpha)) {
    conf.int[i,] <- cronbach.alpha.CI(alpha[i], n[i], items[i], conf.level=conf.level)
  }
  conf.int
}

check.alternative <- function(alternative) {
  switch(tolower(substr(alternative,1,1)),
    t="two.sided",
    g="greater",
    l="less",
    stop("The parameter 'alternative' must be either 'two.sided', 'greater', or 'less' (or just the initial letter).")
  )
}

print.alternative <- function(object) {
  alphas <- paste(paste(paste("a", 1:(length(object@alpha) - 1), sep=""), collapse=", "), if(length(object@alpha) > 2) ",", " and a", length(object@alpha), sep="")

  switch(class(object),
    cocron.n.coefficients={
      paste("Null hypothesis: ", alphas, " are equal\n",
      "Alternative hypothesis: ", alphas, " are not equal", sep="")
    },
    cocron.two.coefficients={
      paste("Null hypothesis: a1 is equal to a2\n",
      "Alternative hypothesis: a1", switch(object@alternative, two.sided=" is not equal to ", greater=" is greater than ", less=" is less than "),
      "a2 (", switch(object@alternative, two.sided="two", "one"), "-sided)", sep="")
    }
  )
}

print.test.statistic <- function(object) {
  s <- {
    if(!is.null(object@statistic)) {
      evaluation <- {
        if(object@p.value <= object@los) "Null hypothesis rejected"
        else "Null hypothesis retained"
      }
      df <- {
        if(is.null(object@df)) ""
        else if(length(object@df) == 1) paste(", df = ", object@df, sep="")
        else if(length(object@df) == 2) paste(", df1 = ", object@df[1], ", df2 = ", object@df[2], sep="")
      }
      paste(object@distribution, " = ", fround(object@statistic), df, ", p-value = ", fround(object@p.value), "\n", evaluation, sep="")
    } else NULL
  }

  s
}

setClass("cocron.n.coefficients",
  representation(
    alpha="numeric",
    n="numeric",
    items="numeric",
    dep="logical",
    r="matrix",
    conf.int="matrix",
    los="numeric",
    alternative="character",
    conf.level="numeric",
    df="numeric",
    statistic="numeric",
    distribution="character",
    p.value="numeric"
  )
)

setClass("cocron.two.coefficients",
  representation(
    alpha="numeric",
    n="numeric",
    dep="logical",
    r="numeric",
    los="numeric",
    alternative="character",
    df="numeric",
    statistic="numeric",
    distribution="character",
    p.value="numeric"
  )
)

print.dependence <- function(object) {
  if(object@dep) "The coefficients are based on dependent groups"
  else "The coefficients are based on independent groups"
}

print.cocron <- function(object) {
  two.or.n <- switch(class(object),
    cocron.two.coefficients="two",
    cocron.n.coefficients="n"
  )

  cat("\n  Compare ", two.or.n, " alpha coefficients\n\n",
  "Comparison between: ", paste(paste("a", 1:length(object@alpha), " = ", round(object@alpha,4), sep=""), collapse=", "), "\n",
  print.dependence(object), "\n", sep="")
  if(class(object) == "cocron.n.coefficients") cat(object@conf.level * 100, "% confidence intervals: ", paste(paste("CI", 1:length(object@alpha), " = ", apply(object@conf.int, 1, function(x) paste(fround(x,4), collapse=" ")), sep=""), collapse=", "), "\n", sep="")
  cat("Group sizes: ", paste(paste("n", 1:length(object@n), " = ", object@n, sep=""), collapse=", "), "\n", sep="")
  if(class(object) == "cocron.n.coefficients") cat("Item count: ", paste(paste("i", 1:length(object@items), " = ", object@items, sep=""), collapse=", "), "\n", sep="")
  cat(print.alternative(object), "\n",
  "Level of significance: ", object@los, "\n\n", sep="")
  cat(print.test.statistic(object), "\n", sep="")
}

setMethod("show", "cocron.n.coefficients",
  function(object) {
    print.cocron(object)
  }
)

setMethod("show", "cocron.two.coefficients",
  function(object) {
    print.cocron(object)
  }
)
