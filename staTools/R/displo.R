#' Discrete Powerlaw Object
#'
#' This function allows to create a discrete powerlaw object to analyze.
#' @param x A vector containing the observations.
#' @param summary Logical, whether print a summary with some information concerning data. By default is set to TRUE.
#' @keywords discrete powerlaw
#' @export displo
#' @examples
#' data(moby)
#' x = moby
#' o = displo(x)

displo = function(x, summary = TRUE)
{
  o = new.env()

  # check data
  x = floor(x)
  x = x[x>0]

  o$x = sort(x)
  o$nx = length(x)
  o$ux = sort(unique(x))
  o$nux = length(o$ux)
  o$p = 1 - cdf(o$x)$y
  o$xmin = min(x)
  o$xmax = max(x)
  o$alpha = numeric()
  o$sigma = numeric()
  o$xmins = NULL
  o$alphas = NULL
  t = as.data.frame(table(x))
  o$freq = t$Freq
  o$cumfreq = cumsum(t$Freq)

  # fit exponential
  o$fitexp = list()
  o$fitexp[["rate"]] = 1/mean(o$x)

  # fit poisson
  o$fitpois = list()
  o$fitpois[["lambda"]] = mean(o$x)

  # fit lognormal
  o$fitlnorm = list()
  o$fitlnorm[["meanlog"]] = sum(log(o$x))/o$nx
  o$fitlnorm[["sdlog"]] = sqrt(sum((log(o$x) - o$fitlnorm[["meanlog"]])^2)/o$nx)

  # print summary
  if (summary)
  {
    cat("\nDiscrete Powerlaw Object",
        "\n************************",
         "\nn =", o$nx,
         "\nn (unique) =", o$nux,
         "\nxmin =", o$xmin,
         "\nxmax =", o$xmax,
         "\nalpha =", o$alpha,
         "\n************************"
        )
  }
  return(o)
}
