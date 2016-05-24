#' @title Confidence Intervals (for one population mean)

#' @description An app to investigate how sample size and confidence level affect the width of a confidence interval.  
#' A sample is drawn from the input population and a confidence interval for
#' the population mean is calculated.  The kernel density plot for the population
#' and the histogram for each new sample are plotted, along with the confidence
#' interval.  Summary information is output to the console to tally the number of times
#' the computed confidence interval covers the true population mean and how many times
#' it misses.  There is an option to draw 100 or 1000 samples at a time.
#' 
#' @rdname CIMean
#' @usage CIMean(form,data)
#' @param form a formula of the form ~var.  
#' @param data A data frame from which var is drawn.
#' @return Graphical and numerical output
#' @export
#' @author Rebekah Robinson \email{rebekah_robinson@@georgetowncollege.edu}
#' @examples
#' \dontrun{
#' if (require(manipulate)) CIMean(~height,data=imagpop)
#' }
CIMean <- function (form, data) 
{
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  nprev = NULL
  lprev = NULL
  tally = c(0, 0)
  beginning = TRUE
  results = matrix(0, 2, 2, dimnames = list(c("Count", "Proportion"), 
                                            c("Covers", "Misses")))
  prsd = with(data, ParseFormula(form))
  varname = as.character(prsd$rhs)
  pop = data[, varname]
  N = length(pop)
  min.samp.size = 5
  max.samp.size = min(50, 0.01 * N)
  if (0.01 * N < 50) {
    print("Population size is too small to use with this app!")
  }
  pop.dens = density(pop)
  max.dens = pop.dens$y[which.max(pop.dens$y)]
  ymax = 1.4 * max.dens
  mu = mean(pop)
  stdev = sd(pop)
  manipulate(n = slider(min.samp.size, max.samp.size, step = 1, initial = 5, label = "Sample Size n"),
             conf.level = slider(60, 99, initial = 95, label = "Confidence Level"),
             sim.reps = picker("One at a time", "100 at a time", "1000 at a time", label = "Number of Repititions"), 
{
  current.sliders = list(n, conf.level)
  prev.sliders = list(nprev, lprev)
  if (!identical(current.sliders, prev.sliders)) {
    tally = c(0, 0)
  }
  if (sim.reps == "One at a time") {
    samp = sample(pop, n, replace = FALSE)
    color = ifelse(beginning, "transparent", "lightblue")
    border = ifelse(beginning, "transparent", "black")
    hist(samp, freq = FALSE, main = paste("Confidence Interval for Population Mean of", 
                                          varname), xlab = varname, col = color, border = border, 
         xlim = c(min(pop), max(pop)), 
         ylim = c(0, ymax))
    lines(density(pop), col = "red")
    abline(v = mu, col = "red", lwd = 2)
    if (!beginning) {
      conf = conf.level/100
      xbar = mean(samp)
      t.input = conf + ((1 - conf)/2)
      t = qt(t.input, df = n - 1)
      se = (sd(samp)/sqrt(n)) * (sqrt((N - n)/(N - 
                                                 1)))
      margin = t * se
      ci = c(xbar - margin, xbar + margin)
      segments(x0 = ci[1], y0 = 0, x1 = ci[2], y1 = 0, 
               col = "yellow", lwd = 2)
      points(xbar, 0, col = "blue", pch = 20)
      if (ci[1] < mu & ci[2] > mu) {
        tally[1] = tally[1] + 1
      }
      else {
        tally[2] = tally[2] + 1
      }
    }
    lprev <<- conf.level
    nprev <<- n
    if (!beginning) {
      if (sum(tally) == 0) {
        results[1, ] = round(tally, 0)
        results[2, ] = c(0, 0)
      }
      else {
        results[1, ] = round(tally, 0)
        results[2, ] = round(c((tally[1]/sum(tally)), 
                               (tally[2]/sum(tally))), 4)
      }
      print(results)
    }
    if (beginning) {
      cat("Density curve for population appears in red.\n")
      cat("Histogram of sample will appear in light blue.\n")
    }
    beginning = FALSE
  }
  if (sim.reps == "100 at a time") {
    conf = conf.level/100
    t.input = conf + ((1 - conf)/2)
    t = qt(t.input, df = n - 1)
    xbar = numeric(100)
    ci.left = numeric(100)
    ci.right = numeric(100)
    for (i in 1:100) {
      samp = sample(pop, n, replace = FALSE)
      
      se = (sd(samp)/sqrt(n)) * (sqrt((N - n)/(N - 
                                                 1)))
      margin = t * se
      xbar[i] = mean(samp)
      ci.left[i] = xbar[i] - margin
      ci.right[i] = xbar[i] + margin
      if (ci.left[i] < mu & ci.right[i] > mu) {
        tally[1] = tally[1] + 1
      }
      else {
        tally[2] = tally[2] + 1
      }
    }
    lprev <<- conf.level
    nprev <<- n
    if (sum(tally) == 0) {
      results[1, ] = round(tally, 0)
      results[2, ] = c(0, 0)
    }
    else {
      results[1, ] = round(tally, 0)
      results[2, ] = round(c((tally[1]/sum(tally)), 
                             (tally[2]/sum(tally))), 4)
    }
    print(results)
  }
  if (sim.reps == "1000 at a time") {
    conf = conf.level/100
    t.input = conf + ((1 - conf)/2)
    t = qt(t.input, df = n - 1)
    xbar = numeric(1000)
    ci.left = numeric(1000)
    ci.right = numeric(1000)
    for (i in 1:1000) {
      samp = sample(pop, n, replace = FALSE)
      
      se = (sd(samp)/sqrt(n)) * (sqrt((N - n)/(N - 
                                                 1)))
      margin = t * se
      xbar[i] = mean(samp)
      ci.left[i] = xbar[i] - margin
      ci.right[i] = xbar[i] + margin
      if (ci.left[i] < mu & ci.right[i] > mu) {
        tally[1] = tally[1] + 1
      }
      else {
        tally[2] = tally[2] + 1
      }
    }
    lprev <<- conf.level
    nprev <<- n
    if (sum(tally) == 0) {
      results[1, ] = round(tally, 0)
      results[2, ] = c(0, 0)
    }
    else {
      results[1, ] = round(tally, 0)
      results[2, ] = round(c((tally[1]/sum(tally)), 
                             (tally[2]/sum(tally))), 4)
    }
    print(results)  
  }
})
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("conf.level","sim.reps"))
