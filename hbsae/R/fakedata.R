#' Generate artificial dataset for demonstration and testing purposes.
#'
#' @export
#' @param M number of areas.
#' @param meanNarea mean number of population units per area.
#' @param sW within area standard deviation.
#' @param sB between area standard deviation.
#' @param sBx random slope standard deviation.
#' @param samplingFraction sampling fraction used to draw a random sample from the population units.
#' @return List containing sample (sam), population totals (Xpop), and true population means for four
#'         target variables (mY0, mY1, mY2, mY3).
generateFakeData <- function(M=50, meanNarea=1000, sW=100, sB=50, sBx=0.5, samplingFraction=0.1) {
  N.area <- pmax(1, trunc(meanNarea * rgamma(M, 1)))
  pop <- data.frame(area = as.factor(rep.int(1:M, N.area)))  # variable defining small areas
  N <- nrow(pop)  # population size
  M2 <- 5  # number of broad areas
  pop$area2 <- as.factor(rep.int(sample(1:M2, M, replace=TRUE), N.area))  # generate broad area variable
  pop$x <- abs(rnorm(N, sqrt(N:1), sqrt(N:1)/4))  # generate quantitive covariate
  pop$x.area <- (rowsum(pop$x, pop$area)/N.area)[pop$area, ]  # area-level covariate
  # generate target variables
  beta.x <- 1  # fixed effect for x
  beta.xarea <- 0.2  # fixed effect for x.area
  beta.area2 <- c(-50, 0, 50, 100, 150)  # fixed effect for area2
  re <- rnorm(M, sd=sB)  # generate random area effects
  re.x <- rnorm(M, sd=sBx)  # generate random slopes
  e <- rnorm(N, sd=sW)  # residual errors
  pop$y0 <- beta.x * pop$x + beta.xarea * pop$x.area + beta.area2[pop$area2] + re[pop$area] + e  # y0 generated according to basic unit-level model
  pop$y1 <- beta.x * pop$x + beta.xarea * pop$x.area + beta.area2[pop$area2] + re[pop$area] + e*(sqrt(pop$x) / mean(sqrt(pop$x)))  # include unit-level heteroscedasticity
  pop$y2 <- beta.x * pop$x + beta.xarea * pop$x.area + beta.area2[pop$area2] + re[pop$area]*(pop$x.area / mean(pop$x.area)) + e  # include area-level heteroscedasticity
  pop$y3 <- beta.x * pop$x + beta.xarea * pop$x.area + beta.area2[pop$area2] + re[pop$area] + re.x[pop$area] * pop$x + e  # include random slopes

  # true area means
  mY0 <- rowsum(pop$y0, pop$area)/N.area
  mY1 <- rowsum(pop$y1, pop$area)/N.area
  mY2 <- rowsum(pop$y2, pop$area)/N.area
  mY3 <- rowsum(pop$y3, pop$area)/N.area

  # build population totals
  Xpop <- rowsum(model.matrix(~ x + x.area + area2 + area2*x, pop), pop$area)

  # draw a random sample, on which population inferences will be based
  sam <- pop[sample(1:N, trunc(samplingFraction * N)), ]
  
  list(sam=sam, Xpop=Xpop, mY0=mY0, mY1=mY1, mY2=mY2, mY3=mY3)
}
