# From "Belenky, Greg" <belenky@wsu.edu>
rtost <- function (x, y = NULL, alpha = 0.05, epsilon = 0.31, tr = 0.2, ...) {
  tt <- yuen.t.test(x, y, conf.level = (1 - 2 * alpha), tr = tr, ...)
  if (tt$conf.int[2] < epsilon & tt$conf.int[1] > -epsilon)
    result <- "rejected"
  else result <- "not rejected"
  mean.diff <- ifelse(length(tt$estimate) == 1, as.numeric(tt$estimate),
                      as.numeric(tt$estimate[1] - tt$estimate[2]))
  se.diff <- as.numeric((tt$conf.int[2] - tt$conf.int[1])/
                        qt(1 - alpha, df = tt$parameter))/2
  pv <- as.numeric(pt((epsilon - abs(mean.diff))/se.diff, tt$parameter,
                      lower.tail = FALSE))
  if (pv < 0.5)
    test.ci <- yuen.t.test(x, y, conf.level = (1 - 2 * pv), tr = tr, ...)$conf.int
  else test.ci <- NULL
  return(list(mean.diff = mean.diff, se.diff = se.diff, alpha = alpha,
              ci.diff = tt$conf.int, df = tt$parameter, epsilon = epsilon,
              result = result, p.value = pv, check.me =test.ci))
}


