jury <- c(780, 117, 114, 384, 58)
chisq.test(jury, p = c(.54, .18, .12, .15, .01))
xchisq.test(jury, p = c(.54, .18, .12, .15, .01)) # to list expected counts

