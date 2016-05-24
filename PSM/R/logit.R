`logit` <-
function(x,xmin=0,xmax=1)
  log((x-xmin)/(xmax-x))

