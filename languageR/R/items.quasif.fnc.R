`items.quasif.fnc` <-
function(dat) {
  items = unique(dat[ , c("Item", "SOA")])
  items = items[order(items$Item, items$SOA), ]
  if ("RT" %in% colnames(dat)) {
     means = tapply(dat$RT, dat$Item, mean)
  } else {
     means = tapply(dat$RTsim, dat$Item, mean)
  }
  items$MeanRT = as.vector(means)
	model = lm(MeanRT ~ SOA, items)
  p = summary(model)$coef["SOAshort", "Pr(>|t|)"]
  return(list(p = p, data = items, model = model))
}

