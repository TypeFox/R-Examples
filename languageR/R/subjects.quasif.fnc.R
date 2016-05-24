`subjects.quasif.fnc` <-
function(dat) {
  subjects = unique(dat[ , c("Subject", "SOA")])
  subjects = subjects[order(subjects$Subject, subjects$SOA), ]
  if ("RT" %in% colnames(dat)) {
   means = tapply(dat$RT, list(dat$Subject, dat$SOA), mean)
  } else {
   means = tapply(dat$RTsim, list(dat$Subject, dat$SOA), mean)
  }
  subjects$MeanRT = as.vector(t(means))
	model = lm(MeanRT ~ SOA * Subject, data = subjects)
  x = anova(model)
  p = 1 - pf(x["SOA", "Mean Sq"]/x["SOA:Subject", "Mean Sq"], 
	           x["SOA", "Df"], x["SOA:Subject", "Df"])
  return(list(p = p, data = subjects, model = x))
}

