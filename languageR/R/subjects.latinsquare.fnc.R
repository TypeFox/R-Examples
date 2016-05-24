`subjects.latinsquare.fnc` <-
function(dat) {
  subjects = unique(dat[ ,c("Group", "Subject", "SOA", "List")])
  subjects = subjects[order(subjects$Subject, subjects$List), ]
  if ("RT" %in% colnames(dat)) {
    subjects$MeanRT = as.vector(t(tapply(dat$RT, 
      list(dat$Subject, dat$List), mean)))
  } else {
    subjects$MeanRT = as.vector(t(tapply(dat$RTsim, 
      list(dat$Subject, dat$List), mean)))
  }
  subjects.lm = lm(MeanRT ~ Group/Subject + SOA*List, 
    data = subjects)
  x = anova(subjects.lm)
  p = 1 - pf(x["SOA", "Mean Sq"]/x["SOA:List", "Mean Sq"], 
    x["SOA", "Df"], x["SOA:List", "Df"])
  return(list(p = p, data = dat, model = subjects.lm))
}

