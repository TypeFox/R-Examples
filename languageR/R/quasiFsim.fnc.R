`quasiFsim.fnc` <-
function(dat) {
  if ("RT" %in% colnames(dat)) {
     dat.lm = lm(RT ~ SOA + Item + Subject + SOA:Subject + Item:Subject, 
	            data = dat)
  } else {
     dat.lm = lm(RTsim ~ SOA + Item + Subject + SOA:Subject + Item:Subject, 
	            data = dat)
  }

  x = anova(dat.lm)
  qF = quasiF.fnc(x["SOA","Mean Sq"], x["Item:Subject", "Mean Sq"],
         x["SOA:Subject", "Mean Sq"], x["Item", "Mean Sq"],
         x["SOA","Df"], x["Item:Subject", "Df"],
         x["SOA:Subject", "Df"], x["Item", "Df"])
  return(list(p = qF$p,  data = dat, model = dat.lm, quasiF = qF))
}

