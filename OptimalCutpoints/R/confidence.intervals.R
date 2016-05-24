confidence.intervals <-
function(Se, Sp, PPV, NPV, DLR.Positive, DLR.Negative, pop.prev, n, control, conf.level) { 
	TP <- Se*n$d
	FP <-(1-Sp)*n$h
	TN <- Sp*n$h
	FN <-(1-Se)*n$d
	
	# Sensitivity and Specificity	 
 	if (control$ci.SeSp == "Exact") {
 		ci.Se <- ci.exact(x = FN, y = TP, accuracy.measure = "Sensitivity", conf.level = conf.level)
 		ci.Sp <- ci.exact(x = FP, y = TN, accuracy.measure = "Specificity", conf.level = conf.level)
    } else if (control$ci.SeSp == "Quadratic") {
 		ci.Se <- ci.quadratic(TP, FN, accuracy.measure = "Sensitivity", conf.level = conf.level)
 		ci.Sp <- ci.quadratic(TN, FP, accuracy.measure = "Specificity", conf.level = conf.level)
	} else if (control$ci.SeSp == "Wald") {
 		ci.Se <- ci.wald(FN, TP, accuracy.measure = "Sensitivity", measure = Se, n$d, conf.level = conf.level)      
 		ci.Sp <- ci.wald(FP, TN, accuracy.measure = "Specificity", measure = Sp, n$h, conf.level = conf.level)
	} else if (control$ci.SeSp == "AgrestiCoull") {
 		ci.Se <- ci.AgrestiCoull(measure = Se, n$d, conf.level = conf.level)
 		ci.Sp <- ci.AgrestiCoull(measure = Sp, n$h, conf.level = conf.level)
	} else if (control$ci.SeSp == "RubinSchenker") {
 		ci.Se <- ci.RubinSchenker(TP, n$d, conf.level = conf.level) 
 		ci.Sp <- ci.RubinSchenker(TN, n$h, conf.level = conf.level)
 	}

 	# PPV and NPV
 	if (control$ci.PV == "Exact") {
 		ci.PPV <- ci.exact(x = FN, y = TP, accuracy.measure = "Positive Predictive Value", z = FP, t = TN, conf.level = conf.level)
 		ci.NPV <- ci.exact(x = FN, y = TP, accuracy.measure = "Negative Predictive Value", z = FP, t = TN, conf.level = conf.level)
 	} else if (control$ci.PV == "Quadratic") {
 		ci.PPV <- ci.quadratic(TP, FP, accuracy.measure = "Positive Predictive Value", conf.level = conf.level)
 		ci.NPV <- ci.quadratic(TN, FN, accuracy.measure = "Negative Predictive Value", conf.level = conf.level)
 	} else if (control$ci.PV == "Wald") {
 		ci.PPV <- ci.wald(TP, FP, accuracy.measure = "Positive Predictive Value", measure = PPV, TP+FP, conf.level = conf.level)
 		ci.NPV <- ci.wald(TN, FN, accuracy.measure = "Negative Predictive Value", measure = NPV, TN+FN, conf.level = conf.level)
 	} else if (control$ci.PV == "AgrestiCoull") {
 		ci.PPV <- ci.AgrestiCoull(measure = PPV, TP+FP, conf.level = conf.level)
 		ci.NPV <- ci.AgrestiCoull(measure = NPV, TN+FN, conf.level = conf.level)
 	} else if (control$ci.PV == "RubinSchenker") {
 		ci.PPV <- ci.RubinSchenker(TP, TP+FP, conf.level = conf.level)
 		ci.NPV <- ci.RubinSchenker(TN, TN+FN, conf.level = conf.level)
 	} else if (control$ci.PV == "Transformed") {
 		ci.PPV <- list(ci = 1/(1+((1-pop.prev)/(pop.prev*ci.transformed(Se, 1-Sp, n, conf.level = conf.level)$ci))))
 		ci.NPV <- list(ci = 1/(1 + (pop.prev/(1-pop.prev))*ci.transformed(1-Se, Sp, n, conf.level = conf.level)$ci))	
 	} else if (control$ci.PV == "NotTransformed") {  	
 	   	ci.PPV <- list(ci = 1/(1+((1-pop.prev)/(pop.prev*ci.NotTransformed(Se, 1-Sp, DLR.Positive, n, conf.level)$ci))))
 		ci.NPV <- list(ci = 1/(1+(pop.prev/(1-pop.prev))*ci.NotTransformed(1-Se, Sp, DLR.Negative, n, conf.level)$ci))
	} else if (control$ci.PV == "GartNam") {  	
 		ci.PPV <- list(ci = 1/(1+((1-pop.prev)/(pop.prev*ci.GartNam(Se, 1-Sp, n, conf.level)$ci))))
 		ci.NPV <- list(ci = 1/(1+((pop.prev*ci.GartNam(1-Se, Sp, n, conf.level)$ci)/(1-pop.prev))))
 	}
 	
 	# DLRs
 	if (control$ci.DLR == "Transformed") {
 		ci.DLR.positive <- ci.transformed(Se, 1-Sp, n, conf.level = conf.level)
 		ci.DLR.negative <- ci.transformed(1-Se, Sp, n, conf.level = conf.level) 
 	}     
    	else if (control$ci.DLR == "NotTransformed") {
 		ci.DLR.positive <- ci.NotTransformed(Se, 1-Sp, DLR.Positive, n, conf.level = conf.level) 
 		ci.DLR.negative <- ci.NotTransformed (1-Se, Sp, DLR.Negative, n, conf.level = conf.level) 
 	} 
    	else if (control$ci.DLR == "GartNam") {
 		ci.DLR.positive <- ci.GartNam(Se, 1-Sp, n, conf.level = conf.level)
 		ci.DLR.negative <- ci.GartNam(1-Se, Sp, n, conf.level = conf.level)
 	}
 	  
 	res <- list(ci.Se = ci.Se$ci, ci.Sp = ci.Sp$ci, ci.PPV = ci.PPV$ci, ci.NPV = ci.NPV$ci, ci.DLR.positive = ci.DLR.positive$ci, ci.DLR.negative = ci.DLR.negative$ci)
 	res
}
