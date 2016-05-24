	ci.exact <-
	function(x, y, accuracy.measure, z, t, conf.level) {
		F1 <- suppressWarnings(qf(1-((1-conf.level)/2), 2*(x+1), 2*y))
		F2 <- suppressWarnings(qf(1-((1-conf.level)/2), 2*(y+1), 2*x))
		      
		if (accuracy.measure == "Negative Predictive Value" | accuracy.measure == "Positive Predictive Value") {
	      	F3 <- suppressWarnings(qf(1-((1-conf.level)/2), 2*(z+1), 2*t))
	      	F4 <- suppressWarnings(qf(1-((1-conf.level)/2), 2*(t+1), 2*z))
	    }	    
	    if (accuracy.measure == "Sensitivity" | accuracy.measure == "Specificity") {
	      	ll <- y/(((x+1)*F1)+y)
	      	ul <- ((y+1)*F2)/(x+(y+1)*F2)
	    }
	    if (accuracy.measure == "Positive Predictive Value") {   
	      	ll <- y/(y+(z+1)*F3)
	      	ul <- ((y+1)*F2)/(z+(y+1)*F2)
	    }           
	    if(accuracy.measure == "Negative Predictive Value") {
	      	ll <- t/(t+(x+1)*F1)
	      	ul <-((t+1)*F4)/(x+(t+1)*F4)
		}
		if(((accuracy.measure == "Sensitivity" | accuracy.measure == "Specificity") & any(is.nan(F1), is.nan(F2)))
		   || (accuracy.measure == "Positive Predictive Value" & any(is.nan(F3), is.nan(F2)))
		   || (accuracy.measure == "Negative Predictive Value" & any(is.nan(F1), is.nan(F4)))) {
			warning(paste(accuracy.measure, " CI: \"Exact\" method may not be valid for some values (see Help Manual).\n", sep = ""), call. = FALSE, immediate. = TRUE)
		}
		res <- list (ci = matrix(c(ll,ul), ncol = 2))           
	}
