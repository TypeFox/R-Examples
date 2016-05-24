pvalFmla <-
function(data, i, covariable=NULL, family=binomial) {
		if(length(unique(data[,i])) == 1) { 
			pval <- 1
      	} else {
			form=CreateFormula(data,i,covariable)
			if(length(covariable)==0) { 
				fit <- glm(formula=form, data=data, family=family)
				pval <- anova(fit, test="LRT")[2,5]
				return(pval)
			}
			fit <- glm(formula=form, data=data.frame(data, covariable), family=family)
		  	pval <- anova(fit, test="LRT")[2,5]
		}
		return(pval)
}
