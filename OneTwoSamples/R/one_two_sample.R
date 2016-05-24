one_two_sample <- function(x, y = NULL, mu = c(Inf, Inf), sigma = c(-1, -1), var.equal = FALSE, ratio = 1, side = 0, alpha = 0.05){

## side and alternative: one-to-one correspondence
if (side<0) alternative = "less"
else if (side==0) alternative = "two.sided"
else alternative = "greater"

graphics.off()
if (is.null(y)){ ## One sample
	attr(x, "DNAME") = deparse(substitute(x))
	result = one_sample(x, mu = mu[1], sigma = sigma[1], side = side, alpha = alpha)
}
else{ ## Two samples
	cat("Interval estimation and test of hypothesis\n")
	## Test whether x and y are from the normal populations
	resx = shapiro.test(x); show(resx)
	resy = shapiro.test(y); show(resy)
	n1 = length(x); n2 = length(y)
	if (resx$p.value>alpha & resy$p.value>alpha){
		cat("\nx and y are both from the normal populations.\n")
		cat("\nx: descriptive statistics, plot, interval estimation and test of hypothesis\n")
		one_sample_x = one_sample(x, mu = mu[1], sigma = sigma[1], side = side, alpha = alpha)
		
		cat("\ny: descriptive statistics, plot, interval estimation and test of hypothesis\n")
		one_sample_y = one_sample(y, mu = mu[2], sigma = sigma[2], side = side, alpha = alpha)
		
		cat("\nInterval estimation and test of hypothesis of mu1-mu2\n")
		if (all(sigma>=0)){ ## sigma1, sigma2 are known
			cat("\nInterval estimation: interval_estimate5()\n")
			a = interval_estimate5(x, y, sigma = sigma, side = side, alpha = alpha); show(a)
			cat("\nTest of hypothesis: mean_test2()\n")
			b = mean_test2(x, y, sigma = sigma, side = side); show(b)
		}
		else{
			if (var.equal ==  TRUE){ ## sigma1=sigma2=sigma unknown
				cat("\nInterval estimation and test of hypothesis: t.test()\n")
				a = b = t.test(x, y, alternative = alternative, var.equal = var.equal, conf.level = 1-alpha); show(a)
			}
			else{ ## sigma1 != sigma2 unknown
				cat("\nInterval estimation and test of hypothesis: t.test()\n")
				a = b = t.test(x, y, alternative = alternative, var.equal = var.equal, conf.level = 1-alpha); show(a)
			}
		}

				
		cat("\nInterval estimation and test of hypothesis of sigma1^2/sigma2^2\n")
		if (all(mu < Inf)){ ## mu1, mu2 are known
			cat("Interval estimation: interval_var4()\n")
			c = interval_var4(x, y, mu = mu, side = side, alpha = alpha); show(c)
			cat("Test of hypothesis: var_test2()\n")
			d = var_test2(x, y, mu = mu, side = side); show(d)
		}
		else{ ## mu1 or mu2 is unknown
			cat("Interval estimation and test of hypothesis: var.test()\n")
			c = d = var.test(x, y, ratio = ratio, alternative = alternative, conf.level = 1-alpha); show(c)
		}
	}
	else {
		cat("\nx or y is not from the normal population.\n")
		one_sample_x = one_sample(x, mu = mu[1], sigma = sigma[1], side = side, alpha = alpha)
		one_sample_y = one_sample(y, mu = mu[2], sigma = sigma[2], side = side, alpha = alpha)
		a = b = c = d = NULL
	}
	
	if (n1 == n2){
		cat("n1 == n2\n")
		cat("\nTest whether x and y are from the same population\n") 
		cat("H0: x and y are from the same population (without significant difference)\n")
		cat("ks.test(x,y)\n")
		res.ks = ks.test(x,y); show(res.ks)
		cat("binom.test(sum(x<y), length(x))\n")
		res.binom = binom.test(sum(x<y), length(x)); show(res.binom)
		cat("wilcox.test(x, y, alternative = alternative, paired = TRUE)\n")
		res.wilcox = wilcox.test(x, y, alternative = alternative, paired = TRUE); show(res.wilcox)
		
		cat("\nFind the correlation coefficient of x and y\n") 
		cat("H0: rho = 0 (x, y uncorrelated)\n")
		cor.pearson = cor.test(x, y, alternative = alternative, method = "pearson", conf.level = 1-alpha); show(cor.pearson)
		cor.kendall = cor.test(x, y, alternative = alternative, method = "kendall", conf.level = 1-alpha); show(cor.kendall)
		cor.spearman = cor.test(x, y, alternative = alternative, method = "spearman", conf.level = 1-alpha); show(cor.spearman)
	}
	else{
		cat("n1 != n2\n")
		cat("\nTest whether x and y are from the same population\n")
		cat("H0: x and y are from the same population (without significant difference)\n")
		cat("ks.test(x,y)\n")
		res.ks = ks.test(x,y); show(res.ks)
		res.binom = NULL
		cat("wilcox.test(x, y, alternative = alternative)\n")
		res.wilcox = wilcox.test(x, y, alternative = alternative); show(res.wilcox)
		cor.pearson = cor.kendall = cor.spearman = NULL
	}
	
	result = list(one_sample_x = one_sample_x, one_sample_y = one_sample_y, 
			  mu1_mu2_interval = a, mu1_mu2_hypothesis = b, 
			  sigma_ratio_interval = c, sigma_ratio_hypothesis = d,
			  res.ks = res.ks, res.binom = res.binom, res.wilcox = res.wilcox,
			  cor.pearson = cor.pearson, cor.kendall = cor.kendall, cor.spearman = cor.spearman)
}
}