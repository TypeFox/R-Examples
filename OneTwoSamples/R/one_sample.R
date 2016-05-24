one_sample <- function(x, mu = Inf, sigma = -1, side = 0, alpha = 0.05){

## 
## Descriptive statistics, plot
## 
DNAME = if (!is.null(attr(x, "DNAME"))) attr(x, "DNAME") else deparse(substitute(x))
cat(paste("quantile of ", DNAME, sep = "")); cat("\n"); print(quantile(x))
cat(paste("data_outline of ", DNAME, sep = "")); cat("\n"); print(data_outline(x))

## Test whether x is from the normal population
res = shapiro.test(x)
res$data.name = DNAME
show(res)
if (res$p.value>alpha){
	cat("\nThe data is from the normal population.\n")
}
else {
	cat("\nThe data is not from the normal population.\n")
}

## Histograms with density estimation curve and normal density curve
w<-seq(min(x),max(x),length.out = 51)
Vector = c(density(x)$y, dnorm(w, mean(x), sd(x)))
ylim = c(min(Vector), max(Vector))

dev.new(); hist(x, freq = FALSE, ylim = ylim, main = paste("Histogram of" , DNAME), xlab = DNAME)
lines(density(x),col="blue",lty = 1)
lines(w, dnorm(w, mean(x), sd(x)), col="red",lty = 2)
leg.txt = c("Density estimation curve","Normal density curve")
legend("topleft",legend = leg.txt,lty = 1:2,col = c('blue','red'))

## Empirical cumulative distribution function (ECDF) vs normal cdf
dev.new(); plot(ecdf(x),verticals = TRUE, do.p = FALSE, main = paste("ecdf(" , DNAME, ")", sep = ""), xlab = DNAME, ylab = paste("Fn(" , DNAME, ")", sep = ""))
w<-seq(min(x),max(x),length.out = 51)
lines(w, pnorm(w, mean(x), sd(x)), col="red")

## QQ plot
dev.new(); qqnorm(x); qqline(x)

## 
## Interval estimation and test of hypothesis
## 
## side, alternative, H0sign, H1sign: one-to-one correspondence
if (side<0) {
	alternative = "less"
	H0sign = ">="
	H1sign = "<"
}
else if (side==0) {
	alternative = "two.sided"
	H0sign = "="
	H1sign = "!="
}
else {
	alternative = "greater"
	H0sign = "<="
	H1sign = ">"
}

if (res$p.value>alpha){
	cat("\nThe data is from the normal population.\n")
	
	cat("\nInterval estimation and test of hypothesis of mu\n")
	mu0 = if (mu < Inf) mu else 0
	if (sigma >= 0){ ## sigma known
		cat("Interval estimation: interval_estimate4()\n")
		a = interval_estimate4(x, sigma = sigma, side = side, alpha = alpha); show(a)
		cat("Test of hypothesis: mean_test1()\n")
		cat("H0: mu", H0sign, mu0, "    "); cat("H1: mu", H1sign, mu0, "\n")
		b = mean_test1(x, mu = mu0, sigma = sigma, side = side); show(b)
	}
	else{ ## sigma unkown
		cat("Interval estimation and test of hypothesis: t.test()\n")
		cat("H0: mu", H0sign, mu0, "    "); cat("H1: mu", H1sign, mu0, "\n")
		a = b = t.test(x, mu = mu0, alternative = alternative, conf.level = 1-alpha); a$data.name = DNAME; show(a)
	}
	
	cat("\nInterval estimation and test of hypothesis of sigma\n")
	cat("Interval estimation: interval_var3()\n")
	c = interval_var3(x, mu = mu, side = side, alpha = alpha); show(c)
	cat("Test of hypothesis: var_test1()\n")
	cat("H0: sigma2", H0sign, sigma^2, "    "); cat("H1: sigma2", H1sign, sigma^2, "\n")
	d = var_test1(x, sigma2 = sigma^2, mu = mu, side = side); show(d)
	
	result = list(mu_interval = a, mu_hypothesis = b, sigma_interval = c, sigma_hypothesis = d)
}
else {
	cat("\nThe data is not from the normal population.\n")
	## sigma known and unknown
	cat("\nInterval estimation of mu: interval_estimate3()\n")
	a = interval_estimate3(x, sigma = sigma, alpha = alpha); show(a)
	cat("\nTest of hypothesis of mu: NA\n")
	cat("\nInterval estimation and test of hypothesis of sigma: NA\n")
	result = list(mu_interval = a, mu_hypothesis = NA, sigma_interval = NA, sigma_hypothesis = NA)
}
}

