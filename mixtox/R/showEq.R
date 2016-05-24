showEq <- function(eq){
	# default Non-linear least square fitting algorithm is the Gauss-Newton algorithm
	print('######################## function selection ###########################')
	#print(c("Hill", "Hill_two", "Hill_three", "Hill_four", "Hill_six", "Weibull", "Weibull_three", "Weibull_four", "Logit", "Logit_three",
	#		"Logit_four", "BCW(Box-Cox-Weibull)", "BCL(Box-Cox-Logit)", "GL(Generalized Logit)", "Brain_Consens", "BCV",
	#		"Biphasic"))
	#eq <- readline('input equation name: ')
	if(missing(eq)) eq <- 'sigmoid'
	fun <- switch(eq,
		# Howard GJ, Webster TF. 2009. Generalized concentration addition: A method for examining mixtures containing partial agonists. J. Theor. Biol. 259:469-477
		Hill = print("Hill: y ~ 1 / (1 + (Alpha / x)^Beta)"),
		# For Hill equation: Alpha = EC50; Beta = m(Hill coefficient); Gamma = Top; Delta = Bottom
		# Hill function with slope parameter 1. Alpha is EC50 here.
		Hill_two = print("Hill_two: y ~ Beta * x / (Alpha + x)"),
		Hill_three = print("Hill_three: y ~ Gamma /(1 + (Alpha / x)^Beta)"),
		Hill_four = print("Hill_four: y ~ Delta + (Gamma - Delta) / (1 + (Alpha / x)^Beta)"),
		Hill_six = print("Hill_six: y ~ (Gamma / (1 + (Alpha / x)^Beta)) * (Gamma_one / (1 + (Alpha_one / x)^Beta_one))"),
		Weibull = print("Weibull: y ~ 1 - exp(-exp(Alpha + Beta * log10(x)))"),
		Weibull_three = print("Weibull_three: y ~ Gamma * (1 - exp(-exp(Alpha + Beta * log10(x))))"),
		Weibull_four = print("Weibull_four: y ~ Gamma + (Delta - Gamma) * exp(-exp(Alpha + Beta * log10(x)))"),
		Logit = print("Logit: y ~ 1/(1 + exp((-Alpha)- Beta * log10(x)))"),
		Logit_three = print("Logit_three: y ~ Gamma / (1 + exp((-Alpha) - Beta * log10(x)))"),
		Logit_four = print("Logit_four: y ~ Delta + (Gamma - Delta) / (1 + exp((-Alpha) - Beta * log10(x)))"),
		BCW = print("BCW(Box-Cox-Weibull): y ~ 1 - exp(-exp(Alpha + Beta * ((x^Gamma - 1) / Gamma)))"),
		BCL = print("BCL(Box-Cox-Logit): y ~ (1 + exp(-Alpha - Beta *((x^Gamma - 1) / Gamma)))^(-1)"),
		GL = print("GL(Generalized Logit): y ~ 1 / (1 + exp(-Alpha - Beta * log10(x)))^Gamma"),
		# An equation to describe dose responses where there isstimulatin of growth at low doses. 1989. Weed Research.
		Brain_Consens = print("Brain_Consens: y ~ 1 - (1 + Alpha * x) / (1 + exp(Beta * Gamma) * x^Beta)"),
		# Vanewijk, P. H. and Hoekstra, J.A. Calculation of the EC50 and its confidence interval when subtoxic stimulus is present. 1993, Ecotoxicol. Environ. Saf.
		BCV = print("BCV: y ~ 1 - Alpha * (1 + Beta * x) / (1 + (1 + 2 * Beta * Gamma) * (x / Gamma)^Delta)"),
		# Cedergreen, N., Ritz, C., Streibig, J.C., 2005. Improved empirical models describing hormesis. Environ. Toxicol. Chem. 24, 3166-3172
		Cegergreen = print("Cedergreen: y ~ 1 - (1 + Alpha * exp(-1 / (x^Beta))) / (1 + exp(Gamma * (log(x) - log(Delta))))"),
		# Beckon, W. et.al. 2008. A general approach to modeling biphasic relationships. Environ. Sci. Technol. 42, 1308~1314.
		#Beckon = print("Beckon: y ~ (Alpha + (1 - Alpha / (1 + (Beta / x)^Gamma))) / (1 + (x / Delta)^Epsilon)"),
		# Zhu X-W, et.al . 2013. Modeling non-monotonic dose-response relationships: Model evaluation and hormetic quantities exploration. Ecotoxicol. Environ. Saf. 89:130-136;
		Biphasic = print("Biphasic: y ~ Alpha - Alpha / (1 + 10^((x - Beta) * Gamma)) + (1 - Alpha) / (1 + 10^((Delta - x) * Epsilon))"),
		
		sigmoid = print(c("Hill", "Hill_two", "Hill_three", "Hill_four", "Weibull", "Weibull_three", "Weibull_four", "Logit", "Logit_three",
			"Logit_four", "BCW(Box-Cox-Weibull)", "BCL(Box-Cox-Logit)", "GL(Generalized Logit)")),
		hormesis = 	print(c("Brain_Consens", "BCV", "Biphasic", "Hill_six"))
	)
}
