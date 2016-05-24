# TODO: Add comment
# 
# Author: wiebke
###############################################################################


#### one sample model ####

onesamp.model <- function() {
	return(list(model=list(typ="onesamp")))
}

mutoss.onesamp.model <- function() { return(new(Class="MutossMethod",
					label="One-sample test",
					callFunction="onesamp.model",
					output=c("model"),
					info="<h2>One sample model</h2>
							<p>The input for this one sample model is a data matrix whose columns represent the samples and the rows represent the multiple endpoints. 
							E.g. for genomics this would be a gene matrix, where each row gives the expression for a single gene.</p> 
							<p> For the next steps you have the following choices:</p>
							<p> Either marginal hypotheses tests (if robust Wilcoxon otherwise t-test) could be performed on each row of the 
							data matrix to obtain raw p-values which then need to be adjusted for multiplicity to control a chosen error rate.</p>
							<p> Or resampling based methods could be performed based on Dudoit and van der Laan (2007) to obtain adjusted p-values which control the 
							FWER. Afterwards it is possible to use augmentation procedures to get adjusted p-values for control of FDR, FDX or gFWER.</p> 
							<h3>Reference:</h3>
							<ul>
							<li>Dudoit, S. and van der Laan, M.J. (2007). <i>Mulitple Testing Procedures and Applications to Genomics.</i> Springer Series in Statistics.</li>
							</ul>",
					parameters=list(
					)
			)) }

onesamp.marginal <- function(data, robust, alternative, psi0) {
	result <- NULL
	if (robust) {
		result <- apply(data, 1, function(x) {wilcox.test(x, alternative=alternative, mu=psi0)$p.value} )
	} else {
		result <- apply(data, 1, function(x) {t.test(x ,alternative=alternative, mu=psi0)$p.value} )
	}
	return(list(pValues=result))
}

mutoss.onesamp.marginal.model <- function(model) {
	return("typ" %in% names(model) && model$typ == "onesamp")
}


mutoss.onesamp.marginal <- function() { return(new(Class="MutossMethod",
					label="One-sample test",
					callFunction="onesamp.marginal",
					output=c("pValues"),
					info="<h2>Marginal one sample test.</h2>
							<p>The robust version uses the Wilcoxon-Mann-Whitney test, otherwise a t-test will be performed.</p> 
							<h3>Reference:</h3>
							<ul>
							<li>Wilcoxon, F. (1945). <i>Individual Comparisons by Ranking Methods.</i> Biometrics Bulletin 1:80-83.</li>\n\
							<li>Mann, H. and Whitney, D. (1947). <i>On a test of whether one of two random variables is stochastically larger 
							than the other.</i> Annals of Mathematical Statistics 18:50-60</li>\n\
							<li>Student (1908).<i>The probable error of a mean.</i> Biometrika, 6(1):1-25.</li>\n 
							</ul>",
					parameters=list(
							data=list(type="ANY"),
							robust=list(type="logical", label="Robust statistic"),
							alternative=list(type="character", label="Alternative", choices=c("two.sided", "less", "greater")),
							psi0=list(type="numeric", label="Hypothesized null value", default=0)
					)
			)) }



#### two sample model ####


twosamp.model <- function(classlabel) {
	classlabel <- as.vector(classlabel)
	return(list(model=list(typ="twosamp", classlabel=classlabel)))
}

mutoss.twosamp.model <- function() { return(new(Class="MutossMethod",
					label="Two sample test",
					callFunction="twosamp.model",
					output=c("model"),
					info="<h2>Two sample model</h2>
							<p>The input for this one sample model is a data matrix whose columns represent the samples and the rows represent the multiple endpoints. 
							E.g. for genomics this would be a gene matrix, where each row gives the expression for a single gene.</p> 
							<p> Furthermore, a classlabel needs to be provided to distinguish the two sample groups.
							<p> For the next steps you have the following choices:</p>
							<p> Either marginal hypotheses tests (if robust Wilcoxon otherwise t-test) could be performed on each row of the 
							data matrix to obtain raw p-values which then need to be adjusted for multiplicity to control a chosen error rate.</p>
							<p> Or resampling based methods could be performed based on Dudoit and van der Laan (2007) to obtain adjusted p-values which control the 
							FWER. Afterwards it is possible to use augmentation procedures to get adjusted p-values for control of FDR, FDX or gFWER.</p> 
							<h3>Reference:</h3>
							<ul>
							<li>Dudoit, S. and van der Laan, M.J. (2007). <i>Mulitple Testing Procedures and Applications to Genomics.</i> Springer Series in Statistics.</li>
							</ul>",
					parameters=list(
							classlabel=list(type="RObject", label="classlabel")
					)
			)) }

twosamp.marginal <- function(data, model, robust, alternative, psi0, equalvar) {
	label <- as.numeric(as.factor(model$classlabel))
	result <- NULL
	if (robust) {
		result <- apply(data, 1, function(x) {wilcox.test(x=x[ ,label==1], y=x[label==2], alternative=alternative, mu=psi0)$p.value} )
	} else {
		result <- apply(data, 1, function(x) {t.test(x=x[ ,label==1], y=x[label==2], alternative=alternative, mu=psi0, equal.var=equalvar)$p.value} )
	}
	return(list(pValues=result))
}

mutoss.twosamp.marginal.model <- function(model) {
	return("typ" %in% names(model) && model$typ == "twosamp")
}

mutoss.twosamp.marginal <- function() { return(new(Class="MutossMethod",
					label="Two sample test",
					callFunction="twosamp.marginal",
					output=c("pValues"),
					info="<h2></h2>
							<p>The robust version uses the Wilcoxon-Mann-Whitney test, otherwise a two-sample t-test will be performed.</p> 
							<h3>Reference:</h3>
							<ul>
							<li>Wilcoxon, F. (1945). <i>Individual Comparisons by Ranking Methods.</i> Biometrics Bulletin 1:80-83.</li>\n\
							<li>Mann, H. and Whitney, D. (1947). <i>On a test of whether one of two random variables is stochastically larger 
							than the other.</i> Annals of Mathematical Statistics 18:50-60</li>\n
							</ul>",
					parameters=list(
							data=list(type="ANY"),
							model=list(type="ANY"),
							robust=list(type="logical", label="Robust statistic"),
							alternative=list(type="character", label="Alternative", choices=c("two.sided", "less", "greater")),
							psi0=list(type="numeric", label="Hypothesized null value", default=0),
							equalvar=list(type="logical", label="Equal variance")
					)
			)) }


### paired sample model ###

paired.model <- function(classlabel) {
	classlabel <- as.vector(classlabel)
	return(list(model=list(typ="pairedsamp", classlabel=classlabel)))
}

mutoss.paired.model <- function() { return(new(Class="MutossMethod",
					label="Paired sample test",
					callFunction="paired.model",
					output=c("model"),
					info="<h2>Paired sample test</h2>
							<p>The robust version uses the Wilcoxon signed rank test, otherwise a paired t-test will be performed.</p> 
							<p>The input for this paired sample model is a data matrix whose columns represent the samples and the rows represent the multiple endpoints. 
							E.g. for genomics this would be a gene matrix, where each row gives the expression for a single gene.</p> 
							<p> Furthermore, a classlabel needs to be provided to distinguish the two paired groups. The arrangement of group indices does not matter, as long
							as the columns are arranged in the same corresponding order between groups. For example, if group 1 is code as 0 and group 2 is 
							coded as 1, for 3 pairs of data, it does not matter if the classlabel is coded as (0,0,0,1,1,1) or (1,1,1,0,0,0) or (0,1,0,1,0,1)
							or (1,0,1,0,1,0), the paired differences between groups will be calculated as group2 - group1. 
							</p>
							<p> For the next steps you have the following choices:</p>
							<p> Either marginal hypotheses tests (if robust Wilcoxon otherwise t-test) could be performed on each row of the 
							data matrix to obtain raw p-values which then need to be adjusted for multiplicity to control a chosen error rate.</p>
							<p> Or resampling based methods could be performed based on Dudoit and van der Laan (2007) to obtain adjusted p-values which control the 
							FWER. Afterwards it is possible to use augmentation procedures to get adjusted p-values for control of FDR, FDX or gFWER.</p> 
							
							
							<h3>Reference:</h3>
							<ul>
							<li>Dudoit, S. and van der Laan, M.J. (2007). <i>Mulitple Testing Procedures and Applications to Genomics.</i> Springer Series in Statistics.</li>
							</ul>",
					parameters=list(
							classlabel=list(type="RObject", label="classlabel")
					)
			)) }

paired.marginal <- function(data, model, robust, alternative, psi0, equalvar) {
	label <- as.numeric(as.factor(model$classlabel))
	result <- NULL
	if (robust) {
		result <- apply(data, 1, function(x) {wilcox.test(x=x[label==1], y=x[label==2], alternative=alternative, mu=psi0, paired=TRUE, var.equal=equalvar)$p.value} )
	} else {
		result <- apply(data, 1, function(x) {t.test(x=x[label==1], y=x[label==2], alternative=alternative, mu=psi0, paired=TRUE, var.equal=equalvar)$p.value} )
	}
	return(list(pValues=result))
}

mutoss.paired.marginal <- function() { return(new(Class="MutossMethod",
					label="Paired sample test",
					callFunction="paired.marginal",
					output=c("pValues"),
					info="<h2></h2>
							<p>The robust version uses the Wilcoxon test, otherwise a paired t-test will be performed.</p> 
							<p>A vector of classlabels needs to be provided to distinguish the two paired groups. The arrangement of group indices does not matter, as long
							as the columns are arranged in the same corresponding order between groups. For example, if group 1 is code as 0 and group 2 is 
							coded as 1, for 3 pairs of data, it does not matter if the classlabel is coded as (0,0,0,1,1,1) or (1,1,1,0,0,0) or (0,1,0,1,0,1)
							or (1,0,1,0,1,0), the paired differences between groups will be calculated as group2 - group1. </p>
							<p>You could either choose a valid R object to load as classlabels or you could provide it manually by inserting e.g. c(0,1,0,1,0,1) or rep(c(0,1), each=5) or rep(c(0,1), 5). </p>",
					parameters=list(
							data=list(type="ANY"),
							model=list(type="ANY"),
							robust=list(type="logical", label="Robust statistic"),
							alternative=list(type="character", label="Alternative", choices=c("two.sided", "less", "greater")),
							psi0=list(type="numeric", label="Hypothesized null value", default=0),
							equalvar=list(type="logical", label="Equal variance")
					)
			)) }

mutoss.paired.marginal.model <- function(model) {
	return("typ" %in% names(model) && model$typ == "pairedsamp")
}

### f test model ###


ftest.model <- function(classlabel) {
	classlabel <- as.vector(classlabel)
	return(list(model=list(typ="ftest", classlabel=classlabel)))
}

mutoss.ftest.model <- function() { return(new(Class="MutossMethod",
					label="F test",
					callFunction="ftest.model",
					output=c("model"),
					info="<h2>F test</h2>
							<p>The input for this F test model is a data matrix whose columns represent the samples and the rows represent the multiple endpoints. 
							E.g. for genomics this would be a gene matrix, where each row gives the expression for a single gene.</p> 
							<p> Furthermore, a classlabel needs to be provided to distinguish k sample groups.</p>
							<p> For the next steps you have the following choices:</p>
							<p> Either marginal hypotheses tests (if robust Kruskal-Wallis test, otherwise F-test) could be performed on each row of the 
							data matrix to obtain raw p-values which then need to be adjusted for multiplicity to control a chosen error rate.</p>
							<p> Or resampling based methods could be performed based on Dudoit and van der Laan (2007) to obtain adjusted p-values which control the 
							FWER. Afterwards it is possible to use augmentation procedures to get adjusted p-values for control of FDR, FDX or gFWER.</p> 
							<h3>Reference:</h3>
							<ul>
							<li>Dudoit, S. and van der Laan, M.J. (2007). <i>Mulitple Testing Procedures and Applications to Genomics.</i> Springer Series in Statistics.</li>
							</ul>",
					parameters=list(
							classlabel=list(type="RObject", label="classlabel")
					)
			)) }

ftest.marginal <- function(data, model, robust) {
	label <- as.numeric(as.factor(model$classlabel))
	result <- NULL
	if (robust) {
		result <- apply(data, 1, function(x) {kruskal.test(x=x, g=label)$p.value} )
	} else {
		result <- apply(data, 1, function(x) { out=x
					anova(lm( out ~ label ))$'Pr(>F)'[1]} )	
	}	
	return(list(pValues=result))
}

mutoss.ftest.marginal <- function() { return(new(Class="MutossMethod",
					label="F test",
					callFunction="ftest.marginal",
					output=c("pValues"),
					info="<h2></h2>
							<p>Robust = Kruskal-Wallis test. Otherwise F-test.</p> 
							<p></p>
							<h3>Reference:</h3>
							<ul>
							<li>Kruskal, W.H. und Wallis, W.A. (1952). <i>Use of ranks in one-criterion variance analysis.</i> JASA, 47:583-621</li>
							</ul>",
					parameters=list(
							data=list(type="ANY"),
							model=list(type="ANY"),
							robust=list(type="logical", label="Robust statistic")
					)
			)) }

mutoss.ftest.marginal.model <- function(model) {
	return("typ" %in% names(model) && model$typ == "ftest")
}
