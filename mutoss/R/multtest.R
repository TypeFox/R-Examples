#######################################################################
########### here are the multtest methods - marginal siblings can be found in marginal.R ######################


### one sample model with multtest ###

onesamp.multtest <- function(data, alternative, robust, psi0, alpha, nulldist, B=1000, method, seed=12345) {
	result <- MTP(X=data, W = NULL, Y = NULL, Z = NULL, Z.incl = NULL, Z.test = NULL, 
			na.rm = TRUE, test = "t.onesamp", robust = robust, 
			standardize = TRUE, alternative = alternative, psi0 = psi0, 
			typeone = "fwer", k = 0, q = 0.1, fdr.method = "restricted", 
			alpha = alpha, smooth.null = FALSE, nulldist = nulldist, 
			B = B, ic.quant.trans = FALSE, MVN.method = "mvrnorm", 
			penalty = 1e-06, method = method, get.cr = FALSE, get.cutoff = FALSE, 
			get.adjp = TRUE, keep.nulldist = FALSE, keep.rawdist = FALSE, 
			seed = seed, cluster = 1, type = NULL, dispatch = NULL, marg.null = NULL, 
			marg.par = NULL, keep.margpar = TRUE, ncp = NULL, perm.mat = NULL, 
			keep.index = FALSE, keep.label = FALSE)
	return(list(adjPValues=result@adjp, rejected=as.vector(result@reject)))
}

mutoss.onesamp.multtest.model <- function(model) {
	return("typ" %in% names(model) && model$typ == "onesamp")
}

mutoss.onesamp.multtest <- function() { return(new(Class="MutossMethod",
					label="Resampling-based one sample test",
					errorControl="FWER",
					callFunction="onesamp.multtest",
					output=c("adjPValues", "rejected"),
					info="<h2>Resampling-based one sample test</h2>
							<p>There are different choices of resampling methods available for estimating the joint test statistics null distribution:\n\
							<ul>
							<li>\"boot.cs\": non-parametric bootstrap with centering and scaling</li>\n\
							<li>\"boot.ctr\": centered-only bootstrap distribution</li>\n\
							<li>\"boot.qt\": quantile transformed bootstrap distribution. the default marginal t-distribution with n-1 degree of freedom is used.</li>\n\
							<li>\"perm\": permutation distribution (refering to the Westfall and Young procedures)</li>\n\
							<li>\"ic\": under GUI construction (available at (library(multtest)) </li>\n
							</ul>
							</p> 
							<p>There are four adjustment methods to control the FWER:
							<ul>
							<li>\"sd.minP\": step-down common-quantile procedure based on the minima of unadjusted p-values</li>
							<li>\"sd.maxT\": step-down common-cut-off procedure based on the maxima of test statistics</li>
							<li>\"ss.minP\": single-step common-quantile procedure</li>
							<li>\"ss.maxT\": single-step common-cut-off procedure</li>
							</ul>
							</p>
							<p> The default number of bootstrap iterations (or number of permutations if resampling method is \"perm\") is 1000. This can be reduced to increase the speed
							of computations, at a cost to precision. However, it is recommended to use a large number of resampling iterations, e.g. 10,000. </p>
							<p> The robust version uses the Wilcoxon-Mann-Whitney test, otherwise a t-test will be performed.</p>
							<h3>Reference:</h3>
							<ul>
							<li>Dudoit, S. and van der Laan, M.J. (2007). <i>Mulitple Testing Procedures and Applications to Genomics.</i> Springer Series in Statistics.</li>\n\
							<li>Westfall, P.H. and Young, S.S. (1993). <i>Resampling-Based Multiple Testing. Examples and Methods for p-value adjustment. </i>Wiley Series in Probability and Mathematical Statistics. </li>\n   
							</ul>",
					parameters=list(
							data=list(type="ANY"),
							model=list(type="ANY"),
							robust=list(type="logical", label="Robust statistic"),
							alternative=list(type="character", label="Alternative", choices=c("two.sided", "less", "greater")),
							psi0=list(type="numeric", label="Hypothesized null value", default=0),
							alpha=list(type="numeric"),
							nulldist=list(type="character", label="Resampling Method", choices=c("boot.cs", "boot.ctr", "boot.qt", "perm")),
							B=list(type="numeric", label="Number of Resampling Iterations", default=1000),
							method=list(type="character", label="Adjustment Method", choices=c("sd.minP","sd.maxT","ss.minP","ss.maxT"),
									seed=list(type="ANY", default=12345)
							)
					))) }


### paired sample model with multtest ###


paired.multtest <- function(data, model, alternative, robust, psi0, alpha, nulldist, B=1000, method, seed=12345) {
	result <- MTP(X=data, W = NULL, Y = model$classlabel, Z = NULL, Z.incl = NULL, Z.test = NULL, 
			na.rm = TRUE, test = "t.pair", robust = robust, 
			standardize = TRUE, alternative = alternative, psi0 = psi0, 
			typeone = "fwer", k = 0, q = 0.1, fdr.method = "restricted", 
			alpha = alpha, smooth.null = FALSE, nulldist = nulldist, 
			B = B, ic.quant.trans = FALSE, MVN.method = "mvrnorm", 
			penalty = 1e-06, method = method, get.cr = FALSE, get.cutoff = FALSE, 
			get.adjp = TRUE, keep.nulldist = FALSE, keep.rawdist = FALSE, 
			seed = seed, cluster = 1, type = NULL, dispatch = NULL, marg.null = NULL, 
			marg.par = NULL, keep.margpar = TRUE, ncp = NULL, perm.mat = NULL, 
			keep.index = FALSE, keep.label = FALSE)
	return(list(adjPValues=result@adjp, rejected=as.vector(result@reject)))
}

mutoss.paired.multtest.model <- function(model) {
	return("typ" %in% names(model) && model$typ == "pairedsamp")
}

mutoss.paired.multtest <- function() { return(new(Class="MutossMethod",
					label="Resampling-based paired sample test",
					errorControl="FWER",
					callFunction="paired.multtest",
					output=c("adjPValues", "rejected"),
					info="<h2>Resampling-based paired test</h2>
							<p>There are different choices of resampling methods available for estimating the joint test statistics null distribution:\n\
							<ul>
							<li>\"boot.cs\": non-parametric bootstrap with centering and scaling</li>\n\
							<li>\"boot.ctr\": centered-only bootstrap distribution</li>\n\
							<li>\"boot.qt\": quantile transformed bootstrap distribution. the default marginal t-distribution with n-1 degree 
							of freedom is used, where n is the number of pairs.</li>\n\
							<li>\"perm\": permutation distribution (refering to the Westfall and Young procedures)</li>\n\
							<li>\"ic\": under GUI construction (available at (library(multtest)) </li>\n
							</ul>
							</p> 
							<p>There are four adjustment methods to control the FWER:
							<ul>
							<li>\"sd.minP\": step-down common-quantile procedure based on the minima of unadjusted p-values</li>
							<li>\"sd.maxT\": step-down common-cut-off procedure based on the maxima of test statistics</li>
							<li>\"ss.minP\": single-step common-quantile procedure</li>
							<li>\"ss.maxT\": single-step common-cut-off procedure</li>
							</ul>
							</p>
							<p> The default number of bootstrap iterations (or number of permutations if resampling method is \"perm\") is 1000. This can be reduced to increase the speed
							of computations, at a cost to precision. However, it is recommended to use a large number of resampling iterations, e.g. 10,000. </p>
							<p> The robust version uses the Wilcoxon-Mann-Whitney test, otherwise a t-test will be performed.</p>
							<h3>Reference:</h3>
							<ul>
							<li>Dudoit, S. and van der Laan, M.J. (2007). <i>Mulitple Testing Procedures and Applications to Genomics.</i> Springer Series in Statistics.</li>\n\
							<li>Westfall, P.H. and Young, S.S. (1993). <i>Resampling-Based Multiple Testing. Examples and Methods for p-value adjustment. </i>Wiley Series in Probability and Mathematical Statistics. </li>\n   
							</ul>",
					parameters=list(
							data=list(type="ANY"),
							model=list(type="ANY"),
							robust=list(type="logical", label="Robust statistic"),
							alternative=list(type="character", label="Alternative", choices=c("two.sided", "less", "greater")),
							psi0=list(type="numeric", label="Hypothesized null value", default=0),
							alpha=list(type="numeric"),
							nulldist=list(type="character", label="Resampling Method", choices=c("boot.cs", "boot.ctr", "boot.qt", "perm")),
							B=list(type="numeric", label="Number of Resampling Iterations", default=1000),
							method=list(type="character", label="Adjustment Method", choices=c("sd.minP","sd.maxT","ss.minP","ss.maxT"),
									seed=list(type="ANY", default=12345))
					)
			)) }

			
### two sample model with multtest ####			
			
twosamp.multtest <- function(data, model, alternative, robust, psi0, equalvar, alpha, nulldist, B=1000, method, seed=12345) {
	if (equalvar) {
		result <- MTP(X=data, W = NULL, Y = model$classlabel, Z = NULL, Z.incl = NULL, Z.test = NULL, 
				na.rm = TRUE, test = "t.twosamp.equalvar", robust = robust, 
				standardize = TRUE, alternative = alternative, psi0 = psi0, 
				typeone = "fwer", k = 0, q = 0.1, fdr.method = "restricted", 
				alpha = alpha, smooth.null = FALSE, nulldist = nulldist, 
				B = B, ic.quant.trans = FALSE, MVN.method = "mvrnorm", 
				penalty = 1e-06, method = method, get.cr = FALSE, get.cutoff = FALSE, 
				get.adjp = TRUE, keep.nulldist = FALSE, keep.rawdist = FALSE, 
				seed = seed, cluster = 1, type = NULL, dispatch = NULL, marg.null = NULL, 
				marg.par = NULL, keep.margpar = TRUE, ncp = NULL, perm.mat = NULL, 
				keep.index = FALSE, keep.label = FALSE)
	}
	else {
		result <- MTP(X=data, W = NULL, Y = model$classlabel, Z = NULL, Z.incl = NULL, Z.test = NULL, 
				na.rm = TRUE, test = "t.twosamp.unequalvar", robust = robust, 
				standardize = TRUE, alternative = alternative, psi0 = psi0, 
				typeone = "fwer", k = 0, q = 0.1, fdr.method = "restricted", 
				alpha = alpha, smooth.null = FALSE, nulldist = nulldist, 
				B = B, ic.quant.trans = FALSE, MVN.method = "mvrnorm", 
				penalty = 1e-06, method = method, get.cr = FALSE, get.cutoff = FALSE, 
				get.adjp = TRUE, keep.nulldist = FALSE, keep.rawdist = FALSE, 
				seed = seed, cluster = 1, type = NULL, dispatch = NULL, marg.null = NULL, 
				marg.par = NULL, keep.margpar = TRUE, ncp = NULL, perm.mat = NULL, 
				keep.index = FALSE, keep.label = FALSE)
	}
	return(list(adjPValues=result@adjp, rejected=as.vector(result@reject)))
}

mutoss.twosamp.multtest.model <- function(model) {
	return("typ" %in% names(model) && model$typ == "twosamp")
}

mutoss.twosamp.multtest <- function() { return(new(Class="MutossMethod",
					label="Resampling-based two sample test",
					errorControl="FWER",
					callFunction="twosamp.multtest",
					output=c("adjPValues", "rejected"),
					info="<h2>Resampling-based two sample test</h2>
							<p>There are different choices of resampling methods available for estimating the joint test statistics null distribution:\n\
							<ul>
							<li>\"boot.cs\": non-parametric bootstrap with centering and scaling</li>\n\
							<li>\"boot.ctr\": centered-only bootstrap distribution</li>\n\
							<li>\"boot.qt\": quantile transformed bootstrap distribution. the default marginal t-distribution with n-1 degree of freedom is used.</li>\n\
							<li>\"perm\": permutation distribution (refering to the Westfall and Young procedures)</li>\n\
							<li>\"ic\": under GUI construction (available at (library(multtest)) </li>\n
							</ul>
							</p> 
							<p>There are four adjustment methods to control the FWER:
							<ul>
							<li>\"sd.minP\": step-down common-quantile procedure based on the minima of unadjusted p-values</li>
							<li>\"sd.maxT\": step-down common-cut-off procedure based on the maxima of test statistics</li>
							<li>\"ss.minP\": single-step common-quantile procedure</li>
							<li>\"ss.maxT\": single-step common-cut-off procedure</li>
							</ul>
							</p>
							<p> The default number of bootstrap iterations (or number of permutations if resampling method is \"perm\") is 1000. This can be reduced to increase the speed
							of computations, at a cost to precision. However, it is recommended to use a large number of resampling iterations, e.g. 10,000. </p>
							<p> The robust version uses the Wilcoxon-Mann-Whitney test, otherwise a t-test will be performed.</p>
							<h3>Reference:</h3>
							<ul>
							<li>Dudoit, S. and van der Laan, M.J. (2007). <i>Mulitple Testing Procedures and Applications to Genomics.</i> Springer Series in Statistics.</li>\n\
							<li>Westfall, P.H. and Young, S.S. (1993). <i>Resampling-Based Multiple Testing. Examples and Methods for p-value adjustment. </i>Wiley Series in Probability and Mathematical Statistics. </li>\n   
							</ul>",
					parameters=list(
							data=list(type="ANY"),
							model=list(type="ANY"),
							robust=list(type="logical", label="Robust statistic"),
							alternative=list(type="character", label="Alternative", choices=c("two.sided", "less", "greater")),
							psi0=list(type="numeric", label="Hypothesized null value", default=0),
							equalvar=list(type="logical", label="Equal variance"),
							alpha=list(type="numeric"),
							nulldist=list(type="character", label="Resampling Method", choices=c("boot.cs", "boot.ctr", "boot.qt", "perm")),
							B=list(type="numeric", label="Number of Resampling Iterations", default=1000),
							method=list(type="character", label="Adjustment Method", choices=c("sd.minP","sd.maxT","ss.minP","ss.maxT"),
									seed=list(type="ANY", default=12345))
					)
			)) }



##################### F test #############################################
###########################################################################


ftest.multtest <- function(data, model, robust, alpha, nulldist, B=1000, method, seed=12345) {
	result <- MTP(X=data, W = NULL, Y = model$classlabel, Z = NULL, Z.incl = NULL, Z.test = NULL, 
			na.rm = TRUE, test = "f", robust = robust, 
			standardize = TRUE, typeone = "fwer", 
			alpha = alpha, smooth.null = FALSE, nulldist = nulldist, 
			B = B, ic.quant.trans = FALSE, MVN.method = "mvrnorm", 
			penalty = 1e-06, method = method, get.cr = FALSE, get.cutoff = FALSE, 
			get.adjp = TRUE, keep.nulldist = FALSE, keep.rawdist = FALSE, 
			seed = seed, cluster = 1, type = NULL, dispatch = NULL, marg.null = NULL, 
			marg.par = NULL, keep.margpar = TRUE, ncp = NULL, perm.mat = NULL, 
			keep.index = FALSE, keep.label = FALSE)
	return(list(adjPValues=result@adjp, rejected=as.vector(result@reject)))
}


mutoss.ftest.multtest <- function() { return(new(Class="MutossMethod",
					label="Resampling-based F test",
					errorControl="FWER",
					callFunction="ftest.multtest",
					output=c("adjPValues", "rejected"),
					info="<h2>Resampling-based F test</h2>
							<p>There are different choices of resampling methods available for estimating the joint test statistics null distribution:\n\
							<ul>
							<li>\"boot.cs\": non-parametric bootstrap with centering and scaling</li>\n\
							<li>\"boot.ctr\": centered-only bootstrap distribution</li>\n\
							<li>\"boot.qt\": quantile transformed bootstrap distribution. the default marginal F-distribution with df1=k-1, df2=n-k for k gruops is used.</li>\n\
							<li>\"perm\": permutation distribution (refering to the Westfall and Young procedures)</li>\n\
							<li>\"ic\": under GUI construction (available at (library(multtest)) </li>\n
							</ul>
							</p> 
							<p>There are four adjustment methods to control the FWER:
							<ul>
							<li>\"sd.minP\": step-down common-quantile procedure based on the minima of unadjusted p-values</li>
							<li>\"sd.maxT\": step-down common-cut-off procedure based on the maxima of test statistics</li>
							<li>\"ss.minP\": single-step common-quantile procedure</li>
							<li>\"ss.maxT\": single-step common-cut-off procedure</li>
							</ul>
							</p>
							<p> The default number of bootstrap iterations (or number of permutations if resampling method is \"perm\") is 1000. This can be reduced to increase the speed
							of computations, at a cost to precision. However, it is recommended to use a large number of resampling iterations, e.g. 10,000. </p>
							<p> The robust version uses the Kruskal-Wallis test, otherwise a F test will be performed.</p>
							
							<h3>Reference:</h3>
							<ul>
							<li>Dudoit, S. and van der Laan, M.J. (2007). <i>Mulitple Testing Procedures and Applications to Genomics.</i> Springer Series in Statistics.</li>\n\
							<li>Westfall, P.H. and Young, S.S. (1993). <i>Resampling-Based Multiple Testing. Examples and Methods for p-value adjustment. </i>Wiley Series in Probability and Mathematical Statistics. </li>\n   
							</ul>",
					parameters=list(
							data=list(type="ANY"),
							model=list(type="ANY"),
							robust=list(type="logical", label="Robust statistic"),
							alpha=list(type="numeric"),
							nulldist=list(type="character", label="Resampling Method", choices=c("boot.cs", "boot.ctr", "boot.qt", "perm")),
							B=list(type="numeric", label="Number of Resampling Iterations", default=1000),
							method=list(type="character", label="Adjustment Method", choices=c("sd.minP","sd.maxT","ss.minP","ss.maxT"),
									seed=list(type="ANY", default=12345)
							)))
	)}

mutoss.ftest.multtest.model <- function(model) {
	return("typ" %in% names(model) && model$typ == "ftest")
}