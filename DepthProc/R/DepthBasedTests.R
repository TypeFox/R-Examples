
################## Multivariate Wilcoxon Test #########################
## TODO:
## - new object (S3 for compability with other tests)
## - more examples

#' @title Depth based multivariate Wilcoxon test for a scale difference.
#' @export
#' 
#' @param x data matrix
#' @param y data matrix
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param ... arguments passed to depth function(e.g. method)
#' 
#' @details
#' 
#' Having two samples  \eqn{ {X}^{n} }  and  \eqn{ {Y}^{m} }  using any depth function, we can compute depth values in a combined sample  \eqn{ {Z}^{n+m} }  =  \eqn{ {X}^{n}\cup {Y}^{m} } , assuming the empirical distribution calculated basing on all observations, or only on observations belonging to one of the samples  \eqn{ {X}^{n} }  or  \eqn{ {Y}^{m}. }   
#' 
#' For example if we observe  \eqn{ {X}_{l}'s }  depths are more likely to cluster tightly around the center of the combined sample, while  \eqn{ {Y}_{l}'s }  depths are more likely to scatter outlying positions, then we conclude  \eqn{ {Y}^{m} }  was drawn from a distribution with larger scale.  
#' 
#' Properties of the DD plot based statistics in the i.i.d setting were studied by Li \& Liu (2004). Authors proposed several DD-plot based statistics and presented bootstrap arguments for their consistency and good effectiveness in comparison to Hotelling  \eqn{ T^2 }  and multivariate analogues of Ansari-Bradley and Tukey-Siegel statistics. Asymptotic distributions of depth based multivariate Wilcoxon rank-sum test statistic under the null and general alternative hypotheses were obtained by Zuo \& He (2006). Several properties of the depth based rang test involving its unbiasedness was critically discussed by Jureckova \& Kalina (2012).  Basing on DD-plot object, which is available within the \pkg{DepthProc} it is possible to define several multivariate generalizations of one-dimensional rank and order statistics in an easy way. These generalizations cover well known {Wilcoxon rang-sum statistic}.
#' 
#'   The depth based multivariate Wilcoxon rang sum test is especially useful for the multivariate scale changes detection and was introduced among other by Liu \& Singh (2003) and intensively studied by Jureckowa \& Kalina (2012) and Zuo \& He (2006) in the i.i.d. setting.
#' 
#'   For the samples  \eqn{ {{{X}}^{m}}=\{{{{X}}_{1}},...,{{{X}}_{m}}\} }  ,  \eqn{ {{{Y}}^{n}}=\{{{{Y}}_{1}},...,{{{Y}}_{n}}\} } , their  \eqn{ d_{1}^{X},...,d_{m}^{X} }  ,  \eqn{ d_{1}^{Y},...,d_{n}^{Y} } ,  depths w.r.t. a combined sample  \eqn{ {{Z}}={{{X}}^{n}}\cup {{{Y}}^{m}} }  the Wilcoxon statistic is defined as  \eqn{ S=\sum\limits_{i=1}^{m}{{{R}_{i}}}}, where  \eqn{ {R}_{i} }  denotes the rang of the i-th observation,  \eqn{ i=1,...,m }  in the combined sample  \eqn{ R({{{y}}_{l}})=  \#\left\{ {{{z}}_{j}}\in {{{Z}}}:D({{{z}}_{j}},{{Z}})\le D({{{y}}_{l}},{{Z}}) \right\}, l=1,...,m. }  
#' 
#' The distribution of  \eqn{ S }  is symmetric about   \eqn{ E(S)=1/2m(m{+}n{+1)} } , its variance is   \eqn{ {{D}^{2}}(S)={1}/{12}\;mn(m+n+1)}.
#' 
#' @references
#' 
#' Jureckova J, Kalina J (2012). Nonparametric multivariate rank tests and their unbiasedness. Bernoulli, 18(1), 229-251.
#' Li J, Liu RY (2004). New nonparametric tests of multivariate locations and scales using data depth. Statistical Science, 19(4), 686-696.
#' Liu RY, Singh K (1995). A quality index based on data depth and multivariate rank tests. Journal of American Statistical Association, 88, 252-260.
#' Zuo Y, He X (2006). On the limiting distributions of multivariate depth-based rank sum statistics and related tests. The Annals of Statistics, 34, 2879-2896.
#' 
#' @examples
#' 
#' 
#' x = mvrnorm(100, c(0,0), diag(2))
#' y = mvrnorm(100, c(0,0), diag(2)*1.4)
#' mWilcoxonTest(x,y)
#' mWilcoxonTest(x,y, method = "LP")
#' 
#' #EXAMPLE 2
#' data(under5.mort)
#' data(inf.mort)
#' data(maesles.imm)
#' data2011=na.omit(cbind(under5.mort[,22],inf.mort[,22],maesles.imm[,22]))
#' data1990=na.omit(cbind(under5.mort[,1],inf.mort[,1],maesles.imm[,1]))
#' mWilcoxonTest(data2011,data1990)
#' 
#' 
mWilcoxonTest = function(x, y, alternative = "two.sided", ...)
{
  total = rbind(x, y)
  dep_x  = depth(x,total, ...)
  dep_y  = depth(y,total, ...)
  test_res = wilcox.test(dep_x,dep_y, alternative = alternative)
  
  
  test_res$null.value = 1
  test_res$method = "Multivariate Wilcoxon test for equality of dispersion"
  force(names(test_res$null.value) <- "dispersion ratio")
  test_res
}


################

