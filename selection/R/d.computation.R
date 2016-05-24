##' Compute Cohen's d
##'
##' Compute Cohen's d
##'	
##' Cohen's d is used to assess average group differences, regardless of sample size. It is consequently a measure of effect size. 
##' @param y The y variable (quantitative) on which group differences are to be assessed.
##' @param group The grouping variable
##' @param short Logical. If true, it will only output the value of d. Otherwise, it also outputs the mean differences, the pooled variance, and
##' the pair of means. 
##' @param target.group The group to reference.
##' @return If \code{short=FALSE}, the mean differences, the pooled variance, and
##' the pair of means. Otherwise, just the d.
##' @author Dustin Fife
##' @export
##' @examples
##' data(iris)
##' d.computation(iris$Petal.Width, iris$Species, short=FALSE)
d.computation = function(y,group, short=TRUE, target.group=NULL){
	levs = unique(group)
	if (!is.null(target.group)){
		g1 = which(as.character(group)==target.group);g2 = which(as.character(group)!=target.group)		
	} else {
		g1 = which(as.character(group)==levs[1]);g2 = which(as.character(group)!=levs[2])		
	}
	n1 = length(g1); n2=length(g2)
	var1 = var(y[g1]);var2=var(y[g2])
	sp = sqrt(((n1-1)*var1 + (n2-1)*var2)/(n1+n2-2))
	xbar1 = mean(y[g1]); xbar2=mean(y[g2])
	d = (xbar1-xbar2)/sp
	
	if (short){
		d
	} else {
		list (d=d, delta=xbar1-xbar2, sigma=sp, xbar1=xbar1, xbar2=xbar2)
	}
}