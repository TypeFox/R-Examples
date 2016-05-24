

CIpvBI <- function(x1, x0, pr, conf.level=0.95,
 alternative=c("two.sided", "less", "greater"), B=5000, shapes1=c(1,1), shapes0=c(1,1), ...)
{

if(!is.integer(x0) & !is.numeric(x0)){stop("x0 must be an integer vector")}
if(length(x0)!=2){stop("x0 must be an integer vector with two elements")}
if(any(x0<0)){stop("elements in x0 must be non-negative")}


if(!is.integer(x1) & !is.numeric(x1)){stop("x1 must be an integer vector")}
if(length(x1)!=2){stop("x1 must be an integer vector with two elements")}
if(any(x1<0)){stop("elements in x1 must be non-negative")}


if(!is.numeric(pr) | length(pr)!=1){stop("pr must be an single numeric value")}
if(pr>=1 | pr<=0){stop("pr must be an single numeric value between 0 and 1")}

if(!is.numeric(conf.level) | length(conf.level)!=1){stop("conf.level must be an single numeric value")}
if(conf.level>=1 | conf.level<=0){stop("conf.level must be an single numeric value between 0 and 1")}

if(!is.numeric(shapes0) | length(shapes0)!=2){stop("shapes0 must be numeric vector")}
if(any(shapes0<=0)){stop("shapes0 must be a vector containing two positive numbers")}

if(!is.numeric(shapes1) | length(shapes1)!=2){stop("shapes1 must be numeric vector")}
if(any(shapes1<=0)){stop("shapes1 must be a vector containing two positive numbers")}

alternative <- match.arg(alternative)

# Equation (2) in Stamey and Holt, (2010), p. 103

Scdj <- rbeta(n=B, shape1 = x1[1] + shapes1[1], shape2 = x1[2] + shapes1[2]) 
Ccdj <- rbeta(n=B, shape1 = x0[2] + shapes0[2], shape2 = x0[1] + shapes0[1]) 

ppvj <- (Scdj*pr)/(Scdj*pr + (1-Ccdj)*(1-pr))
npvj <- (Ccdj*(1-pr))/((1-Scdj)*pr + (Ccdj)*(1-pr))

# corresponding estimates
estScd <- (x1[1] + shapes1[1])/(x1[1] + shapes1[1] + x1[2] + shapes1[2])
estCcd <- (x0[2] + shapes0[2])/(x0[2] + shapes0[2] + x0[1] + shapes0[1])

estppv <- (estScd*pr)/(estScd*pr + (1-estCcd)*(1-pr))
estnpv <- (estCcd*(1-pr))/((1-estScd)*pr + (estCcd)*(1-pr))


switch(alternative,
two.sided={
probs <- c((1-conf.level)/2, 1-(1-conf.level)/2)
cippv <- quantile(x=ppvj, probs=probs,...)
cinpv <- quantile(x=npvj, probs=probs,...)
},

less={
probs <-  conf.level
cippv <- c(0, quantile(x=ppvj, probs=probs,...))
cinpv <- c(0, quantile(x=npvj, probs=probs,...))

},
greater={
probs <- 1-conf.level
cippv <- c(quantile(x=ppvj, probs=probs, ...),1)
cinpv <- c(quantile(x=npvj, probs=probs, ...),1)
}
)

conf.int <- rbind("NPV"=cinpv, "PPV" = cippv)
colnames(conf.int)<-c("lower","upper")
estimate <- c(estnpv, estppv)
postmed <- c( median(npvj), median(ppvj))
names(estimate) <- names(postmed) <- c("NPV","PPV")

tab <- rbind( cbind(x1, x0), c(sum(x1), sum(x0)) )
colnames(tab) <- c("True positive", "True negative")
rownames(tab) <- c("Test positive", "Test negative", "Total")

prior <- cbind(shapes1, shapes0)
colnames(prior) <- c("pr.sens","pr.1-spec")
rownames(prior) <- c("a:shape1","b:shape2")

return(list(conf.int=conf.int, estimate=estimate, postmed=postmed, tab=tab, prior=prior))
}



CIpvBII <- function(x1, x0, xpr=NULL, conf.level=0.95,
 alternative=c("two.sided", "less", "greater"), B=5000, shapes1=c(1,1), shapes0=c(1,1), shapespr=c(1,1), ...)
{

if(!is.integer(x0) & !is.numeric(x0)){stop("x0 must be an integer vector")}
if(length(x0)!=2){stop("x0 must be an integer vector with two elements")}
if(any(x0<0)){stop("elements in x0 must be non-negative")}

if(!is.integer(x1) & !is.numeric(x1)){stop("x1 must be an integer vector")}
if(length(x1)!=2){stop("x1 must be an integer vector with two elements")}
if(any(x1<0)){stop("elements in x1 must be non-negative")}

if(!is.null(xpr)){
if(!is.integer(xpr) & !is.numeric(xpr)){stop("xpr must be an integer vector")}
if(length(xpr)!=2){stop("xpr must be an integer vector with two elements")}
if(any(xpr<0)){stop("elements in xpr must be non-negative")}}

if(!is.numeric(conf.level) | length(conf.level)!=1){stop("conf.level must be an single numeric value")}
if(conf.level>=1 | conf.level<=0){stop("conf.level must be an single numeric value between 0 and 1")}

if(!is.numeric(shapes0) | length(shapes0)!=2){stop("shapes0 must be numeric vector")}
if(any(shapes0<=0)){stop("shapes0 must be a vector containing two positive numbers")}

if(!is.numeric(shapes1) | length(shapes1)!=2){stop("shapes1 must be numeric vector")}
if(any(shapes1<=0)){stop("shapes1 must be a vector containing two positive numbers")}

if(!is.numeric(shapespr) | length(shapespr)!=2){stop("shapespr must be numeric vector")}
if(any(shapespr<=0)){stop("shapespr must be a vector containing two positive numbers")}

alternative<-match.arg(alternative)

# Stamey and Holt, (2010), equation (2)ff (p. 103-104)

Scdj <- rbeta(n=B, shape1 = x1[1] + shapes1[1], shape2 = x1[2] + shapes1[2]) 
Ccdj <- rbeta(n=B, shape1 = x0[2] + shapes0[2], shape2 = x0[1] + shapes0[1]) 

if(is.null(xpr)){pcdj <- rbeta(n=B, shape1 = shapespr[1], shape2=shapespr[2] )}else{ pcdj <- rbeta(n=B, shape1 = (xpr[1] + shapespr[1]), shape2 = (xpr[2] + shapespr[2]) )}

ppvj <- (Scdj*pcdj)/(Scdj*pcdj + (1-Ccdj)*(1-pcdj))

npvj <- (Ccdj*(1-pcdj))/((1-Scdj)*pcdj + (Ccdj)*(1-pcdj))


# corresponding estimates
estScd <- (x1[1] + shapes1[1])/(x1[1] + shapes1[1] + x1[2] + shapes1[2])
estCcd <- (x0[2] + shapes0[2])/(x0[2] + shapes0[2] + x0[1] + shapes0[1])

if(is.null(xpr)){estpr <- shapespr[1]/(shapespr[1]+shapespr[2])}else{estpr <- (xpr[1]+shapespr[1])/(xpr[1] + xpr[2] + shapespr[1] + shapespr[2])}

estppv <- (estScd*estpr)/(estScd*estpr + (1-estCcd)*(1-estpr))
estnpv <- (estCcd*(1-estpr))/((1-estScd)*estpr + (estCcd)*(1-estpr))


switch(alternative,
two.sided={
probs <- c((1-conf.level)/2, 1-(1-conf.level)/2)
cippv <- quantile(x=ppvj, probs=probs,...)
cinpv <- quantile(x=npvj, probs=probs,...)
},

less={
probs <-  conf.level
cippv <- c(0, quantile(x=ppvj, probs=probs,...))
cinpv <- c(0, quantile(x=npvj, probs=probs,...))

},
greater={
probs <- 1-conf.level
cippv <- c(quantile(x=ppvj, probs=probs, ...),1)
cinpv <- c(quantile(x=npvj, probs=probs, ...),1)
}
)


conf.int <- rbind("NPV"=cinpv, "PPV" = cippv)
colnames(conf.int)<-c("lower","upper")
estimate <- c(estnpv, estppv)
postmed <- c( median(npvj), median(ppvj))
names(estimate) <- names(postmed) <- c("NPV","PPV")

tab <- rbind( cbind(x1, x0), c(sum(x1), sum(x0)))
colnames(tab) <- c("True positive", "True negative")
rownames(tab) <- c("Test positive", "Test negative", "Total")

prior <- cbind(shapes1, shapes0)
colnames(prior) <- c("pr.sens","pr.1-spec")
rownames(prior) <- c("a:shape1","b:shape2")

if(is.null(xpr)){XPR <- c(0,0)}else{XPR <- xpr}
names(XPR) <- c("obs.pos","obs.neg")
names(shapespr) <- c("a:shape1","b:shape2")


return(list(conf.int=conf.int, estimate=estimate, postmed=postmed, tab=tab, prior=prior, prev.study=XPR, prev.prior=shapespr))
}
