"ci.reg.coef" <-
function(b.j, SE.b.j=NULL, s.Y=NULL, s.X=NULL, N, p, R2.Y_X=NULL, R2.j_X.without.j=NULL, conf.level=.95, R2.Y_X.without.j=NULL, t.value=NULL, alpha.lower=NULL, alpha.upper=NULL, Noncentral=FALSE, Suppress.Statement=FALSE, ...)
{

# Determine if NC was used instead of Noncentral
#tmp <- try(is.null(NC), silent=TRUE)
#if(tmp==TRUE | tmp==FALSE) Noncentral <- NC

if(!is.null(t.value))
{
obs.t <- t.value
}

if(!is.null(b.j) & !is.null(SE.b.j))
{
obs.t <- b.j/SE.b.j
}

if(is.null(SE.b.j))
{
{
if(is.null(t.value)) 
{
if(is.null(s.Y) & is.null(s.X) & Noncentral==TRUE) 
    {
    s.Y <- 1
    s.X <- 1
    }

if(is.null(s.Y) | is.null(s.X)) stop("You need to specify 's.Y' and 's.X'.")

SE.b.j <- sqrt(((1-R2.Y_X)/((1-R2.j_X.without.j)*(N-p-1))))*(s.Y/s.X)
}
if(!is.null(t.value) & !is.null(b.j)) SE.b.j <- b.j/t.value
if(!is.null(t.value) & is.null(b.j))
{
SE.b.j <- sqrt(((1-R2.Y_X)/((1-R2.j_X.without.j)*(N-p-1))))*(s.Y/s.X)
b.j <- obs.t*SE.b.j 
}
}
if(is.null(t.value)) obs.t <- b.j/SE.b.j
}


if(is.null(b.j))
{
b.j <- (((R2.Y_X - R2.Y_X.without.j)/(1-R2.j_X.without.j))^.5)*(s.Y/s.X)
if(is.null(SE.b.j))
{
SE.b.j <- sqrt(((1-R2.Y_X)/((1-R2.j_X.without.j)*(N-p-1))))*(s.Y/s.X)
}
obs.t <- b.j/SE.b.j
}


if(!is.null(conf.level))
{
if(conf.level >= 1 | conf.level <=0) stop("You have not properly specified \'conf.level\'", call. = FALSE)

alpha.lower <- alpha.upper <- (1-conf.level)/2
}
if(is.null(conf.level))
{
if(is.null(alpha.lower) & is.null(alpha.upper)) stop("You need to specify either \'conf.level\', or \'alpha.lower\' and \'alpha.upper\'.", call.=FALSE)
if(alpha.lower > .5 | alpha.lower <0) stop("You have not properly specified \'alpha.lower\' correctly.", call. = FALSE)
if(alpha.upper > .5 | alpha.upper <0) stop("You have not properly specified \'alpha.upper\' correctly.", call. = FALSE)
}

if(Noncentral==FALSE)
{
if(Suppress.Statement!=TRUE) print(paste((1-alpha.lower-alpha.upper)*100," percent CI limits (with corresponding probability) around the jth population regression coefficient calculated using the (central) t-distribution with ", N-p-1, " degrees of freedom follow.", sep=""))
return(list(Lower.Limit.for.beta.j=b.j + qt(alpha.lower, df=N-p-1)*SE.b.j, Prob.Less.Lower=alpha.lower, Upper.Limit.for.beta.j=b.j + qt(1-alpha.upper, df=N-p-1)*SE.b.j, Prob.Greater.Upper=alpha.upper))
}

if(Noncentral==TRUE)
{
NC.t.values <- conf.limits.nct(ncp=obs.t, df=N-p-1, conf.level=NULL, alpha.lower=alpha.lower, alpha.upper=alpha.upper)
if(Suppress.Statement!=TRUE) print(paste((1-alpha.lower-alpha.upper)*100," percent CI limits (with corresponding probability) around the jth population regression coefficient calculated using the noncentral t-distribution with ", N-p-1, " degrees of freedom and noncentrality parameter ", round(obs.t, 4), " follow.", sep=""))
return(list(Lower.Limit.for.beta.j=NC.t.values$Lower*SE.b.j, Prob.Less.Lower=NC.t.values$Prob.Less.Lower, Upper.Limit.for.beta.j=NC.t.values$Upper.Limit*SE.b.j, Prob.Greater.Upper=NC.t.values$Prob.Greater.Upper))
}

}
