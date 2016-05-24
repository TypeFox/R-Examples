Power3Cats <- function(kappa0, kappa1, props, raters=2, alpha=0.05, power=0.80)
{
#Error Checking
if ( (raters != 1) && (raters !=2) && (raters !=3) && (raters != 4) && (raters != 5) && (raters != 6) )
        stop("Sorry, this function is designed for between 2 to 6 raters.")

if (length(props) != 3)
	stop("Sorry, you must enter the anticipated proportions of each of the three categories.")

if ( abs( sum(props) - 1) >= 0.001 )
	stop("Sorry, the three proportions must sum to one.")

for (i in 1:3)
{
if ((props[i] >= 1) || (props[i] <= 0) )
        stop("Sorry, the proportion, props must lie within (0,1).")
}

if ((kappa0 >= 1) || (kappa0 <= 0) || (kappa1 <= 0) || (kappa1 >= 1))
        stop("Sorry, the null and alternative vallues of kappa must lie within (0,1).")

if ((alpha >= 1) || (alpha <= 0) || (power <= 0) || (power >= 1))
        stop("Sorry, the alpha and power must lie within (0,1).")

X <- NULL;
X$kappa0 <- kappa0;
X$kappa1 <- kappa1;
X$props <- props;
X$raters <- raters;
X$alpha <- alpha;
X$power <- power;

#2 raters
if (raters == 2)
{

P1 <- function(k, p)
{
x <- p[1]^2*(-(-p[1]-k+k*p[1])/p[1])
return(x)
}

P2 <- function(k, p)
{
x <- p[2]^2*(-(-p[2]-k+k*p[2])/p[2])
return(x)
}

P3 <- function(k, p)
{
x <- p[3]^2*(-(-p[3]-k+k*p[3])/p[3])
}

P0 <- function(k, p)
{
x <- 1 - P1(k,p) - P2(k,p) - P3(k,p);
return(x)
}


Results <- c(0,0,0,0);

Results[1] <- ( P0(k=kappa1, p = props) - P0(k=kappa0, p = props) )^2/ P0(k=kappa0, p = props)
Results[2] <- ( P1(k=kappa1, p = props) - P1(k=kappa0, p = props) )^2/ P1(k=kappa0, p = props)
Results[3] <- ( P2(k=kappa1, p = props) - P2(k=kappa0, p = props) )^2/ P2(k=kappa0, p = props)
Results[4] <- ( P3(k=kappa1, p = props) - P3(k=kappa0, p = props) )^2/ P3(k=kappa0, p = props)

}


#3 Raters;
if (raters == 3)
{

P1 <- function(k, p)
{
x <- p[1]^3*((p[1]^2+3*k*p[1]-2*k*p[1]^2+2*k^2-3*k^2*p[1]+k^2*p[1]^2)/p[1]^2)/(1+k)
return(x)
}

P2 <- function(k, p)
{
x <- p[2]^3*((p[2]^2+3*k*p[2]-2*k*p[2]^2+2*k^2-3*k^2*p[2]+k^2*p[2]^2)/p[2]^2)/(1+k)
return(x)
}

P3 <- function(k, p)
{
x <- p[3]^3*((p[3]^2+3*k*p[3]-2*k*p[3]^2+2*k^2-3*k^2*p[3]+k^2*p[3]^2)/p[3]^2)/(1+k)
return(x)
}

P0 <- function(k, p)
{
x <- 1 - P1(k,p) - P2(k,p) - P3(k,p);
return(x)
}

Results <- c(0,0,0,0);

Results[1] <- ( P0(k=kappa1, p = props) - P0(k=kappa0, p = props) )^2/ P0(k=kappa0, p = props)
Results[2] <- ( P1(k=kappa1, p = props) - P1(k=kappa0, p = props) )^2/ P1(k=kappa0, p = props)
Results[3] <- ( P2(k=kappa1, p = props) - P2(k=kappa0, p = props) )^2/ P2(k=kappa0, p = props)
Results[4] <- ( P3(k=kappa1, p = props) - P3(k=kappa0, p = props) )^2/ P3(k=kappa0, p = props)
}


#4 Raters;

if (raters == 4)
{

P1 <- function(k, p)
{
x <- -p[1]*(-p[1]^3-6*k*p[1]^2+3*k*p[1]^3-11*k^2*p[1]+12*k^2*p[1]^2-3*k^2*p[1]^3-6*k^3+11*k^3*p[1]-6*k^3*p[1]^2+k^3*p[1]^3)/((1+k)*(1+2*k))
return(x)
}

P2 <- function(k, p)
{
x <- -p[2]*(-p[2]^3-6*k*p[2]^2+3*k*p[2]^3-11*k^2*p[2]+12*k^2*p[2]^2-3*k^2*p[2]^3-6*k^3+11*k^3*p[2]-6*k^3*p[2]^2+k^3*p[2]^3)/((1+k)*(1+2*k))
return(x)
}

P3 <- function(k, p)
{
x <- -p[3]*(-p[3]^3-6*k*p[3]^2+3*k*p[3]^3-11*k^2*p[3]+12*k^2*p[3]^2-3*k^2*p[3]^3-6*k^3+11*k^3*p[3]-6*k^3*p[3]^2+k^3*p[3]^3)/((1+k)*(1+2*k))

return(x)
}

P0 <- function(k, p)
{
x <- 1 - P1(k,p) - P2(k,p) - P3(k,p);
return(x)
}

Results <- c(0,0,0,0);

Results[1] <- ( P0(k=kappa1, p = props) - P0(k=kappa0, p = props) )^2/ P0(k=kappa0, p = props)
Results[2] <- ( P1(k=kappa1, p = props) - P1(k=kappa0, p = props) )^2/ P1(k=kappa0, p = props)
Results[3] <- ( P2(k=kappa1, p = props) - P2(k=kappa0, p = props) )^2/ P2(k=kappa0, p = props)
Results[4] <- ( P3(k=kappa1, p = props) - P3(k=kappa0, p = props) )^2/ P3(k=kappa0, p = props)
}


#5 Raters;

if (raters == 5)
{

P1 <- function(k, p)
{
x <- p[1]*(p[1]^4+10*k*p[1]^3-4*k*p[1]^4+35*k^2*p[1]^2-30*k^2*p[1]^3+6*k^2*p[1]^4+50*k^3*p[1]-70*k^3*p[1]^2+30*k^3*p[1]^3-4*k^3*p[1]^4+24*k^4-50*k^4*p[1]+35*k^4*p[1]^2-10*k^4*p[1]^3+k^4*p[1]^4)/((1+k)*(1+2*k)*(1+3*k))
return(x)
}

P2 <- function(k, p)
{
x <- p[2]*(p[2]^4+10*k*p[2]^3-4*k*p[2]^4+35*k^2*p[2]^2-30*k^2*p[2]^3+6*k^2*p[2]^4+50*k^3*p[2]-70*k^3*p[2]^2+30*k^3*p[2]^3-4*k^3*p[2]^4+24*k^4-50*k^4*p[2]+35*k^4*p[2]^2-10*k^4*p[2]^3+k^4*p[2]^4)/((1+k)*(1+2*k)*(1+3*k))
return(x)
}

P3 <- function(k, p)
{
x <- p[3]*(p[3]^4+10*k*p[3]^3-4*k*p[3]^4+35*k^2*p[3]^2-30*k^2*p[3]^3+6*k^2*p[3]^4+50*k^3*p[3]-70*k^3*p[3]^2+30*k^3*p[3]^3-4*k^3*p[3]^4+24*k^4-50*k^4*p[3]+35*k^4*p[3]^2-10*k^4*p[3]^3+k^4*p[3]^4)/((1+k)*(1+2*k)*(1+3*k))
return(x)
}

P0 <- function(k, p)
{
x <- 1 - P1(k,p) - P2(k,p) - P3(k,p);
return(x)
}

Results <- c(0,0,0,0);

Results[1] <- ( P0(k=kappa1, p = props) - P0(k=kappa0, p = props) )^2/ P0(k=kappa0, p = props)
Results[2] <- ( P1(k=kappa1, p = props) - P1(k=kappa0, p = props) )^2/ P1(k=kappa0, p = props)
Results[3] <- ( P2(k=kappa1, p = props) - P2(k=kappa0, p = props) )^2/ P2(k=kappa0, p = props)
Results[4] <- ( P3(k=kappa1, p = props) - P3(k=kappa0, p = props) )^2/ P3(k=kappa0, p = props)
}

#6 Raters;
if (raters == 6)
{

P1 <- function(k, p)
{
x <- -p[1]*(-p[1]^5-15*k*p[1]^4+5*k*p[1]^5-85*k^2*p[1]^3-10*k^2*p[1]^5+60*k^2*p[1]^4-225*k^3*p[1]^2-90*k^3*p[1]^4+255*k^3*p[1]^3+10*k^3*p[1]^5-274*k^4*p[1]-255*k^4*p[1]^3-5*k^4*p[1]^5+450*k^4*p[1]^2+60*k^4*p[1]^4-120*k^5+274*k^5*p[1]-225*k^5*p[1]^2+85*k^5*p[1]^3-15*k^5*p[1]^4+k^5*p[1]^5)/((1+k)*(1+2*k)*(1+3*k)*(1+4*k))

return(x)
}

P2 <- function(k, p)
{
x <- -p[2]*(-p[2]^5-15*k*p[2]^4+5*k*p[2]^5-85*k^2*p[2]^3+60*k^2*p[2]^4-10*k^2*p[2]^5-225*k^3*p[2]^2+10*k^3*p[2]^5+255*k^3*p[2]^3-90*k^3*p[2]^4-274*k^4*p[2]-5*k^4*p[2]^5+60*k^4*p[2]^4+450*k^4*p[2]^2-255*k^4*p[2]^3-120*k^5+274*k^5*p[2]-225*k^5*p[2]^2+85*k^5*p[2]^3-15*k^5*p[2]^4+k^5*p[2]^5)/((1+k)*(1+2*k)*(1+3*k)*(1+4*k))
return(x)
}

P3 <- function(k, p)
{
x <- -p[3]*(-p[3]^5-15*k*p[3]^4+5*k*p[3]^5-85*k^2*p[3]^3-10*k^2*p[3]^5+60*k^2*p[3]^4-225*k^3*p[3]^2-90*k^3*p[3]^4+255*k^3*p[3]^3+10*k^3*p[3]^5-274*k^4*p[3]-255*k^4*p[3]^3-5*k^4*p[3]^5+450*k^4*p[3]^2+60*k^4*p[3]^4-120*k^5+274*k^5*p[3]-225*k^5*p[3]^2+85*k^5*p[3]^3-15*k^5*p[3]^4+k^5*p[3]^5)/((1+k)*(1+2*k)*(1+3*k)*(1+4*k))
return(x)
}

P0 <- function(k, p)
{
x <- 1 - P1(k,p) - P2(k,p) - P3(k,p);
return(x)
}

Results <- c(0,0,0,0);

Results[1] <- ( P0(k=kappa1, p = props) - P0(k=kappa0, p = props) )^2/ P0(k=kappa0, p = props)
Results[2] <- ( P1(k=kappa1, p = props) - P1(k=kappa0, p = props) )^2/ P1(k=kappa0, p = props)
Results[3] <- ( P2(k=kappa1, p = props) - P2(k=kappa0, p = props) )^2/ P2(k=kappa0, p = props)
Results[4] <- ( P3(k=kappa1, p = props) - P3(k=kappa0, p = props) )^2/ P3(k=kappa0, p = props)
}

X$N <- .powchi(1, alpha, power)[1]/sum(Results)

class(X) <- "Power3Cats";
return(X);

}


#ChiSquare Non-centrality functions, needed for sample size calculations...

.hichi <- function(chival,df,conf)
{
	uc <- c(chival,2*chival,3*chival)
	llim <- (1-conf)/2
	while(pchisq(chival,df,uc[1])<llim) {
		uc <- c(uc[1]/4,uc[1],uc[3])
	}
	while(pchisq(chival,df,uc[3])>llim) {
		uc <- c(uc[1],uc[3],uc[3]+chival)
	}

	diff <- 1
	while(diff > .00001) {
		if(pchisq(chival,df,uc[2])<llim) 
			uc <- c(uc[1],(uc[1]+uc[2])/2,uc[2]) 
			else uc <- c(uc[2],(uc[2]+uc[3])/2,uc[3])
		diff <- abs(pchisq(chival,df,uc[2]) - llim)
		lcdf <- pchisq(chival,df,uc[2])
	}
	c(uc[2],lcdf)
}

.powchi <- function(df,alpha,power)
{
	.hichi(qchisq(1-alpha,df),df,1-(1-power)*2)
}


#Print Method
print.Power3Cats <- function(x, ...)
{
cat("A minimum of", ceiling(x$N), "subjects are required for this study of interobserver agreement. \n")

for (i in 1:length(x$props))
{
if (x$props[i] * x$N < 5)
{
cat("Warning: At least one expected cell count is less than five. \n")
}
}

}

#Summary Method
summary.Power3Cats <- function(object, ...)
{
cat("Power-Based Sample Size Estimation for Studies of Interobserver Agreement with 3 Outcomes", "\n \n")
cat("Assuming:", "\n")
cat("Kappa0:", object$kappa0, "\n")
cat("Kappa1:", object$kappa1, "\n")
cat("Event Proportions:", object$props, "\n")
cat("Type I Error Rate (alpha) = ", object$alpha, " and Power = ", object$power, "\n \n",sep="")

cat("A minimum of", ceiling(object$N), "subjects are required for this study of interobserver agreement if ")
cat(object$raters, "raters are available. \n \n")

for (i in 1:length(object$props))
{
if (object$props[i] * object$N < 5)
{
cat("Warning: At least one expected cell count is less than five. \n")
}
}

}