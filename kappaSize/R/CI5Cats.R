CI5Cats <- function(kappa0, kappaL, kappaU=NA, props, raters=2, alpha=0.05)
{

#Error Checking
if ( (raters != 1) && (raters !=2) && (raters !=3) && (raters != 4) && (raters != 5) && (raters != 6) )
        stop("Sorry, this function is designed for between 2 to 6 raters.")

if (length(props) != 5)
	stop("Sorry, you must enter the anticipated proportions of each of the five categories.")

if ( abs( sum(props) - 1) >= 0.001 )
	stop("Sorry, the three proportions must sum to one.")

for (i in 1:5)
{
if ((props[i] >= 1) || (props[i] <= 0) )
        stop("Sorry, the proportion, props must lie within (0,1).")
}

if ((kappa0 >= 1) || (kappa0 <= 0) || (kappaL <= 0) || (kappaL >= 1) )
        stop("Sorry, the null and lower values of kappa must lie within (0,1).")

if ((!is.na(kappaU)) && ( (kappaU <= 0) || (kappaU >= 1) ) )
	stop("Sorry, the upper value of kappa must be either NA for a one-sided test or within (0,1).");

if ( (kappaL >= kappa0) || ( (!is.na(kappaU)) && (kappaU <= kappa0) ) )
	stop("Remember, kappaL < kappa0 < kappaU...")

if ( (is.na(kappaL)) && (is.na(kappaU)) )
	stop("Sorry, at least one of kappaL or kappaU must be specified.")

if ( (alpha >= 1) || (alpha <= 0) )
        stop("Sorry, the alpha and power must lie within (0,1).")

X <- NULL;
X$kappa0 <- kappa0;
X$kappaL <- kappaL;
X$kappaU <- kappaU;
X$props <- props;
X$raters <- raters;
X$alpha <- alpha;

if (is.na(kappaU))
{
X$ChiCrit <- qchisq((1-2*alpha),1);
}

if ( (!is.na(kappaL)) && (!is.na(kappaU)) )
{
X$ChiCrit <- qchisq((1-alpha),1);
}


#2 raters
if (raters == 2)
{

.CalcIT <- function(rho0, rho1, Pi, n)
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

P4 <- function(k, p)
{
x <- p[4]^2*(-(-p[4]-k+k*p[4])/p[4])
}

P5 <- function(k, p)
{
x <- p[5]^2*(-(-p[5]-k+k*p[5])/p[5])
}

P0 <- function(k, p)
{
x <- 1 - P1(k,p) - P2(k,p) - P3(k,p) -P4(k,p) - P5(k,p);
return(x)
}

results <- c(0,0,0,0,0,0);

results[1] <- ( (n*P0(k=rho0, p=Pi)) - (n*P0(k=rho1, p = Pi)) )^2/ (n*P0(k=rho1, p = Pi))
results[2] <- ( (n*P1(k=rho0, p=Pi)) - (n*P1(k=rho1, p = Pi)) )^2/ (n*P1(k=rho1, p = Pi))
results[3] <- ( (n*P2(k=rho0, p=Pi)) - (n*P2(k=rho1, p = Pi)) )^2/ (n*P2(k=rho1, p = Pi))
results[4] <- ( (n*P3(k=rho0, p=Pi)) - (n*P3(k=rho1, p = Pi)) )^2/ (n*P3(k=rho1, p = Pi))
results[5] <- ( (n*P4(k=rho0, p=Pi)) - (n*P4(k=rho1, p = Pi)) )^2/ (n*P4(k=rho1, p = Pi))
results[6] <- ( (n*P5(k=rho0, p=Pi)) - (n*P5(k=rho1, p = Pi)) )^2/ (n*P5(k=rho1, p = Pi))

return(sum(results,na.rm=TRUE))

}

}


#3 Raters;
if (raters == 3)
{

.CalcIT <- function(rho0, rho1, Pi, n)
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


P4 <- function(k, p)
{
x <- p[4]^3*((p[4]^2+3*k*p[4]-2*k*p[4]^2+2*k^2-3*k^2*p[4]+k^2*p[4]^2)/p[4]^2)/(1+k)
return(x)
}

P5 <- function(k, p)
{
x <- p[5]^3*((p[5]^2+3*k*p[5]-2*k*p[5]^2+2*k^2-3*k^2*p[5]+k^2*p[5]^2)/p[5]^2)/(1+k)
return(x)
}

P0 <- function(k, p)
{
x <- 1 - P1(k,p) - P2(k,p) - P3(k,p) -P4(k,p) - P5(k,p);
return(x)
}

results <- c(0,0,0,0,0,0);

results[1] <- ( (n*P0(k=rho0, p=Pi)) - (n*P0(k=rho1, p = Pi)) )^2/ (n*P0(k=rho1, p = Pi))
results[2] <- ( (n*P1(k=rho0, p=Pi)) - (n*P1(k=rho1, p = Pi)) )^2/ (n*P1(k=rho1, p = Pi))
results[3] <- ( (n*P2(k=rho0, p=Pi)) - (n*P2(k=rho1, p = Pi)) )^2/ (n*P2(k=rho1, p = Pi))
results[4] <- ( (n*P3(k=rho0, p=Pi)) - (n*P3(k=rho1, p = Pi)) )^2/ (n*P3(k=rho1, p = Pi))
results[5] <- ( (n*P4(k=rho0, p=Pi)) - (n*P4(k=rho1, p = Pi)) )^2/ (n*P4(k=rho1, p = Pi))
results[6] <- ( (n*P5(k=rho0, p=Pi)) - (n*P5(k=rho1, p = Pi)) )^2/ (n*P5(k=rho1, p = Pi))

return(sum(results,na.rm=TRUE))

}

}
#4 Raters;

if (raters == 4)
{

.CalcIT <- function(rho0, rho1, Pi, n)
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

P4 <- function(k, p)
{
x <- -p[4]*(-p[4]^3-6*k*p[4]^2+3*k*p[4]^3-11*k^2*p[4]+12*k^2*p[4]^2-3*k^2*p[4]^3-6*k^3+11*k^3*p[4]-6*k^3*p[4]^2+k^3*p[4]^3)/((1+k)*(1+2*k))
return(x)
}

P5 <- function(k, p)
{
x <- -p[5]*(-p[5]^3-6*k*p[5]^2+3*k*p[5]^3-11*k^2*p[5]+12*k^2*p[5]^2-3*k^2*p[5]^3-6*k^3+11*k^3*p[5]-6*k^3*p[5]^2+k^3*p[5]^3)/((1+k)*(1+2*k))
return(x)
}

P0 <- function(k, p)
{
x <- 1 - P1(k,p) - P2(k,p) - P3(k,p) - P4(k,p) - P5(k,p);
return(x)
}

results <- c(0,0,0,0,0,0);

results[1] <- ( (n*P0(k=rho0, p=Pi)) - (n*P0(k=rho1, p = Pi)) )^2/ (n*P0(k=rho1, p = Pi))
results[2] <- ( (n*P1(k=rho0, p=Pi)) - (n*P1(k=rho1, p = Pi)) )^2/ (n*P1(k=rho1, p = Pi))
results[3] <- ( (n*P2(k=rho0, p=Pi)) - (n*P2(k=rho1, p = Pi)) )^2/ (n*P2(k=rho1, p = Pi))
results[4] <- ( (n*P3(k=rho0, p=Pi)) - (n*P3(k=rho1, p = Pi)) )^2/ (n*P3(k=rho1, p = Pi))
results[5] <- ( (n*P4(k=rho0, p=Pi)) - (n*P4(k=rho1, p = Pi)) )^2/ (n*P4(k=rho1, p = Pi))
results[6] <- ( (n*P5(k=rho0, p=Pi)) - (n*P5(k=rho1, p = Pi)) )^2/ (n*P5(k=rho1, p = Pi))

return(sum(results,na.rm=TRUE))

}
}


#5 Raters;

if (raters == 5)
{

.CalcIT <- function(rho0, rho1, Pi, n)
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

P4 <- function(k, p)
{
x <- p[4]*(p[4]^4+10*k*p[4]^3-4*k*p[4]^4+35*k^2*p[4]^2-30*k^2*p[4]^3+6*k^2*p[4]^4+50*k^3*p[4]-70*k^3*p[4]^2+30*k^3*p[4]^3-4*k^3*p[4]^4+24*k^4-50*k^4*p[4]+35*k^4*p[4]^2-10*k^4*p[4]^3+k^4*p[4]^4)/((1+k)*(1+2*k)*(1+3*k))
return(x)
}

P5 <- function(k, p)
{
x <- p[5]*(p[5]^4+10*k*p[5]^3-4*k*p[5]^4+35*k^2*p[5]^2-30*k^2*p[5]^3+6*k^2*p[5]^4+50*k^3*p[5]-70*k^3*p[5]^2+30*k^3*p[5]^3-4*k^3*p[5]^4+24*k^4-50*k^4*p[5]+35*k^4*p[5]^2-10*k^4*p[5]^3+k^4*p[5]^4)/((1+k)*(1+2*k)*(1+3*k))
return(x)
}

P0 <- function(k, p)
{
x <- 1 - P1(k,p) - P2(k,p) - P3(k,p) - P4(k,p) - P5(k,p);
return(x)
}

results <- c(0,0,0,0,0,0);

results[1] <- ( (n*P0(k=rho0, p=Pi)) - (n*P0(k=rho1, p = Pi)) )^2/ (n*P0(k=rho1, p = Pi))
results[2] <- ( (n*P1(k=rho0, p=Pi)) - (n*P1(k=rho1, p = Pi)) )^2/ (n*P1(k=rho1, p = Pi))
results[3] <- ( (n*P2(k=rho0, p=Pi)) - (n*P2(k=rho1, p = Pi)) )^2/ (n*P2(k=rho1, p = Pi))
results[4] <- ( (n*P3(k=rho0, p=Pi)) - (n*P3(k=rho1, p = Pi)) )^2/ (n*P3(k=rho1, p = Pi))
results[5] <- ( (n*P4(k=rho0, p=Pi)) - (n*P4(k=rho1, p = Pi)) )^2/ (n*P4(k=rho1, p = Pi))
results[6] <- ( (n*P5(k=rho0, p=Pi)) - (n*P5(k=rho1, p = Pi)) )^2/ (n*P5(k=rho1, p = Pi))

return(sum(results,na.rm=TRUE))
}
}


#6 Raters;
if (raters == 6)
{

.CalcIT <- function(rho0, rho1, Pi, n)
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

P4 <- function(k, p)
{
x <- -p[4]*(-p[4]^5-15*k*p[4]^4+5*k*p[4]^5-85*k^2*p[4]^3-10*k^2*p[4]^5+60*k^2*p[4]^4-225*k^3*p[4]^2-90*k^3*p[4]^4+255*k^3*p[4]^3+10*k^3*p[4]^5-274*k^4*p[4]-255*k^4*p[4]^3-5*k^4*p[4]^5+450*k^4*p[4]^2+60*k^4*p[4]^4-120*k^5+274*k^5*p[4]-225*k^5*p[4]^2+85*k^5*p[4]^3-15*k^5*p[4]^4+k^5*p[4]^5)/((1+k)*(1+2*k)*(1+3*k)*(1+4*k))
return(x)
}

P5 <- function(k, p)
{
x <- -p[5]*(-p[5]^5-15*k*p[5]^4+5*k*p[5]^5-85*k^2*p[5]^3-10*k^2*p[5]^5+60*k^2*p[5]^4-225*k^3*p[5]^2-90*k^3*p[5]^4+255*k^3*p[5]^3+10*k^3*p[5]^5-274*k^4*p[5]-255*k^4*p[5]^3-5*k^4*p[5]^5+450*k^4*p[5]^2+60*k^4*p[5]^4-120*k^5+274*k^5*p[5]-225*k^5*p[5]^2+85*k^5*p[5]^3-15*k^5*p[5]^4+k^5*p[5]^5)/((1+k)*(1+2*k)*(1+3*k)*(1+4*k))
return(x)
}

P0 <- function(k, p)
{
x <- 1 - P1(k,p) - P2(k,p) - P3(k,p) - P4(k,p) - P5(k,p);
return(x)
}

results <- c(0,0,0,0,0,0);

results[1] <- ( (n*P0(k=rho0, p=Pi)) - (n*P0(k=rho1, p = Pi)) )^2/ (n*P0(k=rho1, p = Pi))
results[2] <- ( (n*P1(k=rho0, p=Pi)) - (n*P1(k=rho1, p = Pi)) )^2/ (n*P1(k=rho1, p = Pi))
results[3] <- ( (n*P2(k=rho0, p=Pi)) - (n*P2(k=rho1, p = Pi)) )^2/ (n*P2(k=rho1, p = Pi))
results[4] <- ( (n*P3(k=rho0, p=Pi)) - (n*P3(k=rho1, p = Pi)) )^2/ (n*P3(k=rho1, p = Pi))
results[5] <- ( (n*P4(k=rho0, p=Pi)) - (n*P4(k=rho1, p = Pi)) )^2/ (n*P4(k=rho1, p = Pi))
results[6] <- ( (n*P5(k=rho0, p=Pi)) - (n*P5(k=rho1, p = Pi)) )^2/ (n*P5(k=rho1, p = Pi))

return(sum(results,na.rm=TRUE))
}
}

n <- 10;

if ( (!is.na(kappaL)) && (!is.na(kappaU)) )
{
resultsl <- 0;
resultsu <- 0;

while ( (abs(resultsl - 0.001) < X$ChiCrit) || (abs(resultsu - 0.001) < X$ChiCrit) )
{
n <- n + 1;
resultsl <- .CalcIT(rho0=kappa0, rho1 = kappaL, Pi=props, n=n)
resultsu <- .CalcIT(rho0=kappa0, rho1 = kappaU, Pi=props, n=n)

if (is.infinite(resultsu))
{
resultsu <- 0
}

if (is.infinite(resultsl))
{
resultsl <- 0
}

}

}


if (is.na(kappaU))
{
resultsl <- 0;

while (abs(resultsl - 0.001) < X$ChiCrit) 
{
n <- n + 1;
resultsl <- .CalcIT(rho0=kappa0, rho1 = kappaL, Pi=props, n=n)

if (is.infinite(resultsl))
{
resultsl <- 0
}

}

}

X$n <- n;

class(X) <- "CI5Cats";
return(X);

}


#Print Method
print.CI5Cats <- function(x, ...)
{
cat("A minimum of", max(x$n, x$n), "subjects are required for this study of interobserver agreement. \n")


for (i in 1:length(x$props))
{
if (x$props[i] * x$n < 5)
{
cat("Warning: At least one expected cell count is less than five. \n")
}
}

}

#Summary Method
summary.CI5Cats <- function(object, ...)
{
cat("Confidence-Interval Based Sample Size Estimation for \n")
cat("Studies of Interobserver Agreement with 5 Categorical Outcomes \n \n")
cat("Assuming:", "\n")
cat("Kappa0:", object$kappa0, "\n")
cat("KappaL:", object$kappaL, "\n")
cat("KappaU:", object$kappaU, "\n")
cat("Event Proportion:", object$props, "\n")
cat("Type I Error Rate (alpha) = ", object$alpha, "\n \n")

if (is.na(object$kappaU))
{
cat("A minimum of", ceiling(object$n), "subjects are required to ensure the lower \n")
cat("confidence limit is at least ", object$kappaL, ". \n \n", sep="")


for (i in 1:length(object$props))
{
if (object$props[i] * object$n < 5)
{
cat("Warning: At least one expected cell count is less than five. \n")
}
}

}

if ( (!is.na(object$kappaL)) && (!is.na(object$kappaU)) )
{
cat("A minimum of", ceiling(object$n), "subjects are required to ensure the lower \n")
cat("confidence limit is at least", object$kappaL, "and the upper confidence \n")
cat("limit does not exceed ", object$kappaU, ". \n \n", sep="")


for (i in 1:length(object$props))
{
if (object$props[i] * object$n < 5)
{
cat("Warning: At least one expected cell count is less than five. \n")
}
}

}

}