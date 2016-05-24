CIBinary <- function(kappa0, kappaL, kappaU=NA, props, raters=2, alpha=0.05)
{
#Error Checking
if ( (raters != 1) && (raters !=2) && (raters !=3) && (raters != 4) && (raters != 5) && (raters != 6) )
        stop("Sorry, this function is designed for between 2 to 6 raters.")

if (length(props) == 1)
{
if ((props >= 1) || (props <= 0) )
        stop("Sorry, the proportion, props must lie within (0,1).")
}

if (length(props) == 2)
{
if ( abs( sum(props) - 1) >= 0.001 )
	stop("Sorry, the two proportions must sum to one.")

for (i in 1:2)
{
if ((props[i] >= 1) || (props[i] <= 0) )
        stop("Sorry, the proportion, props must lie within (0,1).")
}
props <- props[1];
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

P0 <- function(r, p)
{
x <- (1- p)^2 + r*p*(1 - p) 
return(x)
}

P1 <- function(r, p)
{
x <- 2*(1 - r)*p*(1 - p) 
return(x)
}

P2 <- function(r, p)
{
x <- p^2 + r*p*(1 - p) 
return(x)
}

Results <- c(0,0,0)
Results[1] <- ( (n*P0(r=rho0, p = Pi)) - (n*P0(r=rho1, p = Pi)) )^2/ (n*P0(r=rho1, p = Pi))
Results[2] <- ( (n*P1(r=rho0, p = Pi)) - (n*P1(r=rho1, p = Pi)) )^2/ (n*P1(r=rho1, p = Pi))
Results[3] <- ( (n*P2(r=rho0, p = Pi)) - (n*P2(r=rho1, p = Pi)) )^2/ (n*P2(r=rho1, p = Pi))

return(sum(Results,na.rm=TRUE))

}

}


#3 Raters;
if (raters == 3)
{

.CalcIT <- function(rho0, rho1, Pi, n)
{

P0 <- function(r, p)
{
x <- (1- p)^3 + p*r*( (1 - p)^2 + (1 - p) )
return(x)
}

P1 <- function(r, p)
{
x <- 3*p*(1 - r)*(1 - p)^2
return(x)
}

P2 <- function(r, p)
{
x <- 3*p^2*(1 - p)*(1 - r)
return(x)
}

P3 <- function(r, p)
{
x <- p^3 + r*p*(1 - p^2)
return(x)
}

Results <- c(0,0,0,0);
Results[1] <- ( n*P0(r=rho0, p = Pi) - n*P0(r=rho1, p = Pi) )^2/ (n*P0(r=rho1, p = Pi))
Results[2] <- ( n*P1(r=rho0, p = Pi) - n*P1(r=rho1, p = Pi) )^2/ (n*P1(r=rho1, p = Pi))
Results[3] <- ( n*P2(r=rho0, p = Pi) - n*P2(r=rho1, p = Pi) )^2/ (n*P2(r=rho1, p = Pi))
Results[4] <- ( n*P3(r=rho0, p = Pi) - n*P3(r=rho1, p = Pi) )^2/ (n*P3(r=rho1, p = Pi))

return(sum(Results,na.rm=TRUE))

}

}
#4 Raters;

if (raters == 4)
{


.CalcIT <- function(rho0, rho1, Pi, n)
{

P0 <- function(r, p)
{
x <- p^4 -r*p^4 -4*p^3 +4*r*p^3 +6*p^2 -6*r*p^2 -4*p +3*r*p +1
return(x)
}

P1 <- function(r, p)
{
x <- (4*(1-3*p-r+3*r*p+3*p^2-3*r*p^2-p^3+r*p^3))*p
return(x)
}

P2 <- function(r, p)
{
x <- -(6*(-1+r+2*p-2*r*p-p^2+r*p^2))*p^2
return(x)
}

P3 <- function(r, p)
{
x <- (4*(1-r-p+r*p))*p^3
return(x)
}

P4 <- function(r,p)
{
x <- -(-p^3-r+r*p^3)*p
return(x)
}

Results <- c(0,0,0,0,0);
Results[1] <- ( n*P0(r=rho0, p = Pi) - n*P0(r=rho1, p = Pi) )^2/ (n*P0(r=rho1, p = Pi))
Results[2] <- ( n*P1(r=rho0, p = Pi) - n*P1(r=rho1, p = Pi) )^2/ (n*P1(r=rho1, p = Pi))
Results[3] <- ( n*P2(r=rho0, p = Pi) - n*P2(r=rho1, p = Pi) )^2/ (n*P2(r=rho1, p = Pi))
Results[4] <- ( n*P3(r=rho0, p = Pi) - n*P3(r=rho1, p = Pi) )^2/ (n*P3(r=rho1, p = Pi))
Results[5] <- ( n*P4(r=rho0, p = Pi) - n*P4(r=rho1, p = Pi) )^2/ (n*P4(r=rho1, p = Pi))

return(sum(Results,na.rm=TRUE))

}
}


#5 Raters;

if (raters == 5)
{


.CalcIT <- function(rho0, rho1, Pi, n)
{

P0 <- function(r, p)
{
x <- -p^5+r*p^5+5*p^4-5*r*p^4+10*r*p^3-10*p^3-10*r*p^2+10*p^2+4*r*p-5*p+1
return(x)
}
P1 <- function(r, p)
{
x <- -(5*(-1+4*p+r-4*r*p-6*p^2+6*r*p^2+4*p^3-4*r*p^3-p^4+r*p^4))*p
return(x)
}

P2 <- function(r, p)
{
x <- (10*(1-r-3*p+3*r*p+3*p^2-3*r*p^2-p^3+r*p^3))*p^2
return(x)
}

P3 <- function(r, p)
{
x <- -(10*(-1+r+2*p-2*r*p-p^2+r*p^2))*p^3
return(x)
}

P4 <- function(r,p)
{
x <- (5*(1-r-p+r*p))*p^4
return(x)
}

P5 <- function(r,p)
{
x <- -(-p^4-r+r*p^4)*p
return(x)
}

Results <- c(0,0,0,0,0,0);
Results[1] <- ( n*P0(r=rho0, p = Pi) - n*P0(r=rho1, p = Pi) )^2/ (n*P0(r=rho1, p = Pi))
Results[2] <- ( n*P1(r=rho0, p = Pi) - n*P1(r=rho1, p = Pi) )^2/ (n*P1(r=rho1, p = Pi))
Results[3] <- ( n*P2(r=rho0, p = Pi) - n*P2(r=rho1, p = Pi) )^2/ (n*P2(r=rho1, p = Pi))
Results[4] <- ( n*P3(r=rho0, p = Pi) - n*P3(r=rho1, p = Pi) )^2/ (n*P3(r=rho1, p = Pi))
Results[5] <- ( n*P4(r=rho0, p = Pi) - n*P4(r=rho1, p = Pi) )^2/ (n*P4(r=rho1, p = Pi))
Results[6] <- ( n*P5(r=rho0, p = Pi) - n*P5(r=rho1, p = Pi) )^2/ (n*P5(r=rho1, p = Pi))

return(sum(Results,na.rm=TRUE))
}
}
#6 Raters;
if (raters == 6)
{

.CalcIT <- function(rho0, rho1, Pi, n)
{

P0 <- function(r, p)
{
x <- p^6-r*p^6+6*r*p^5-6*p^5-15*r*p^4+15*p^4+20*r*p^3-20*p^3-15*r*p^2+15*p^2+5*r*p-6*p+1
return(x)
}

P1 <- function(r, p)
{
x <- (6*(1-5*p-r+5*r*p+10*p^2-10*r*p^2-10*p^3+10*r*p^3+5*p^4-5*r*p^4-p^5+r*p^5))*p
return(x)
}

P2 <- function(r, p)
{
x <- -(15*(-1+r+4*p-4*r*p-6*p^2+6*r*p^2+4*p^3-4*r*p^3-p^4+r*p^4))*p^2
return(x)
}

P3 <- function(r, p)
{
x <- (20*(1-r-3*p+3*r*p+3*p^2-3*r*p^2-p^3+r*p^3))*p^3
return(x)
}
0
P4 <- function(r,p)
{
x <- -(15*(-1+r+2*p-2*r*p-p^2+r*p^2))*p^4
return(x)
}

P5 <- function(r,p)
{
x <- (6*(1-r-p+r*p))*p^5
return(x)
}

P6 <- function(r,p)
{
x <- -(-p^5-r+r*p^5)*p
return(x)
}

Results <- c(0,0,0,0,0,0,0);

Results[1] <- ( n*P0(r=rho0, p = Pi) - n*P0(r=rho1, p = Pi) )^2/ (n*P0(r=rho1, p = Pi))
Results[2] <- ( n*P1(r=rho0, p = Pi) - n*P1(r=rho1, p = Pi) )^2/ (n*P1(r=rho1, p = Pi))
Results[3] <- ( n*P2(r=rho0, p = Pi) - n*P2(r=rho1, p = Pi) )^2/ (n*P2(r=rho1, p = Pi))
Results[4] <- ( n*P3(r=rho0, p = Pi) - n*P3(r=rho1, p = Pi) )^2/ (n*P3(r=rho1, p = Pi))
Results[5] <- ( n*P4(r=rho0, p = Pi) - n*P4(r=rho1, p = Pi) )^2/ (n*P4(r=rho1, p = Pi))
Results[6] <- ( n*P5(r=rho0, p = Pi) - n*P5(r=rho1, p = Pi) )^2/ (n*P5(r=rho1, p = Pi))
Results[7] <- ( n*P6(r=rho0, p = Pi) - n*P6(r=rho1, p = Pi) )^2/ (n*P6(r=rho1, p = Pi))

return(sum(Results,na.rm=TRUE))
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

class(X) <- "CIBinary";
return(X);

}


#Print Method
print.CIBinary <- function(x, ...)
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
summary.CIBinary <- function(object, ...)
{
cat("Confidence-Interval Based Sample Size Estimation for \n")
cat("Studies of Interobserver Agreement with a Binary Outcome \n \n")
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