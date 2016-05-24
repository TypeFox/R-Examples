PowerBinary <- function(kappa0, kappa1, props, raters=2, alpha=0.05, power=0.80)
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
Results[1] <- ( P0(r=kappa1, p = props) - P0(r=kappa0, p = props) )^2/ P0(r=kappa0, p = props)
Results[2] <- ( P1(r=kappa1, p = props) - P1(r=kappa0, p = props) )^2/ P1(r=kappa0, p = props)
Results[3] <- ( P2(r=kappa1, p = props) - P2(r=kappa0, p = props) )^2/ P2(r=kappa0, p = props)
}


#3 Raters;
if (raters == 3)
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
Results[1] <- ( P0(r=kappa1, p = props) - P0(r=kappa0, p = props) )^2/ P0(r=kappa0, p = props)
Results[2] <- ( P1(r=kappa1, p = props) - P1(r=kappa0, p = props) )^2/ P1(r=kappa0, p = props)
Results[3] <- ( P2(r=kappa1, p = props) - P2(r=kappa0, p = props) )^2/ P2(r=kappa0, p = props)
Results[4] <- ( P3(r=kappa1, p = props) - P3(r=kappa0, p = props) )^2/ P3(r=kappa0, p = props)
}


#4 Raters;

if (raters == 4)
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
Results[1] <- ( P0(r=kappa1, p = props) - P0(r=kappa0, p = props) )^2/ P0(r=kappa0, p = props)
Results[2] <- ( P1(r=kappa1, p = props) - P1(r=kappa0, p = props) )^2/ P1(r=kappa0, p = props)
Results[3] <- ( P2(r=kappa1, p = props) - P2(r=kappa0, p = props) )^2/ P2(r=kappa0, p = props)
Results[4] <- ( P3(r=kappa1, p = props) - P3(r=kappa0, p = props) )^2/ P3(r=kappa0, p = props)
Results[5] <- ( P4(r=kappa1, p = props) - P4(r=kappa0, p = props) )^2/ P4(r=kappa0, p = props)
}


#5 Raters;

if (raters == 5)
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
Results[1] <- ( P0(r=kappa1, p = props) - P0(r=kappa0, p = props) )^2/ P0(r=kappa0, p = props)
Results[2] <- ( P1(r=kappa1, p = props) - P1(r=kappa0, p = props) )^2/ P1(r=kappa0, p = props)
Results[3] <- ( P2(r=kappa1, p = props) - P2(r=kappa0, p = props) )^2/ P2(r=kappa0, p = props)
Results[4] <- ( P3(r=kappa1, p = props) - P3(r=kappa0, p = props) )^2/ P3(r=kappa0, p = props)
Results[5] <- ( P4(r=kappa1, p = props) - P4(r=kappa0, p = props) )^2/ P4(r=kappa0, p = props)
Results[6] <- ( P5(r=kappa1, p = props) - P5(r=kappa0, p = props) )^2/ P5(r=kappa0, p = props)
}

#6 Raters;
if (raters == 6)
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

Results[1] <- ( P0(r=kappa1, p = props) - P0(r=kappa0, p = props) )^2/ P0(r=kappa0, p = props)
Results[2] <- ( P1(r=kappa1, p = props) - P1(r=kappa0, p = props) )^2/ P1(r=kappa0, p = props)
Results[3] <- ( P2(r=kappa1, p = props) - P2(r=kappa0, p = props) )^2/ P2(r=kappa0, p = props)
Results[4] <- ( P3(r=kappa1, p = props) - P3(r=kappa0, p = props) )^2/ P3(r=kappa0, p = props)
Results[5] <- ( P4(r=kappa1, p = props) - P4(r=kappa0, p = props) )^2/ P4(r=kappa0, p = props)
Results[6] <- ( P5(r=kappa1, p = props) - P5(r=kappa0, p = props) )^2/ P5(r=kappa0, p = props)
Results[7] <- ( P6(r=kappa1, p = props) - P6(r=kappa0, p = props) )^2/ P6(r=kappa0, p = props)
}

X$N <- .powchi(1, alpha, power)[1]/sum(Results)

class(X) <- "PowerBinary";
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
print.PowerBinary <- function(x, ...)
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
summary.PowerBinary <- function(object, ...)
{
cat("Power-Based Sample Size Estimation for Studies of Interobserver Agreement with a Binary Outcome", "\n \n")
cat("Assuming:", "\n")
cat("Kappa0:", object$kappa0, "\n")
cat("Kappa1:", object$kappa1, "\n")
cat("Event Proportion:", object$props, "\n")
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