mr.cv <-function(data, A, structural.cost, epsilon, sampling.cost, pilot=FALSE, m0=4, gamma=.49, verbose=FALSE)
{

# Internal functions (could be stand alone if one desired by simply 
# putting them in the workspace directly as their own functions).
########################################################################################################################
V2.n<-function(data)
{
if(is.data.frame(data)) 
{
data <- as.vector(data)
}
n <- length(data)
t1<-(sqrt(var(data)))^4/(4*mean(data)^4)
t2<-mu_4.n(data)/(4*mean(data)^4)
t3<-var(data)/(2*mean(data)^2)
t4<-mu_3.n(data)/(mean(data)^3)
return(t1+t2+t3-t4)
}


mu_4.n <- function(data)
{
if(is.data.frame(data)) 
{
data <- as.vector(data)
}
n <- length(data)
y2 <- data^2
y3 <- data^3
y4 <- data^4
w1 <- n^2*sum((data-mean(data))^4)
w2 <- (-2*n+3)*sum(y4)
w3 <- (8*n-12)*mean(data)*sum(y3)
w4 <- (-6+9/n)*(sum(y2))^2
return((w1+w2+w3+w4)/((n-1)*(n-2)*(n-3)))
}

mu_3.n <- function(data)
{
if(is.data.frame(data)) 
{
data <- as.vector(data)
}
n <- length(data)
w <- sum((data-mean(data))^3)
return(n*w/((n-1)*(n-2)))
}
########################################################################################################################

if(!missing(A))
{
if(!missing(structural.cost)) stop("Because you specified \'A\' directly, you should not also specify \'structural.cost\'.")
if(!missing(epsilon)) stop("Because you specified \'A\' directly, you should not also specify \'epsilon\'.")
if(A<=0) stop("A should be a non-zero positive value")
}

if(missing(A))
{
if(missing(structural.cost)) stop("Because you did not specificy \'A\' directly, you must specify \'structural.cost\'.")
if(missing(epsilon)) stop("Because you did not specificy \'A\' directly, you must specify \'epsilon\'.")
A <- structural.cost/(epsilon^2) 
if(A<=0) stop("A should be a non-zero positive value")
}

if(pilot==FALSE)
{
if(is.data.frame(data)) 
{
data <- as.vector(data)
}
n <- length(data)

Criterion <- (A/sampling.cost)*(V2.n(data)+n^(-2*gamma))

CV <- cv(mean=mean(data), sd = sqrt(var(data))) 

Rk<-(1/n)*A*V2.n(data)+(sampling.cost*n) # This is the risk function. 

if(verbose==FALSE) Outcome <- rbind(list(Risk=Rk, N=n,cv=CV, "Is.Satisfied?"=(n^2 >= Criterion)))
if(verbose==TRUE) Outcome <- rbind(list(Risk=Rk, N=n, cv=CV, V2.n=V2.n(data), mu_3.n=mu_3.n(data), mu_4.n=mu_4.n(data), Criterion=Criterion, Is.Satisfied=n^2 >= Criterion))
}

if(pilot==TRUE)
{
Outcome <- c(Pilot.SS=max(m0, ceiling((A/sampling.cost)^(1/((2+2*gamma))))))
}

return(Outcome)
}
