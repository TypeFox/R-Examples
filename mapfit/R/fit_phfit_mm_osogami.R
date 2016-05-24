
mixtureofPH <- function(tau1, T1, xi1, tau2, T2, xi2, p) {

tau <- c(p*tau1, (1-p)*tau2)

dim1 <- dim(T1)[1]
dim2 <- dim(T2)[2]
T <- rbind(cbind(as.matrix(T1), matrix(0, dim1, dim2)),
		cbind(matrix(0, dim2, dim1), as.matrix(T2)))
xi <- c(xi1, xi2)
list(tau=tau, T=T, xi=xi)
}

convolutionofPH <- function(tau1, T1, xi1, tau2, T2, xi2) {
	tau0 <- 1 - sum(tau1)
	if (tau0 == 0) {
		tau <- c(tau1, rep(0, length(tau2)))
	} else {
		tau <- c(tau1, tau0*tau2)
	}

	dim1 <- dim(T1)[1]
	t <- matrix(xi1, dim1, 1)
  dim2 <- dim(T2)[1]    
	T <- rbind(
		cbind(T1, t %*% matrix(tau2,1,dim2)),
		cbind(matrix(0, dim2, dim1), T2))
  xi <- c(rep(0,length(xi1)), xi2)
  list(tau=tau, T=T, xi=xi)
}

isPH2 <- function(m2, m3) {

##%% If the input is not even a PH, return 0;
if( (m2 <= 1) || (m3 <= m2) ) {
	return(FALSE)
}

##%% Main:
if( m2 > 2 )  {
 if( m3 > 3/2 * m2 ) {
     return(TRUE)
 }
} else if( m2 < 2 ) {
 if( (m2 >= 3/2) && (m3 <= 6*(m2-1)/m2) && (m3 >= (9*m2-12+3*(2-m2)*sqrt(2*(2-m2)))/m2) ) {
     return(TRUE)
 }
}
return(FALSE)
}

matching3PH2 <- function(mean, m2, m3, epsilon=0.001) {

## If the input moment is zero, return error.

if( mean == 0 ) {
	stop(sprintf("Illegal inputs: %f %f %f", mean, m2, m3))
}

##%% If there is no 2-phase PH distribution with the input moments,
##%% calculate the normalized moments of the "closest" 2-phase PH distribution.

if( m2 < 3/2 ) {

  m2 <- 3/2
  m3 <- 2

} else if( m2 >= 2 - epsilon^2 ) {

  if( abs(m2-2) < epsilon^2 ) {
    m2 <- (1+epsilon) * 2
  }

  if( m3 <= 3/2 * m2 + epsilon^2 ) {
    m3 <- (1+epsilon) * 3/2 * m2
  }

} else {

  if( m3 > 6*(m2-1)/m2 ) {
   m3 <- 6*(m2-1)/m2 ## m2?
   m3 <- 2*m2 - 1
  } else if( m3 < (9*m2-12+3*(2-m2)*sqrt(2*(2-m2))) / m2 ) {
   m3 <- (9*m2-12+3*(2-m2)*sqrt(2*(2-m2))) / m2 ## m2?
   m3 <- 2*m2 - 1
  }

}

u <- (6-2*m3) / (3*m2-2*m3)
v <- (12-6*m2) / (m2*(3*m2-2*m3))

sqrtD <- sqrt(u^2-4*v)

mu1 <- (u+sqrtD)/2
mu2 <- (u-sqrtD)/2
p <- mu2*(mu1-1)/mu1

MU1 <- mu1 / mean
MU2 <- mu2 / mean

tau <- c(1, 0)
T <- matrix(c(-MU1, 0, p*MU1, -MU2), 2, 2)
xi <- c((1-p)*MU1,MU2)

list(tau=tau, T=T, xi=xi)
}

matching3EC <- function(mean, m2, m3, epsilon=0.001) {

if( isPH2(m2,m3) ) {

  return(matching3PH2(mean, m2, m3, epsilon=epsilon))

}

## determined the mass probability at zero.

if( (m3 > 2*m2-1) && (1/(m2-1)-floor(1/(m2-1)) < epsilon^2) ) {

  p <- (m2*m2+2*m2-1) / (2*m2*m2)

} else if( m3 < 2*m2-1 ) {

  p <- 1 / (2*m2-m3)

} else {

  p <- 1

}

mean <- mean / p
m2 <- m2 * p
m3 <- m3 * p


## determine the number of phases, n

 if( (m3 == 2*m2 - 1) && (m2 <= 2) ) {

  n <- floor(m2/(m2-1) + epsilon^2)

 } else {

  n <- floor(m2/(m2-1) + 1 - epsilon^2)

 }


## calculate the mean and the normalized moments, (meanX,m2X,m3X), of 2-phase PH Coxian, X.

 m2X <- ((n-3)*m2-(n-2)) / ((n-2)*m2-(n-1))
 meanX <- mean / ((n-2)*m2X-(n-3))
 alpha <- (n-2)*(m2X-1)*(n*(n-1)*m2X^2-n*(2*n-5)*m2X+(n-1)*(n-3))
 beta <- ((n-1)*m2X-(n-2)) * ((n-2)*m2X-(n-3))^2
 m3X <- (beta*m3-alpha) / m2X

## calculate the rate of Erlang

 lambda <- 1 / ((m2X-1)*meanX)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%% When X is exponential, EC is just an Erlang %%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if( (m2X == 2) && (m3X == 3) ) {

  mu <- 1 / meanX
  T <- matrix(0, n-1, n-1)
  if (n > 2) {
  for (i in 1:(n-2)) {
   T[i,i] <- -mu
   T[i,i+1] <- mu
  }
  }
  T[n-1,n-1] <- -mu

  tau <- numeric(n-1)
  tau[1] <- p

  xi <- numeric(n-1)
  xi[n-1] <- mu

  list(tau=tau, T=T, xi=xi)

}

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%% When X is not exponential (2-phase PH) %%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phX <- matching3PH2(meanX, m2X, m3X, epsilon=epsilon)

T <- rbind(cbind(matrix(0, n-2,n-2), matrix(0, n-2,2)),
	cbind(matrix(0, 2, n-2), phX$T))
if (n > 2) {
for (i in 1:(n-2)) {
  T[i,i] <- -lambda
  T[i,i+1] <- lambda
}
}

tau <- numeric(n)
tau[1] <- p

xi <- c(rep(0,n-2), phX$xi)

list(tau=tau, T=T, xi=xi)
}

matching3PH <- function(momentsofG, max.phase = 50, epsilon = 0.001) {

## If the input moment is zero, return error.

if( (momentsofG[1] == 0) || (momentsofG[2] == 0) ) {
	stop(sprintf("Illegal inputs: %f %f %f", momentsofG[1], momentsofG[2], momentsofG[3]))
}

## Calculate the normalized moments of the input, G.

m2 <- momentsofG[2] / momentsofG[1]^2;
m3 <- momentsofG[3] / (momentsofG[1] * momentsofG[2]);

## Any PH distribution has m2 > 1 and m3 > m2.
## If there is no PH for the input, return error.

if( (m2 <= 1 + epsilon^2) || (m3 <= m2 + epsilon^2) ) {
	stop(sprintf("Illegal inputs: %f %f %f", momentsofG[1], momentsofG[2], momentsofG[3]))
}

## If the necessary number of phases exceeds the limit, MAX_PHASE,
## calculate the normalized moments of the "closest" PH distribution.

smallestM2 <- (max.phase/(max.phase-1) + (max.phase-1)/(max.phase-2)) / 2;
if( m2 < smallestM2 ) {
  m2 <- smallestM2
}

smallestM3 <- (max.phase+2)/(max.phase+1) * m2;
if( m3 < smallestM3 ) {
  m3 <- smallestM3
}

## If there is no PH distribution with the input moments,
## calculate the normalized moments of the "closest" PH distribution.

if( (m3 < 2*m2 - 1) && (abs(m3 - 3/2 * m2) < epsilon^2) ) {

  m3 <- (1+epsilon) * m3

} else if ( m3 > 2*m2 - 1 ) {

 r <- m3 / (2*m2-1)
 n <- floor(1/(m2-1))

 if( abs(1/(m2-1) - n) < epsilon^2 ) {
   m2 <- 1 + 1 / (n*(1-epsilon))
   m3 <- r * (2*m2-1)
 }

}

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%% The matching PH is exponential %%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if( (abs(m2-2) < epsilon^2) && (abs(m3-3) < epsilon^2) ) {

  tau <- 1
  T <- - 1 / momentsofG[1]
  xi <- 1 / momentsofG[1]
  return(list(tau=tau, T=T, xi=xi))

}

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%% The matching PH is 2-phase PH %%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if( isPH2(m2,m3) ) {
  return(matching3PH2(momentsofG[1], m2, m3, epsilon=epsilon))
}

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%% The matching PH is EC (Erlang-Coxian) %%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if( m3 >= 2*m2 - 1 ) {
  return(matching3EC(momentsofG[1], m2, m3, epsilon=epsilon))
}

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%% The matching PH is "p * EC + (1-p) * EXP" %%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k <- floor((2*m2-m3)/(m3-m2) + epsilon^2)

if( m3 >= ((k+1)*m2 + (k+4)) / (2*(k+2)) * m2 ) {

  y <- (2-m2) / (4*(3/2-m3/m2));
  p <- (2-m2)^2 / ((2-m2)^2+4*(2*m2-1-m3))
  m2X <- 2 * y
  m3X <- 2 * m2X - 1

  phX <- matching3EC(1, m2X, m3X, epsilon=epsilon)
  ph <- mixtureofPH(phX$tau, phX$T, phX$xi, 1, matrix(-1/y,1,1), 1/y, p)

  ## normalize
  mean <- sum(solve(a=-t(ph$T), ph$tau))
  T <- ph$T * mean / momentsofG[1]
  xi <- ph$xi * mean / momentsofG[1]

  return(list(tau=ph$tau, T=T, xi=xi))
}

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%% The matching PH is "EC + EXP"   %%
##%% EC has mass probability at zero %%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if( abs(m2-2) < epsilon^2 ) {

  y <- (m3-2*(k+3)/(k+2)) / (3-m3)
  m2X <- 2*(1+y)

} else {

  A <- m2 * ((m3-3) - 2*(k+3)/(k+2)*(m2-2))
  B <- m2 * sqrt((m3-3)^2 + 8*(k+3)/(k+2)*(m2-2)*(3/2-m3/m2))
  C <- 2 * (k+3)/(k+2) * (m2-2)^2
  y <- (A+B) / C
  m2X <- (1+y) * (m2*(1+y)-2*y)

##  %% This is a magic spell for Matlab, who doesn't like unnecessary precision...
##%%  m2X = round(m2X * 1000) / 1000;

}
  
m3X <- (k+3)/(k+2) * m2X

phX <- matching3EC(1, m2X, m3X, epsilon=epsilon)

ph <- convolutionofPH(phX$tau, phX$T, phX$xi, 1, matrix(-1/y,1,1), 1/y)

## normalize
mean <- sum(solve(a=-t(ph$T), ph$tau))
T <- ph$T * mean / momentsofG[1]
xi <- ph$xi * mean / momentsofG[1]

return(list(tau=ph$tau, T=T, xi=xi))

}

