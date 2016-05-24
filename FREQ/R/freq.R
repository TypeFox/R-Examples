freq <-
function(fi){

# Data construction
i <- c (1:length (fi))
ifi <- i * fi

xdata <- cbind (i,fi, ifi)

###############################
##### PSE Geometric Distribution ####
###############################

paramgeom <- c("M(t+1)", "q", "p", "N", "SE", "lower_symm_CI", "upper_symm_CI", "lower_asymm_CI", "upper_asymm_CI")

GeomDist <- function (i, fi, ifi) {
M <- sum (fi)# number of captured individuals
q <- (sum (ifi)-sum (fi))/(sum (ifi) - 1)# parameter of the distribution
p <- 1-q# 1-q=p
N <- (sum (fi))/(q)# estimated population size
var <- (N*(N-sum (fi))^2)/(sum (fi))^2# variance of N
se <- sqrt (var)# SE of N
l <- N-1.96*se# lower symmetric confidence intervall
u <- N+1.96*se# upper symmetric confidence intervall
f0 <- N-M# individuals not captured
C <- exp (1.96 * sqrt(log(1+(var/f0^2))))# log-normal distribution of the individuals not captured
al <- M + f0/C# lower asymmetric confidence intervall
au <- M + f0 * C# upper asymmetric confidence intervall
result <- c(M,
round (q, digits = 4),
round (p, digits = 4), 
round (N, digits = 2), 
round (se, digits = 2), 
round (l, digits = 2), 
round (u, digits = 2),
round (al, digits = 2),
round (au, digits = 2))
resmat <- data.frame (cbind (paramgeom, result))
return (resmat)
}
res_geomdist <- as.vector(GeomDist (i,fi,ifi) [c(1:2,4:9),2])

# expected frequencies out of the geometric distribution
expected_frequencies <- vector ("numeric", length (fi))
for (j in 1:length (fi)) {
q <- (sum (ifi)-sum (fi))/(sum (ifi) - 1)
p <- 1-q
n <- (sum (fi))/(q)
expected_frequencies [j] <- n*p*q^(i [j])
expection_geom <- cbind (i,fi,expected_frequencies)
}
expection_geom

GeomDist (i,fi,ifi)


###############################
##### PSE Poisson Distribution ####
###############################

# precondition: all 7 basic assumptions are true, especially equal catchability

# calculation of lambda, average number of catches per individual incl. all individuals not captured
# calculation of x as average number of catches per individual CAPTURED
# x is denoted as a in the following

a <- (sum(ifi))/(sum(fi))

# a = la/(1-e^(-la))
# f(la) = 0 = a-a*e^(-la)-la
# for an overview in R:
# curve (a-a*exp(-x)-x)
# than you can simply calculate the root with
# la <- uniroot (function (x) a-a*exp(-x)-x, c(0.001,5))$root

# or calculate it using Taylor and
# building the first derivation
# f'(la) =a*e^(-la)-1 
# building the Taylor series
# f~(la) = f(la) + f'(la)/1! * (la-la0) = f'(la) *  (la-la0)f(la)
# assume that the Taylor series is perfect, then f~(la) = f(la)
# only true, if la = la0
# la0 is unknown, so we assume a la and calculate a new la0, and due to the monotonicity of the function this aspires against the root
# la0 = la - f(la)/f'(la), and therefore
# la(n+1) = la(n) - f(la)/f'(la)
# then the R-code would be the following:

# x <- c(1:5000)/1000
# for (s in 1:length(x)) {
#if (x[s+1] > x[s] - ((a-a*exp(-x[s])-x[s])/(a*exp(-x[s]-1)))) {
#next
#}
#else {
#la <- x[s-1]
#break
#}
# }

# but I would prefer the Newton procedure, due to the interval independency
# then f(la) = 0 = a-la/(1-exp(-la))
# because the function is not defined at 1-exp(-la) = 0
# and f'(la) = (la * exp(-la))/(1-exp(-la))^2 - 1/(1-exp(-la))
# la(n+1) = la(n) - f(la)/f'(la)

x <- vector ("numeric")
for (s in 1:10) {
x[1] <- 0.01
x[s+1] <- x[s] - ((a-(x[s])/(1-exp(-x[s])))/((x[s] * exp(-x[s]))/(1-exp(-x[s]))^2 - 1/(1-exp(-x[s]))))
la <- x[s+1]
}

# Poisson distribution

parampois <- c("M(t+1)", "av_number_catches/ind", "N", "SE", "lower_symm_CI", "upper_symm_CI", "lower_asymm_CI", "upper_asymm_CI")

PoisDist <- function (i, fi, ifi, la) {
M <- sum (fi)# number of captured individuals
a <- (sum(ifi))/(sum(fi))# average number of catches per individual captured
N <- (sum (ifi))/(la)# estimated population size
var <- ((N^3)*(1-exp(-(sum(ifi))/N))^3)/((sum (fi))*(sum(ifi))*(1-exp(-(sum(ifi))/N)+((sum(fi))/N)*exp(-(sum(ifi))/N)))# variance of N
se <- sqrt (var)# SE of N
l <- N-1.96*se# lower symmetric confidence intervall
u <- N+1.96*se# upper symmetric confidence intervall
f0 <- N-M# individuals not captured
C <- exp (1.96 * sqrt(log(1+(var/f0^2))))# log-normal distribution of the individuals not captured
al <- M + f0/C# lower asymmetric confidence intervall
au <- M + f0 * C# upper asymmetric confidence intervall
result <- c(M,
round (la, digits = 4),
round (N, digits = 2), 
round (se, digits = 2), 
round (l, digits = 2), 
round (u, digits = 2),
round (al, digits = 2),
round (au, digits = 2))
resmat <- data.frame (cbind (parampois, result))
return (resmat)
}
res_poisdist <- as.vector(PoisDist (i,fi,ifi,la) [,2])

# expected frequencies out of the Poisson distribution
expected_freqPois <- vector ("numeric", length (fi))
for (j in 1:length (fi)) {
expected_freqPois [j] <- (exp(-la))*((sum(fi))/(1-exp(-la)))*((la^j)/factorial(j))
expection_pois <- cbind (i,fi,expected_freqPois)
}
expection_pois

PoisDist(i,fi,ifi,la)

###############################
##### PSE Negative-Binomial Distribution ####
###############################

# assumption of equal catchability is negligible

paramnegbinom <- c("M(t+1)", "av_number_catches/ind", "N", "SE", "lower_symm_CI", "upper_symm_CI", "lower_asymm_CI", "upper_asymm_CI")

NegBinomDist <- function (i, fi, ifi) {
M <- sum (fi)# number of captured individuals
a <- (sum(ifi))/(sum(fi))# average number of catches per individual captured
ssq <- (1/((sum(fi))-1))*((sum(i^2*fi))-((sum(ifi))^2/(sum(fi))))# s squared as variance of the capture frequence of all individuals captured at least once
la <- a-(ssq*fi[1])/(a*(sum(fi)-fi[1]))# lambda as average catches per individual inkl.  individuals not captured
N <- (sum (ifi))/(la)# estimated population size
var <- ((N^3)*(1-exp(-(sum(ifi))/N))^3)/((sum (fi))*(sum(ifi))*(1-exp(-(sum(ifi))/N)+((sum(fi))/N)*exp(-(sum(ifi))/N)))# variance of N
se <- sqrt (var)# SE of N
l <- N-1.96*se# lower symmetric confidence intervall
u <- N+1.96*se# upper symmetric confidence intervall
f0 <- N-M# individuals not captured
C <- exp (1.96 * sqrt(log(1+(var/f0^2))))# log-normal distribution of the individuals not captured
al <- M + f0/C# lower asymmetric confidence intervall
au <- M + f0 * C# upper asymmetric confidence intervall
result <- c(M,
round (la, digits = 4),
round (N, digits = 2), 
round (se, digits = 2), 
round (l, digits = 2), 
round (u, digits = 2),
round (al, digits = 2),
round (au, digits = 2))
resmat <- data.frame (cbind (paramnegbinom, result))
return (resmat)
}
res_negbinomdist <- as.vector(NegBinomDist (i,fi,ifi) [,2])

# expected frequencies out of the Negative Binomial distribution
expected_freqNeg <- vector ("numeric", length (fi))
for (j in 1:length (fi)) {
a <- (sum(ifi))/(sum(fi))
ssq <- (1/((sum(fi))-1))*((sum(i^2*fi))-((sum(ifi))^2/(sum(fi))))
w <- (a/ssq)*(1-(fi[1])/(sum(fi)))
k <- (w*a-(fi[1])/(sum(fi)))*(1-w)^(-1)
expected_freqNeg [j] <- (sum(fi))*((w^k*(1-w)^j)/(1-w^k))*(factorial(k+j-1)/(factorial(j)*factorial(k-1)))
expection_Neg <- cbind (i,fi,expected_freqNeg)
}
expection_Neg

NegBinomDist(i,fi,ifi)

###############################
##### PSE Truncated Geometric Distribution ####
###############################

# table for q, Seber 1982, Tab A4 (582)
tabtrunc <- matrix (NA,nrow=26,ncol=11)
s <- c(1:26)
rownames (tabtrunc) = s
P <- c(0.001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.999)
colnames (tabtrunc) = P
tabtrunc [1,] <- P
tabtrunc [2,] <- c(1.5,1.474,1.444,1.412,1.375,1.333,1.286,1.231,1.167,1.091,1.001)
tabtrunc [3,] <- c(1.999,1.93,1.852,1.767,1.673,1.571,1.462,1.345,1.226,1.108,1.001)
tabtrunc [4,] <- c(2.499,2.369,2.225,2.069,1.904,1.733,1.562,1.396,1.244,1.111,1.001)
tabtrunc [5,] <- c(2.998,2.79,2.563,2.323,2.078,1.839,1.615,1.416,1.248,1.111,1.001)
tabtrunc [6,] <- c(3.497,3.195,2.868,2.533,2.206,1.905,1.642,1.424,1.25,1.111,1.001)
tabtrunc [7,] <- c(3.996,3.582,3.142,2.705,2.298,1.945,1.655,1.427,1.25,1.111,1.001)
tabtrunc [8,] <- c(4.495,3.953,3.387,2.844,2.363,1.969,1.661,1.428,1.25,1.111,1.001)
tabtrunc [9,] <- c(4.993,4.308,3.605,2.955,2.408,1.982,1.664,1.428,1.25,1.111,1.001)
tabtrunc [10,] <- c(5.492,4.647,3.797,3.043,2.439,1.99,1.666,1.429,1.25,1.111,1.001)
tabtrunc [11,] <- c(5.99,4.969,3.966,3.111,2.46,1.995,1.666,1.429,1.25,1.111,1.001)
tabtrunc [12,] <- c(6.488,5.277,4.115,3.165,2.474,1.997,1.666,1.429,1.25,1.111,1.001)
tabtrunc [13,] <- c(6.986,5.569,4.244,3.206,2.483,1.998,1.667,1.429,1.25,1.111,1.001)
tabtrunc [14,] <- c(7.484,5.847,4.356,3.238,2.489,1.999,1.667,1.429,1.25,1.111,1.001)
tabtrunc [15,] <- c(7.981,6.111,4.453,3.262,2.493,2,1.667,1.429,1.25,1.111,1.001)
tabtrunc [16,] <- c(8.479,6.360,4.537,3.280,2.495,2,1.667,1.429,1.25,1.111,1.001)
tabtrunc [17,] <- c(8.976,6.597,4.608,3.294,2.497,2,1.667,1.429,1.25,1.111,1.001)
tabtrunc [18,] <- c(9.473,6.821,4.670,3.304,2.498,2,1.667,1.429,1.25,1.111,1.001)
tabtrunc [19,] <- c(9.970,7.033,4.722,3.312,2.499,2,1.667,1.429,1.25,1.111,1.001)
tabtrunc [20,] <- c(10.467,7.232,4.767,3.317,2.499,2,1.667,1.429,1.25,1.111,1.001)
tabtrunc [21,] <- c(10.963,7.429,4.805,3.322,2.5,2,1.667,1.429,1.25,1.111,1.001)
tabtrunc [22,] <- c(11.460,7.597,4.836,3.325,2.5,2,1.667,1.429,1.25,1.111,1.001)
tabtrunc [23,] <- c(11.956,7.763,4.863,3.327,2.5,2,1.667,1.429,1.25,1.111,1.001)
tabtrunc [24,] <- c(12.425,7.92,4.886,3.329,2.5,2,1.667,1.429,1.25,1.111,1.001)
tabtrunc [25,] <- c(12.948,8.066,4.905,3.330,2.5,2,1.667,1.429,1.25,1.111,1.001)
tabtrunc [26,] <- c(13.444,8.204,4.921,3.331,2.5,2,1.667,1.429,1.25,1.111,1.001)

paramgeom <- c("M(t+1)", "q", "p", "N", "SE", "lower_symm_CI", "upper_symm_CI", "lower_asymm_CI", "upper_asymm_CI")

TruncGeomDist <- function (i, fi, ifi) {
M <- sum (fi)# number of captured individuals
imax <- length(i)# number of secondary periods (occasions)
a <- (sum(ifi))/(sum(fi))# average number of catches per individual captured
for (j in 1:11) {
if (tabtrunc[imax,j] >= a)
{
plower <- (tabtrunc[1,j])
xupper <- (tabtrunc[imax,j])
}
else
{
pupper <- (tabtrunc[1,j])
xlower <- (tabtrunc[imax,j])
break
}
}
p <- plower+((xupper-a)/(xupper-xlower))*(pupper-plower)# exact p
q <- 1-p# q=1-p, exact q
N <- (sum (fi))/(q)# estimated population size
var <- (N*(N-sum (fi))^2)/(sum (fi))^2# variance of N
se <- sqrt (var)# SE of N
l <- N-1.96*se# lower symmetric confidence intervall
u <- N+1.96*se# upper symmetric confidence intervall
f0 <- N-M# individuals not captured
C <- exp (1.96 * sqrt(log(1+(var/f0^2))))# log-normal distribution of the individuals not captured
al <- M + f0/C# lower asymmetric confidence intervall
au <- M + f0 * C# upper asymmetric confidence intervall
result <- c(M,
round (q, digits = 4),
round (p, digits = 4), 
round (N, digits = 2), 
round (se, digits = 2), 
round (l, digits = 2), 
round (u, digits = 2),
round (al, digits = 2),
round (au, digits = 2))
resmat <- data.frame (cbind (paramgeom, result))
return (resmat)
}
res_truncgeomdist <- as.vector(TruncGeomDist (i,fi,ifi) [c(1:2,4:9),2])

# expected frequencies out of the truncated geometric distribution
expected_freqTrunc <- vector ("numeric", length (fi))
for (r in 1:length (fi)) {
imax <- length(i)
a <- (sum(ifi))/(sum(fi))
for (j in 1:11) {
if (tabtrunc[imax,j] >= a)
{
plower <- (tabtrunc[1,j])
xupper <- (tabtrunc[imax,j])
}
else
{
pupper <- (tabtrunc[1,j])
xlower <- (tabtrunc[imax,j])
break
}
}
p <- plower+((xupper-a)/(xupper-xlower))*(pupper-plower)
q <- 1-p
N <- (sum (fi))/(q)
expected_freqTrunc [r] <- n*p*q^(i [r])
expection_trunc <- cbind (i,fi,expected_freqTrunc)
}
expection_trunc

TruncGeomDist (i,fi,ifi)

###############################
#### complete output ##################
###############################

# all expected and meassured values for each occasion
res_expected <- matrix (NA,length(i),5)
colnames (res_expected) <- c("original data", "Geometric Distribution", "Truncated Geometric Distribution", "Poisson Distribution", "Negative Binomial Distribution")
rownames (res_expected) <- i
res_expected [,1] <- fi
res_expected [,2] <- expection_geom [,3]
res_expected [,3] <- expection_trunc [,3]
res_expected [,4] <- expection_pois [,3]
res_expected [,5] <- expection_Neg [,3]

# the calculations and estimations for the population size
res_estimates <- matrix(NA, 8,4)
res_estimates <- as.data.frame(res_estimates)
colnames (res_estimates) <- c("Geometric Distribution", "Truncated Geometric Distribution", "Poisson Distribution", "Negative Binomial Distribution")
rownames (res_estimates) <- c("number of individuals captured","parameter q or lambda", "population size N", "standard error", "lower symmetric 95%-CI", "upper symmetric 95%-CI", "lower asymmetric 95%-CI", "upper asymmetric 95%-CI")
res_estimates [,1] <- res_geomdist
res_estimates [,2] <- res_truncgeomdist
res_estimates [,3] <- res_poisdist
res_estimates [,4] <- res_negbinomdist

message3 <- paste ("If you want to continue with further data, recall freq(c(...)).\n",
"\n",
"\n",
"FREQ 2013",
"\n",
"written by Annegret Grimm :-)\n",
"\n",sep="")


message1 <- paste ("All meassured and expected values:\n")
message2 <- paste ("All estimated values:\n")
cat (message1)
print(res_expected)
cat(paste("\n","\n"))
cat (message2)
print(res_estimates)
cat(paste("\n","\n"))
cat (message3)

}
