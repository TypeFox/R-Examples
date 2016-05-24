Performance.G.Poisson <-
function(SampleSize,cv,Looks,M=1,RR=2){

LLR <- function(cc,uu) {
	if(cc<=uu) x=0
	if(cc>uu) x = (uu-cc) + cc*log(cc/uu)
	x
	}
CV<- cv
T<- SampleSize
L<- T
MinCases<- M
Group<- T/Looks

if(is.numeric(MinCases)==FALSE|MinCases!=round(MinCases)|MinCases<1){stop("M must be a positive integer.",call. =FALSE)}
if(T<=0|is.numeric(T)==FALSE){stop("SampleSize must be a positive number.",call. =FALSE)}
if(is.numeric(Looks)==FALSE|Looks<1|Looks!=round(Looks)){stop("Looks must be a positive integer.",call. =FALSE)}


# T = maximum length of surveillance, defined in terms of expected counts under H0
# CV = Critical value
# Group = Time between looks in a group sequential trial, must be greater than 0
# Relative risk


imax = ceiling(T/Group)					# The number of tests performed, including final time T.
mu = seq(length=imax,from=Group,by=Group)		# An array of the expected counts at each of the tests
								# mu[i] is the commulative expected count at the i'th test
mu[imax]=T							# Sets the expected count at the last test to equal T



absorb = seq(length=imax,from=0,by=0)		# Contains the number of events needed at time mu[i] in order to reject H0
for(i in 1:imax)
	while( LLR(absorb[i],mu[i])<CV ) 
		absorb[i]=absorb[i]+1;

mu=mu*RR

absorb[absorb<MinCases]<- MinCases

p = seq(length=imax*(absorb[imax]+1), from=0, by=0)	# p[i,j] is the probability of having j-1 cases at time mu[i]
									# starting probabilities are all set to zero's
dim(p) = c(imax,absorb[imax]+1)				# i in 1:imax-1 is the rows and j in 1:imax is the column

# Calculating the first row in the p[][] matrix for which there is a chance to reject H0
# --------------------------------------------------------------------------------------

for(s in 1:absorb[1]) p[1,s]=dpois(s-1,mu[1])		# Probability of having s-1 cases at time mu[1]
p[1,absorb[1]+1]=1-ppois(absorb[1]-1,mu[1])			# probability of rejecting H0 at time mu[1]


# Calculating the remaining rows in the p[][] matix
# -------------------------------------------------
if(imax>1){
for(i in 2:imax) {
	for(j in 1:absorb[i])					# This loop calculates the p[][] matix, one column at a time, from left to right
		for(k in 1:min(j,absorb[i-1]))
			p[i,j]=p[i,j]+p[i-1,k]*dpois(j-k,mu[i]-mu[i-1])	# Calculates the standard p[][] cell values
	for(k in 1:absorb[i-1]) p[i,absorb[i]+1]=p[i,absorb[i]+1]+p[i-1,k]*(1-ppois(absorb[i]-k,mu[i]-mu[i-1]))
} # end for i	
          }

# Sums up the probabilities of absorbing states when a signal occurs, to get the alpha level
# ------------------------------------------------------------------------------------------

power=0
for(i in 1:imax) power=power+p[i,absorb[i]+1]					
time=0
for(i in 1:imax) time=time+i*Group*p[i,absorb[i]+1]
signaltime=time/power					
surveillancetime<- time+(1-power)*L

result<- list(power,signaltime,surveillancetime)
names(result)<- c("Power","ESignalTime","ESampleSize")
return(result)


}


