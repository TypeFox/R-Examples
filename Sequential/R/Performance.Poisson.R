Performance.Poisson <-
function(SampleSize,D=0,M=1,cv,RR=2){

# ------------------- INPUT VARIABLE ----------------------------------------------------------
# L = maximum length of surveillance, defined in terms of expected counts under H0
# cv = critical value
# RR = relative risk, RR=1 corresponds to H0
# M = The minimum number of cases for which a signal is allowed to occur
# D = Time < T for first look at the data, defined in terms of the expected counts under H0
# alpha = significance level
T<- SampleSize
L<- T
LLR<- cv
MinCases<- M
Late<- D


####### Tests to verify the validity of the chosen parameters

teste1<- 0

if(T<=0){teste1<- 1; out<- c("SampleSize must be > 0")}
if(teste1==0 & cv<=0){teste1<- 1; out<- c("cv must be >0")}
if(RR<1 & teste1==0){teste1<- 1; out<- c("RR must be >=1.") }
if(teste1==0 & M>100){teste1<- 1; out<- c("M must be a positive integer in the range [1,100]")}
#if(teste1==0 & T>1000){teste1<- 1; out<- c("Use SampleSize<=1000")}


if(Late>T & teste1==0){teste1<- 1; out<- c("D must be <= SampleSize") }
if(Late<0 & teste1==0){teste1<- 1; out<- c("Negative values for D does not make sense. Use 0<=D<=SampleSize.") }
if(M<1 & teste1==0){teste1<- 1; out<- c("M must be a positive integer in the range[1,100].") }

# If the parameters are incorrect in any sense, the code is interrupted and an error message is informed according to the possibilies above
#------------------------------------------------------------------------------------------------------------------------------------------
if(teste1==1){stop(out,call.=FALSE)}

#---------------------------------------------------------------------
# Function that calculates the product log through a recursive formula
#---------------------------------------------------------------------
ProdLog <- function(z) {
	x = z-z^2+1.5*z^3-(8/3)*z^4+(125/24)*z^5-(54/5)*z^6+(16807/720)*z^7
	for(i in 1:10) x = x-(x*exp(x)-z)/(exp(x)+x*exp(x))
	x
	}

c = 1:(2*T)

z = -exp(-1-LLR/c)
mu = -c * ProdLog(z) 		#The expected counts under H0 that is needed to reject the null with i number of adverse events
mtemp = c(0,mu)
mmu = diff(mtemp) 		#The marginal difference of the mu[] vector
RRmmu=mmu*RR			#The marginal difference of the mu[] vector
RRmu=mu*RR				#The expected counts under HA that is needed to reject the null with ii number of adverse events

imin=MinCases
while (mu[imin] < Late) imin=imin+1
if(imin>MinCases) { 
	mu[imin-1]=Late
	mmu[imin]=mu[imin]-Late
	} # end if 

imax=1						
while (mu[imax] < T) imax=imax+1 	# imax is the maximum number of cases that will generate a signal. 
	
if(imin<imax){
				              
# Defining the p[][] matrix
# -------------------------

p = seq(length=(imax-1)*imax, from=0, by=0)
dim(p) = c(imax-1,imax)


# Calculating the first row in the p[][] matrix for which there is a chance to reject H0
# When MinCases=1, there is no skipping, and it is the first row in the matrix (p[1][]).
# --------------------------------------------------------------------------------------

if(imin==MinCases) {
	for(s in 1:imin) p[imin,s]=dpois(s-1,RRmu[imin])		# Probability of having s-1 cases at time mu[MinCases]
	p[imin,imin+1]=1-ppois(imin-1,RRmu[imin])
	} # end if 

if(imin>MinCases) {
	for(s in 1:imin) p[imin-1,s]=dpois(s-1,RRmu[imin-1])		# Probability of having s-1 cases at time mu[imin-1], not rejecting H0
	p[imin-1,imin+1] = 1-ppois(imin-1,RRmu[imin-1])			# Probability of having s+ cases at time mu[imin-1], rejecting H0
	for(s in 1:imin) 								# Probability of having s-1 cases at time mu[imin], not rejectinh H0
		for(k in 1:s) 
			p[imin,s]=p[imin,s]+p[imin-1,k]*dpois(s-k,RRmmu[imin])	
	for(k in 1:imin) 
		p[imin,imin+1] = p[imin,imin+1] + p[imin-1,k]*(1-ppois(imin-k,RRmmu[imin]))
} # end if 


# Calculating the remaining rows in the p[][] matix
# -------------------------------------------------

if(MinCases+1<=imax-1)
for(i in (imin+1):(imax-1)) {
	for(j in 1:(i-1))								# This loop calculates the p[][] matix, one column at a time, from left to right
		for(k in 1:j) 
			p[i,j]=p[i,j]+p[i-1,k]*dpois(j-k,RRmmu[i])	# Calculates the standard p[][] cell values
	for(k in 1:(i-1))
		p[i,i]=p[i,i]+p[i-1,k]*dpois(i-k,RRmmu[i])		# Calculates the diagonal under the absorbing states, which requires a unique formula
	for(k in 1:(i-1)) 
		p[i,i+1]=p[i,i+1]+p[i-1,k]*(1-ppois(i-k,RRmmu[i]))# Calculates the diagonal absorbing states where H0 is rejected
} # end for i	
pp=0
for(k in 1:(imax-1)) pp=pp+p[imax-1,k]*(1-ppois(imax-k,(T-mu[imax-1])*RR)) #Calculates the last probability to signal before time T



# Calculating the power
# ---------------------

power=0
if(imin>MinCases) power=p[imin-1,imin+1]
for(i in imin:(imax-1)) power=power+p[i,i+1]			# Sums up the probabilities when a signal occurs, to get total power
power=power+pp

            


# Calculates the time until a signal occurs
#------------------------------------------

etime=0
if(imin==MinCases)
	for(n in imin:1000) etime=etime+dpois(n,RRmu[imin])*mu[imin]*imin/(n+1) 
if(imin>MinCases) {
	etime=etime+(1-ppois(imin-1,RRmu[imin-1]))*Late  						# (Late=mu[imin-1]) 
	etime=etime+mu[imin]*p[imin,imin+1]
	for(k in 1:imin) 
		for(n in 0:100)
			etime=etime+p[imin-1,k]*dpois(imin-k+1+n,RRmmu[imin])*mmu[imin]*(imin-k+1)/(imin-k+1+n+1)	# Adding expected times to signal, using a beta distribution
	} # end if


if(imin+1<=imax-1)
for(i in (imin+1):(imax-1)) {			
	etime=etime+mu[i-1]*p[i,i+1]
	for(k in 1:(i-1)) 
		for(n in 0:100)
			etime=etime+p[i-1,k]*dpois(i-k+1+n,RRmmu[i])*mmu[i]*(i-k+1)/(i-k+1+n+1)	# Adding expected times to signal, using a beta distribution
	} # end for i

etime=etime+mu[imax-1]*pp
margin=(T-mu[imax-1])

for(k in 1:(imax-1)) 
	for(n in 0:1000)
		etime=etime+p[imax-1,k]*dpois(imax-k+1+n,margin*RR)*margin*(imax-k+1)/(imax-k+1+n+1)	# Adding expected times to signal, using a beta distribution


# The expected time until signal, given that a signal occurs
# ----------------------------------------------------------

signaltime = etime/power


# The expected length of surveillance
# -----------------------------------

surveillancetime = etime+(1-power)*T


               }else{
                     power<- 1-ppois(imax-1,mu[imax]*RR)
                     
                     etime<- sum(dpois(seq(0,imax-1,1),mu[imax-1]*RR)*(mu[imax-1])*(1-ppois(imax-1-seq(0,imax-1,1),RR*(mu[imax]-mu[imax-1]))) )
                     etime<- etime + sum(dpois(seq(0,imax-1,1),mu[imax-1]*RR)*(mu[imax]-mu[imax-1])*(1-ppois(imax-seq(0,imax-1,1),RR*(mu[imax]-mu[imax-1])))/(RR*(mu[imax]-mu[imax-1])) )
                     signaltime<- etime + (1-ppois(imax-1,mu[imax-1]*RR))*mu[imax-1]
                    
                     surveillancetime = signaltime+(1-power)*mu[imax] 

                    } # end if(imin<imax)


# Output assigned as a vector
# ---------------------------

out=matrix(c(power,signaltime,surveillancetime),ncol=3)
colnames(out)<- c("Power","ESignalTime","ESampleSize")
return(out)

}
