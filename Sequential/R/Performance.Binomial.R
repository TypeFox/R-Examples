


#----------ExactBinomial.R-------------------------------------------------------------------

# Version of Feb/2015

# -------------------------------------------------------------------------
# Function produces power and expected signal time for the continuous Sequential Binomial MaxSPRT
# -------------------------------------------------------------------------

Performance.Binomial<- function(N,M=1,cv,z=1,RR=2){
alpha<- 0.05
MinCases<- M
# N = maximum length of surveillance defined in terms of the total number of adverse events
# alpha = desired alpha level
# MinCases = The minimum number of cases for which a signal is allowed to occur, default=1
# z = matching ratio between exposed and unexposed cases
# RR is the relative risk

if(is.numeric(N)==FALSE|N<=0){stop("N must be a positive integer.",call. =FALSE)}
if(is.numeric(M)==FALSE|M<=0){stop("M must be a positive integer smaller than or equal to 'N'.",call. =FALSE)}
if(is.numeric(alpha)==FALSE|alpha<=0|alpha>0.5){stop("alpha must be a number in the '(0,0.5]' interval.",call. =FALSE)}
if(round(N)!=N){stop("N must be a positive integer.",call. =FALSE)}
if(round(M)!=M){stop("M must be a positive integer.",call. =FALSE)}
if(is.numeric(cv)==FALSE){stop("cv must be a positive number.",call. =FALSE)}
if(is.numeric(z)==FALSE|z<=0){stop("z must be a number greater than zero.",call. =FALSE)}
if(is.numeric(RR)==FALSE|RR<=0){stop("RR must be a number greater than zero.",call. =FALSE)}


if(MinCases>N){stop("'MinCases' must be an integer smaller than or equal to N.",call. =FALSE)}
if(MinCases<1){stop("'MinCases' must be an integer greater than zero.",call. =FALSE)}
if(MinCases!=round(MinCases)){stop("'MinCases' must be an integer.",call. =FALSE)}

if(alpha<=0|alpha>0.5){stop("alpha must be a number greater than zero and smaller than 0.5.",call. =FALSE)}

if(N<=0){stop("N must be a positive integer.",call. =FALSE)}
if(N!=round(N)){stop("'N' must be an integer.",call. =FALSE)}

if(z<=0){stop("z must be greater than zero.",call. =FALSE)} 
if(RR<1){stop("RR must be greater than or equal to 1.",call. =FALSE)}
if(cv<=0){stop("cv must be a positive number.",call. =FALSE)}

pst<- 1/(1+z)

# Function that calculates the LLR for a given observed (c) and expected (u) number of cases
#-------------------------------------------------------------------------------------------

LLR <- function(cc,n,z){

       if(cc==n){x = n*log(1+z)}else{
         if(z*cc/(n-cc)<=1){x=0}else{
	       x = cc*log(cc/n)+(n-cc)*log((n-cc)/n)-cc*log(1/(z+1))-(n-cc)*log(z/(z+1))
                                    }
                                  } 	
      	x
	}

# Preparing the information of absorbing states
#-------------------------------------------------------------------------------------------
 
absorb = rep(0,N)		# Contains the number of events needed at time mu[i] in order to reject H0
aux<- rep(0,N)

for(i in 1:N){
	while( LLR(absorb[i],i,z)<cv &absorb[i]<i){ 
		absorb[i]=absorb[i]+1               }
             if(LLR(absorb[i],i,z)>=cv){aux[i]<- 1}
             }
if(MinCases>1){
aux[1:(MinCases-1)]<- 0
              }

absorb[aux==0]<- absorb[absorb[aux==0]]+1

for(i in 1:N){if(absorb[i]<MinCases&i>=MinCases){absorb[i]<- MinCases};if(absorb[i]<MinCases&i<MinCases){absorb[i]<- i+1}}


uc<- absorb-1

ps<- RR/(RR+z)

# Auxiliar functions to run the binomial Markov Chain in a fast way:
func_aux2<- function(j,i){ k<- seq(1,uc[i-1]+1); return(sum(p[i-1,k]*dbinom(j-k,1,ps)))} ; func_aux3<- function(i){ k<- seq(1,uc[i-1]+1); return(sum(p[i-1,k]*(1-pbinom(absorb[i]-k,1,ps))))}
func_aux1<- function(i){ j<- matrix(seq(1,absorb[i]),ncol=1) ; return(apply(j,1,func_aux2,i))}

p<- matrix(0,N,N+2)	# p[i,j] is the probability of having j-1 cases at time mu[i]
                   	# p[i,j] is the probability of having j-1 cases at time mu[i]
				# starting probabilities are all set to zero's
                        # i in 1:imax-1 is the rows and j in 1:imax is the column

# Calculating the first row in the p[][] matrix for which there is a chance to reject H0
# --------------------------------------------------------------------------------------

for(s in 1:absorb[1]){ p[1,s]=dbinom(s-1,1,ps)}		# Probability of having s-1 cases at time mu[1]
p[1,absorb[1]+1]=1-pbinom(absorb[1]-1,1,ps)			# probability of rejecting H0 at time mu[1]


if(N>1){
i<- 1
while(i<N){
i<- i+1
       p[i,1:absorb[i]]<- func_aux1(i) # Calculates the standard p[][] cell values
       p[i,absorb[i]+1]<- func_aux3(i) # Calculates the diagonal absorbing states where H0 is rejected                             
                 
          } # end for i	
       }	
# Sums up the probabilities of absorbing states when a signal occurs, to get the power
# ------------------------------------------------------------------------------------------

power=0
time=0

for(i in 1:N){ power=power+p[i,absorb[i]+1];time= time+i*p[i,absorb[i]+1]}
signaltime<- time/power
surveillancetime<- time + N*(1-power)

result<- list(power,signaltime,surveillancetime)
names(result)<- c("Power","ESignaltime","ESampleSize")
return(result)


} #end function Performance.Binomial







