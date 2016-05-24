
#----------ExactBinomial.R-------------------------------------------------------------------

# -------------------------------------------------------------------------
# Function produces power and expected signal time for the Group Sequential Binomial MaxSPRT
# -------------------------------------------------------------------------

Performance.G.Binomial <- function(N,M=1,cv,z=1,RR=2,GroupSizes){
Groups<- GroupSizes
alpha<- 0.05
MinCases<- M
# N = maximum length of surveillance defined in terms of the total number of adverse events
# alpha = desired alpha level
# MinCases = The minimum number of cases for which a signal is allowed to occur, default=1
# z = matching ratio between exposed and unexposed cases    
# Groups: Vector with the number of adverse events (exposed+unexposed) between two looks at the data, i.e, irregular group sizes. Important: Must sums up N
# RR= Relative risk 

if(length(Groups)==1){
if(is.numeric(Groups)==FALSE){stop("'Groups' must be an integer smaller than or equal to 'N'.",call. =FALSE)}
if(Groups==0){stop("'N' must be a positive integer'.",call. =FALSE)}
if(N/Groups!=round(N/Groups)){stop("'N' must be a multiple of 'Groups'.",call. =FALSE)}
if(Groups>N){stop("'N' must be a multiple of 'Groups'.",call.=FALSE)}
}

if(length(Groups)>1){
if(sum(is.numeric(Groups))==0){stop("'Groups' must be a vector of positive integers.",call. =FALSE)}else{
if(sum(Groups)!=N){stop("'Groups' must sums up 'N'.",call. =FALSE)}
if(sum(Groups<0)>0){stop("The vector 'Groups' must contain only positive integers.",call. =FALSE)}
}
}
if((N<=0|is.numeric(N)==FALSE)){stop("'N' must be a positive integer.",call. =FALSE)}
if(alpha<=0|alpha>0.5||is.numeric(alpha)==FALSE){stop("'alpha' must be a number greater than zero and smaller than 0.5.",call. =FALSE)}
if(z<0||is.numeric(z)==FALSE){stop("'z' must be a number greater than zero.",call. =FALSE)}
if(MinCases>N||is.numeric(MinCases)==FALSE){stop("'M' must be an integer smaller than or equal to N.",call. =FALSE)}
if(MinCases<1){stop("'M' must be an integer greater than zero.",call. =FALSE)}
if(MinCases!=round(MinCases)){stop("'M' must be an integer.",call. =FALSE)}
if(cv<0||is.numeric(cv)==FALSE|length(cv)>1){stop("cv must be a number greater than zero.",call. =FALSE)}
if(RR<1||is.numeric(RR)==FALSE){stop("RR must be a number greater than or equal to 1.",call. =FALSE)}
if(round(N)!=N){stop("'N' must be a positive integer.",call. =FALSE)}
if(round(M)!=M){stop("'M' must be a positive integer.",call. =FALSE)}

pst<- 1/(1+z)


if(length(Groups)>1){an<- Groups%*%(upper.tri(matrix(0,length(Groups),length(Groups)),diag=T)*1)}else{an<- seq(Groups,N,Groups)
                                                                                                      if(max(an)<N){an<- c(an,N)}
                                                                                                     }

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

# absorb[i]: number of acumulated cases (from the exposed period) needed to reject the null at the i-th adverse event 
# aux[i]: has zero entree if LLR(absorb[i],i,z)< cv or has 1 entree otherwise
# an[kk]:  order of the adverse event associated to the kk-th test

absorb = rep(0,N)		# Contains the number of events needed at time mu[i] in order to reject H0
aux<- rep(0,N)

for(i in 1:N){
      if(sum(an==i)>0){
	while( LLR(absorb[i],i,z)<cv &absorb[i]<i){ 
		absorb[i]=absorb[i]+1               }
             if(LLR(absorb[i],i,z)>=cv){aux[i]<- 1}
                      }else{absorb[i]<- i+1}
             }

if(MinCases>1){
aux[1:(MinCases-1)]<- 0
              }


for(i in 1:N){if(absorb[i]<MinCases&i>=MinCases){absorb[i]<- MinCases};if(absorb[i]<MinCases&i<MinCases|aux[i]==0){absorb[i]<- i+1}}

uc<- absorb-1

ps<- RR/(RR+z)

# Auxiliar functions to run the binomial Markov Chain in a fast way:
func_aux2<- function(j,i){ k<- seq(1,uc[i-1]+1); return(sum(p[i-1,k]*dbinom(j-k,1,ps)))} ; func_aux3<- function(i){ k<- seq(1,uc[i-1]+1); return(sum(p[i-1,k]*(1-pbinom(absorb[i]-k,1,ps))))}
func_aux1<- function(i){ j<- matrix(seq(1,absorb[i]),ncol=1) ; return(apply(j,1,func_aux2,i))}

p<- matrix(0,N,N+2)    	# p[i,j] is the probability of having j-1 cases at time mu[i]
									# starting probabilities are all set to zero's


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

# Sums up the probabilities of absorbing states when a signal occurs, to get the alpha level
# ------------------------------------------------------------------------------------------

power=0
time<- 0
for(i in 1:N){ 
                                    power=power+p[i,absorb[i]+1]
                                    time<- time + i*p[i,absorb[i]+1]
                                          
             }					
	signaltime<- time/power
      surveillancetime<- time + N*(1-power)

result<- list(power,signaltime,surveillancetime)
names(result)<- c("Power","ESignalTime","ESampleSize")
return(result)


} #end function Performance.G.Binomial

