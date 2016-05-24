
#----------ExactBinomial.R-------------------------------------------------------------------

## VERSION OF FEB/17 2015

# -------------------------------------------------------------------------------------------
# Function produces critical value for the continuous Sequential Binomial MaxSPRT
# -------------------------------------------------------------------------------------------

CV.Binomial<- function(N,alpha=0.05,M=1,z=1){

MinCases<- M

# N = maximum length of surveillance defined in terms of the total number of adverse events
# alpha = desired alpha level
# MinCases = The minimum number of cases for which a signal is allowed to occur, default=1
# z = matching ratio between exposed and unexposed cases 

if(is.numeric(N)==FALSE){stop("The maximum length of surveillance, 'N', must be a positive integer.",call. =FALSE)} 
if(is.numeric(M)==FALSE){stop("The minimum number of cases, 'M', must be a positive integer.",call. =FALSE)} 
if(is.numeric(alpha)==FALSE){stop("The significance level, 'alpha', must be a number greater than zero and smaller than 0.5.",call. =FALSE)}
if(is.numeric(z)==FALSE){stop("The matching ratio, 'z', must be a positive number.",call. =FALSE)}

if(N!=round(N)){stop("The maximum length of surveillance, 'N', must be a positive integer.",call. =FALSE)}
if(MinCases!=round(MinCases)){stop("The minimum number of cases, 'M', must be a positive integer.",call. =FALSE)}

if(N<=0){stop("The maximum length of surveillance, 'N', must be a positive integer.",call. =FALSE)}
if(M<=0){stop("The minimum number of cases, 'M', must be a positive integer.",call. =FALSE)}
if(z<=0){stop("The matching ratio, 'z', must be a positive number.",call. =FALSE)}
if(M>N){
pst<- 1/(1+z)
if(1-pbinom(MinCases-1,N,pst)<alpha){Mr<- MinCases-1;while(1-pbinom(Mr-1,N,pst)<alpha&Mr>1){Mr<- Mr-1};if(Mr==0){Mr<- 1}
stop(c("For this maximum length of surveillance, the minimum number of cases, 'M', must be a positive integer smaller than or equal to"," ",Mr),call. =FALSE)}
       }



if(alpha<=0|alpha>0.5){stop("The significance level, 'alpha', must be a number greater than zero and smaller than 0.5.",call. =FALSE)}


pst<- 1/(1+z)
if(1-pbinom(MinCases-1,N,pst)<alpha){Mr<- MinCases-1;while(1-pbinom(Mr-1,N,pst)<alpha&Mr>1){Mr<- Mr-1};if(Mr==0){Mr<- 1}
stop(c("For this maximum length of surveillance, the minimum number of cases, 'M', must be a positive integer smaller than or equal to"," ",Mr),call. =FALSE)
                                    }

if(1-pbinom(N-1,N,pst)>alpha){
Nr<- N
while(1-pbinom(Nr-1,Nr,pst)>alpha){Nr<- Nr+1}
stop(c("There is no solution with Type I error probability smaller than"," ",alpha,". Use 'N' of at least"," ",Nr,"."),call. =FALSE)
                             }


# Function that calculates the LLR for a given observed (cc) cases from exposed period and n total of cases
#-------------------------------------------------------------------------------------------

LLR <- function(cc,n,z){

       if(cc==n){x = n*log(1+z)}else{
         if(z*cc/(n-cc)<=1){x=0}else{
	       x = cc*log(cc/n)+(n-cc)*log((n-cc)/n)-cc*log(1/(z+1))-(n-cc)*log(z/(z+1))
                                    }
                                  } 	
      	x
	}

# Auxiliar code that calculates the type I error probability for a given cv
#### Function to calculate the Type I error 

Erro_I<- function(cv,RR=1)
{


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
									# starting probabilities are all set to zero's

# Calculating the first row in the p[][] matrix for which there is a chance to reject H0
# --------------------------------------------------------------------------------------

for(s in 1:absorb[1]){ p[1,s]=dbinom(s-1,1,ps)}		# Probability of having s-1 cases at time mu[1]
p[1,absorb[1]+1]=1-pbinom(absorb[1]-1,1,ps)			# probability of rejecting H0 at time mu[1]


if(N>1){
i<- 1
alphai<- 0
while(i<N&alphai<alpha){
i<- i+1

       p[i,1:absorb[i]]<- func_aux1(i) # Calculates the standard p[][] cell values
       p[i,absorb[i]+1]<- func_aux3(i) # Calculates the diagonal absorbing states where H0 is rejected
	

                             alphai<- alphai+ p[i,absorb[i]+1] 
                 
                       } # end for i	
}

alpha=0
for(i in 1:N){ alpha=alpha+p[i,absorb[i]+1]}	
return(alpha)

}
## Ending fucntion type I error

## Finding the critical value by using bisection strategy

omega<- matrix(0,nrow=1)
for(i in 1:N){j<- matrix(seq(1,i,1),ncol=1) ; omega<- cbind(omega,matrix(apply(j,1,LLR,i,z),nrow=1))}

omega<- omega[order(omega)]
begin<- sum(omega==0)+1
omega<- omega[begin:length(omega)]

i1<- 1 ;  i2<- length(omega) ; im<- round((i1+i2)/2)
aux_extrem<- 0
while(i2-i1>1){
 error1<- Erro_I(omega[im])
      if(error1>alpha){i1<- im}else{i2<- im; resold<- error1; ir<- im;aux_extrem<-1}
 im<- round((i1+i2)/2)
              }

if(aux_extrem==0){cv<- omega[length(omega)];error1<- Erro_I(cv)}else{error1<- resold; cv<- omega[i1] }
 
if(round(cv,5)<=cv){cv<- round(cv+(9*10^{-6}),5)}else{cv<- round(cv,5)}
Type_I_Error<- error1
res<- list(cv,Type_I_Error)
names(res)<- c("cv","Type_I_Error")
return(res)

} #end function CV.Binomial







