


#----------ExactBinomial.R-------------------------------------------------------------------

# Version of Feb/2015

# -------------------------------------------------------------------------
# Function produces critical value for the Group Sequential Binomial MaxSPRT
# -------------------------------------------------------------------------

CV.G.Binomial <- function(N,alpha=0.05,M=1,z=1,GroupSizes){
Groups<- GroupSizes
MinCases<- M

# N = maximum length of surveillance defined in terms of the total number of adverse events
# alpha = desired alpha level
# MinCases = The minimum number of cases for which a signal is allowed to occur, default=1
# z = matching ratio between exposed and unexposed cases 
# Groups: Vector with the number of adverse events (exposed+unexposed) between two looks at the data, i.e, irregular group sizes. Important: Must sums up N
 


alpha1<- alpha
if(length(Groups)==1){
if(is.numeric(Groups)==FALSE){stop("'Groups' must be an integer smaller than or equal to 'N'.",call. =FALSE)}
if(sum(Groups<=0)>0){stop("'Groups' must be a positive integer or a vector with positive integers.",call. =FALSE)}
if(sum(Groups==round(Groups))!=length(Groups)){stop("'Groups' must be a positive integer smaller than or equal to 'N'.",call. =FALSE)}

if(Groups==0){stop("'Groups' must be a positive integer or a vector with positive integers.",call. =FALSE)}
if(N/Groups!=round(N/Groups)){stop("The maximum length of surveillance, 'N', must be a multiple of 'Groups'.",call. =FALSE)}
if(Groups>N){stop("The maximum length of surveillance, 'N', must be a multiple of 'Groups'.",call.=FALSE)}
}

if(length(Groups)>1){
if(sum(is.numeric(Groups))==0){stop("'Groups' must be a vector of positive integers.",call. =FALSE)}else{
if(is.numeric(N)==FALSE){stop("The maximum length of surveillance, 'N', must be a positive integer.",call. =FALSE)}
if(N!=round(N)){stop("The maximum length of surveillance, 'N', must be a positive integer.",call. =FALSE)}
if(sum(Groups<=0)>0){stop("'Groups' must be a positive integer or a vector with positive integers.",call. =FALSE)}
if(sum(Groups==round(Groups))!=length(Groups)){stop("'Groups' must be a positive integer or a vector with positive integers.",call. =FALSE)}
if(sum(Groups)!=N){stop("'Groups' must sum up equal to 'N'.",call. =FALSE)}
}
}



if((is.numeric(M)==FALSE)){stop("The minimum number of cases, 'M', must be a positive integer.",call. =FALSE)}
if((M<=0)){stop("The minimum number of cases, 'M', must be a positive integer.",call. =FALSE)}
if((is.numeric(alpha)==FALSE)){stop("'alpha' must be a number in the '(0,0.5]' interval.",call. =FALSE)}
if((alpha<=0|alpha>0.5)){stop("'alpha' must be a number in the '(0,0.5]' interval.",call. =FALSE)}
if(is.numeric(z)==FALSE){stop("'z' must be a positive number.",call. =FALSE)}
if(z<0){stop("'z' must be a positive number.",call. =FALSE)}
if(MinCases>N){
if(M>1){
if(1-pbinom(MinCases-1,N,pst)<alpha){Mr<- MinCases-1;while(1-pbinom(Mr-1,N,pst)<alpha&Mr>1){Mr<- Mr-1};if(Mr==0){Mr<- 1}
stop(c("For these parameters, 'M' must be of at most"," ",Mr),call. =FALSE)
                                    }
       }
stop(c("The minimum number of cases, 'M', must be a positive integer smaller than or equal to"," ",Mr,"."),call. =FALSE)}
if(MinCases<1){stop("The minimum number of cases, 'M', must be a positive integer.",call. =FALSE)}
if(MinCases!=round(MinCases)){stop("The minimum number of cases, 'M', must be a positive integer.",call. =FALSE)}
if(length(M)>1){stop("The minimum number of cases, 'M', must be a single positive integer.",call. =FALSE)}


pst<- 1/(1+z)
if(M>1){
if(1-pbinom(MinCases-1,N,pst)<alpha){Mr<- MinCases-1;while(1-pbinom(Mr-1,N,pst)<alpha&Mr>1){Mr<- Mr-1};if(Mr==0){Mr<- 1}
stop(c("For these parameters, 'M' must be of at most"," ",Mr,"."),call. =FALSE)
                                    }
       }

if(1-pbinom(N-1,N,pst)>alpha){
Nr<- N
while(1-pbinom(Nr-1,Nr,pst)>alpha){Nr<- Nr+1}
stop(c("For this 'N' there is no solution with prob of Type I error smaller than"," ",alpha,". Use 'N' of at least"," ",Nr,"."),call. =FALSE)
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



if(length(Groups)>1){an<- Groups%*%(upper.tri(matrix(0,length(Groups),length(Groups)),diag=T)*1)}else{an<- seq(Groups,N,Groups)
                                                                                                      if(max(an)<N){an<- c(an,N)}
                                                                                                     }

Erro_I<- function(cv,RR=1){


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

}# end function P_erro_I
################################################################

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

} #end function CV.G.Binomial






