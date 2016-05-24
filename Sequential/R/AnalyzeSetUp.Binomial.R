
# -------------------------------------------------------------------------
# Function to perform the unpredictable binomial MaxSPRT surveillance - Version edited at Jan-15-2015
# -------------------------------------------------------------------------

AnalyzeSetUp.Binomial<- function(name, N,alpha=0.05,zp=1,M=1,AlphaSpendType="Wald",rho="n",title="n")
{

pho<- rho
z<- zp

phoref<- rho

if(AlphaSpendType!="Wald"&AlphaSpendType!="power-type"){stop("Set AlphaSpendType= 'Wald' or AlphaSpendType= 'power-type'.",call. =FALSE)}
if(AlphaSpendType=="power-type"&is.numeric(pho)!=TRUE){stop("Symbols and texts are not applicable for 'rho'. It must be a positive number.",call. =FALSE)}
if(pho<0&AlphaSpendType=="power-type"){stop("rho must be greater than zero or equal to the default (rho='n')",call. =FALSE)}

if(pho=="n"){pho<- 0}

if(AlphaSpendType=="Wald"){pho<- 0}

safedir<- getwd()
if(title== "n"){title<- 0}
name1<- name
address<- choose.dir(default = "", caption = paste("Select the folder where the file '",name,"' is going to be saved."))

address1<- tempdir()
y<- substr(address1,1:4,4) ; i<- 2
while(i<=nchar(address1)-3&y!="Temp"&y!="TEMP"){y<- substr(address1,i:(i+3),i+3);i<- i+1}
address1<- substr(address1,1:(i+3),i+3)

if(paste(address)=="NA"){address<- address1}
address2<- data.frame(c(0))
address2[1,1]<- address
setwd(address1)
write.table(address2,paste(name,"address.txt",sep=""),sep=";")
setwd(address)
name<- paste(name,".","txt",sep="")
if(file.exists(name)==TRUE){
stop(c("There already exists a file called"," ",name1,".
","You may want check if some test has been performed for this monitoring before. 
If you really want to overwrite the existent file, please, go to ",address," 
to delete the file '",name,"'. Alternatively, you can delete that file by using the 
following commands: ", "setwd(","'",getwd(),"'",")","; ", "file.remove(","'",name,"'","), 
then try 'AnalyzeSetUp.Binomial' again."),call. =FALSE)
                           }
MinCases<- M

if( sum(is.numeric(alpha))!=1){stop("Symbols and texts are not applicable for 'alpha'. It must be a number in the (0,0.5) interval.",call. =FALSE)}
if( sum(is.numeric(z))!=1){stop("Symbols and texts are not applicable for 'z'. It must be a number greater than zero.",call. =FALSE)}
if( sum(is.numeric(MinCases))!=1){stop("Symbols and texts are not applicable for 'M'. It must be an integer greater than zero.",call. =FALSE)}
if(is.numeric(N)==FALSE){stop("Symbols and texts are not applicable for 'N'. It must be an integer greater than zero.",call. =FALSE)}

if(z<0){stop("'z' must be a number greater than zero.",call. =FALSE)}
if(N<=0){stop("'N' must be an integer greater than zero.",call. =FALSE)}

if(alpha<=0|alpha>0.5||is.numeric(alpha)==FALSE){stop("'alpha' must be a number greater than zero and smaller than 0.5.",call. =FALSE)}
if(MinCases>N||is.numeric(MinCases)==FALSE){stop("'M' must be an integer smaller than or equal to 'N'.",call. =FALSE)}
if(MinCases<1){stop("'M' must be an integer greater than zero.",call. =FALSE)}
if(MinCases!=round(MinCases)){stop("'M' must be an integer.",call. =FALSE)}


alpha1<- alpha
#posi<- length(total_cases)+1
posi<- 2

rejt<- 0
#----- THE MAXSPRT STATISTIC

LLR <- function(cc,n,z){

       if(cc==n){x = n*log(1+z)}else{
         if(z*cc/(n-cc)<=1){x=0}else{
	       x = cc*log(cc/n)+(n-cc)*log((n-cc)/n)-cc*log(1/(z+1))-(n-cc)*log(z/(z+1))
                                    }
                                  } 	
      	x
	                 }
#--------------------------


# -------------------------------------------------------------------------
# Function produces spending alpha for a flat critical value - continuous binomial MaxSPRT
# -------------------------------------------------------------------------

SalphafLAtcv <- function(N,alpha,MinCases,z) {


# alpha = desired alpha level

# MinCases = The minimum number of cases for which a signal is allowed to occur, default=1
# Group = Time between looks in a group sequential trial, must be greater than 0



Perror_I<- function(cv){

absorb = rep(0,N+2)		# Contains the number of events needed at time mu[i] in order to reject H0
aux<- rep(0,N+2)
for(i in 1:N){
	while( LLR(absorb[i],i,z)<cv &absorb[i]<i){ 
		absorb[i]=absorb[i]+1               }
             if(LLR(absorb[i],i,z)>=cv){aux[i]<- 1}
             }
if(MinCases>1){
aux[1:(MinCases-1)]<- 0
              }

absorb[1:N][aux[1:N]==0]<- absorb[absorb[aux==0]]+1

for(i in 1:N){if(absorb[i]<MinCases&i>=MinCases){absorb[i]<- MinCases};if(absorb[i]<MinCases&i<MinCases){absorb[i]<- i+1}}

uc<- absorb-1

ps<- 1/(1+z)


# Auxiliar functions to run the binomial Markov Chain in a fast way:
func_aux2<- function(j,i){ k<- seq(1,uc[i-1]+1); return(sum(p[i-1,k]*dbinom(j-k,1,ps)))} ; func_aux3<- function(i){ k<- seq(1,uc[i-1]+1); return(sum(p[i-1,k]*(1-pbinom(absorb[i]-k,1,ps))))}
func_aux1<- function(i){ j<- matrix(seq(1,absorb[i]),ncol=1) ; return(apply(j,1,func_aux2,i))}

### VERIFYING IF IT IS POSSIBLE TO SPEND SOME AMOUNT OF ALPHA FOR THE OBSERVED TOTAL CASES

if(sum(aux==0)==length(aux)){

         salpha<- rep(0,N)

         error1=0
              
                            }else{ #opens number 1


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
      

# Sums up the probabilities of absorbing states when a signal occurs, to get the alpha level
# ------------------------------------------------------------------------------------------


         salpha<- rep(0,N)
         for(i in 1:N){ salpha[i]<- p[i,absorb[i]+1]}

              error1=0
              for(i in 1:N){ error1=error1+p[i,absorb[i]+1]} 
              
        

                               }## closes number 1
	
return(list(error1,salpha,absorb))

}# end function P_error_I
##------------------------------------------------
####################################################


omega<- matrix(0,nrow=1)
for(i in 1:N){j<- matrix(seq(1,i,1),ncol=1) ; omega<- cbind(omega,matrix(apply(j,1,LLR,i,z),nrow=1))}

omega<- omega[order(omega)]
begin<- sum(omega==0)+1
omega<- omega[begin:length(omega)]

i1<- 1 ;  i2<- length(omega) ; im<- round((i1+i2)/2)

while(i2-i1>1){
 res<- Perror_I(omega[im]) ; error1<- res[[1]]
      if(error1>alpha){i1<- im}else{i2<- im; resold<- res}
 im<- round((i1+i2)/2)
              }

error1<- res[[1]]
if(error1<alpha){salpha<- res[[2]]; absorb<- res[[3]]}else{salpha<- resold[[2]]; absorb<- resold[[3]]}

result<- list(salpha,absorb)
            
names(result)<- c("salpha","absorb")

return(result)

} #end function 

##############################################################################################################

if(pho==0){
result<- SalphafLAtcv(N,alpha=alpha1,MinCases,z)
sa<- result$salpha
absorb1<- result$absorb
if(sum(sa)==0){stop("Choose larger N. It is not possible to find a solution for the desired alpha with the current N choice.",call. =FALSE)}
sum_sa<- sa%*%(upper.tri(matrix(0,length(sa),length(sa)),diag=T))
            }else{
x<- seq(1/N,by=1/N,1)
sum_sa<- alpha*(x^pho)

y<- qbinom(1-alpha,N,1/(1+z)); y<- y+1 # auxiliar variable
if(y<N){absorb1<- c(seq(2,y),y,seq(y+1,N))}else{absorb1<- c(seq(2,y),y)}; absorb1<- c(absorb1,0,0)

          }



## inputSetUp matrix contains:
# line 1: the index for the order of the test (zero entries if we did not have tests before), N, alpha, z, M, posi, base(the line of p where the looping will start in the next test),title, rejt (the index indicating if and when H0 was rejected), pho
# line 2: sum_sa
# line 3: absorb1
# line 4: cases
# line 5: controls
# line 6: current alpha spent
# line 7: absorb for the unpredictable test
# line 8: matching ratio history

inputSetUp<- as.data.frame(matrix(0,8,max(length(absorb1),length(sum_sa),9)))
inputSetUp[1,]<- 0
inputSetUp[1,1:9]<- c(0,N,alpha,z,M,posi,1,0,pho) 
inputSetUp[2,]<- 0
inputSetUp[2,1:length(sum_sa)]<- sum_sa
inputSetUp[3,]<- absorb1 
inputSetUp[6,]<- 0
inputSetUp[7,]<- 0
inputSetUp[8,]<- 0
write.table(inputSetUp,name)
titlecheck<- data.frame(matrix(0,1,1))
if(title!=0){titlecheck[1,1]<- title}
 
message(c("The parameters were successfully set at '",address,"'."),domain = NULL, appendLF = TRUE)
message(c("The temporary directory of your computer has the address of the directory where the settings information of this sequential analysis is saved.
Thus, do not clean the temporary directory before finishing this sequential analysis."),domain = NULL, appendLF = TRUE)

if(AlphaSpendType=="Wald"&phoref!="n"){message(c("The value of 'rho' is ignored, as it is not used when AlphaSpendType='Wald'."),domain = NULL, appendLF = TRUE)}
if(AlphaSpendType=="power-type"&zp!=1){message(c("The value of 'zp' is ignored, as it is not used when AlphaSpendType='power-type'."),domain = NULL, appendLF = TRUE)}

write.table(titlecheck,paste(name1,"title.txt",sep=""))

setwd(safedir)

} ## end function AnalyzeSetUp.Binomial








