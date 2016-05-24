# -------------------------------------------------------------------------
# Function to perform the unpredictable binomial MaxSPRT surveillance - Version edited at Jan-14-2015
# -------------------------------------------------------------------------


Analyze.Binomial<- function(name,test,z,cases,controls,AlphaSpend="n")
{

name1<- name
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
safedir<- getwd()
address1<- tempdir()
y<- substr(address1,1:4,4) ; i<- 2
while(i<=nchar(address1)-3&y!="Temp"&y!="TEMP"){y<- substr(address1,i:(i+3),i+3);i<- i+1}
address1<- substr(address1,1:(i+3),i+3)
setwd(address1)

if(file.exists(paste(name,"address.txt",sep=""))==FALSE){
stop(c("There is no file called"," '",name,"' yet."," You must run the function 'AnalyzeSetUp.Binomial' first."),call. =FALSE)
                                                        }

address<- as.character(read.table(paste(name,"address.txt",sep=""),sep=";")[1,1])
setwd(address)

name<- paste(name,".","txt",sep="")

if(file.exists(name)==FALSE){
stop(c("There is no file called"," '",name,"' yet."," You must run the function 'AnalyzeSetUp.Binomial' first."),call. =FALSE)
                            }

if( sum(is.numeric(cases))!=1|sum(is.numeric(controls))!=1){stop("Symbols and texts are not applicable for 'cases' and 'controls'. They must be integer numbers or zero.",call. =FALSE)}

if( sum(is.numeric(test))!=1|length(test)>1){stop("Symbols and texts are not applicable for 'test'. It must be an integer greater than zero.",call. =FALSE)}

if(test<0){stop("'test' must be an integer greater than zero.",call. =FALSE)}

if(sum(cases<0)>0|sum(controls<0)>0){stop("The counts 'cases' and 'controls' must be integers greater than or equal to zero.",call. =FALSE)}

if(sum(is.numeric(z))!=1){stop("Symbols and texts are not applicable for 'z'. It must be a number greater than zero.",call. =FALSE)}

if(sum(z<=0)>0){stop("The entries of 'z' must be numbers greater than zero.",call. =FALSE)}

if( sum(is.numeric(cases))==1&sum(is.numeric(controls))==1){if(sum(controls+cases)==0){stop("At least one of the entries of 'cases' + 'controls' must be greater than zero.",call. =FALSE)}}

if(length(controls)!=length(cases)|length(controls)!=length(z)){stop("'cases', 'controls' and 'z' are vectors that must have the same dimension.",call. =FALSE)}


inputSetUp<- read.table(name)
if(inputSetUp[1,1]!=test-1){stop(c("The current test should be"," ",inputSetUp[1,1]+1,". ",
"If you do not have information about previous tests, see the user manual for more details."),call. =FALSE)}



### Organizing the matching ratios subject by subject

cases_new<- cases
controls_new<- controls

cases<- sum(cases)        
controls<- sum(controls)  

z_new<- rep(z,cases_new+controls_new) 

###

## inputSetUp is a matrix containing:
# line 1: zero indicates that we have not performed tests yet. Here we have N, alpha, z, M, posi, base(the line to start running p),title, rejt
# line 2: sum_sa
# line 3: absorb1
# line 4: cases
# line 5: controls
# line 6: current alpha spent
# line 7: absorb for the unpredictable test
# line 8: matching ratio history
# lines 9 to (N+2): p matrix

# zr is the vector of matching ratios for the current test 

N<- inputSetUp[1,2] ; alpha<- inputSetUp[1,3] ; zr<- inputSetUp[1,4] ; z<- inputSetUp[1,4]; M<- inputSetUp[1,5] ; posi<- inputSetUp[1,6] ; MinCases<- M ; base<- inputSetUp[1,7] 
rejt<- inputSetUp[1,8] ; pho<- inputSetUp[1,9]
sum_sa<- inputSetUp[2,1:ncol(inputSetUp)] ; sum_saref<- as.numeric(apply(sum_sa,1,max))  ; sum_sa_used<- inputSetUp[6,1:test]
absorb1<- inputSetUp[3,1:(N+2)]



titlecheck<- paste(name1,"title.txt",sep="")
title<- read.table(titlecheck)
title<- title[1,1]
if(title==0){title<- " "}else{title<- as.character(title)}

## Updating information for future tests.
inputSetUp[4,test]<- cases
inputSetUp[5,test]<- controls

########## Updating the number of rows of matrix p for new cases


if(cases+controls>0){inpaux<- as.data.frame(matrix(0,cases+controls,ncol(inputSetUp)));names(inpaux)<- names(inputSetUp)
inputSetUp<- rbind(inputSetUp,inpaux)
                    }

m<- max(N, sum(inputSetUp[4,1:test])+sum(inputSetUp[5,1:test]))

                      

######### Adjusting the matrix p for the case where the number of cases is greater than N

if(sum(inputSetUp[4,1:test])+sum(inputSetUp[5,1:test])>N){
                                                          dif1<- sum(inputSetUp[4,1:test])+sum(inputSetUp[5,1:test])-N 
                                                          inputSetUp<- cbind(inputSetUp,matrix(0,nrow(inputSetUp),dif1))
                                                          inpaux<- as.data.frame(matrix(0,dif1,ncol(inputSetUp)))
                                                          names(inpaux)<- names(inputSetUp)
                                                          inputSetUp<- rbind(inputSetUp,inpaux)
                                                          #aux<- c(aux,rep(0,dif1))
                                                          m<-  sum(inputSetUp[4,1:test])+sum(inputSetUp[5,1:test])
                                                         }

## Updating information for future tests.
inputSetUp[1,1]<- inputSetUp[1,1]+1

# Updating historical matching ratio case by case (zhc) and test by test(zh) 

zh<- rep(0,test)

if(test==1){inputSetUp[8,1:sum(inputSetUp[4,1]+inputSetUp[5,1])]<- z_new;  zh[1]<- mean(z_new)}else{
     inputSetUp[8,(sum(inputSetUp[4,1:(test-1)]+inputSetUp[5,1:(test-1)])+1):sum(inputSetUp[4,1:test]+inputSetUp[5,1:test])]<- z_new

     for(i in 1:test){zh[i]<-  mean(as.numeric(inputSetUp[8,1:sum(inputSetUp[4,1:i]+inputSetUp[5,1:i])]))}
                                                                                                   }




zhc<- as.numeric(inputSetUp[8,1:sum(inputSetUp[4,1:test]+inputSetUp[5,1:test])])    # historical matching ratio case by case
p<- inputSetUp[9:nrow(inputSetUp),]
absorb<- as.numeric(inputSetUp[7,1:m])


cases<- inputSetUp[4,1:test]
controls<- inputSetUp[5,1:test]

#### Adjusting the target alpha spending

salphaux1<- 0  

if(AlphaSpend=="n"){AlphaSpend<- c(0)}

if( (AlphaSpend!=0)& sum(inputSetUp[4,1:test])+sum(inputSetUp[5,1:test])>0){ # starts *  
if(sum(inputSetUp[4,1:test])+sum(inputSetUp[5,1:test])>N){AlphaSpend<- inputSetUp[1,3]; salphaux1<- 1} 

if(m<=N){

if( sum(is.numeric(AlphaSpend))!=1){stop("Symbols and texts are not applicable for 'AlphaSpend'. If you want to use the default, use 'n'. Otherwise,  'AlphaSpend' must be a positive number smaller than or equal to 'alpha'.",call. =FALSE)}
if( sum(length(AlphaSpend))!=1){stop("'AlphaSpend' must be a single value, not a vector.",call. =FALSE)}
if(AlphaSpend<0){stop("'AlphaSpend' must be a positive number smaller than 'alpha'.",call.=FALSE)}
if(AlphaSpend>alpha){stop(c("'AlphaSpend' must be smaller than or equal to ",alpha,"."),call. =FALSE)}

auxsas<- sum(inputSetUp[4,1:test])+sum(inputSetUp[5,1:test])
auxsas1<- sum(inputSetUp[4,1:(test-1)])+sum(inputSetUp[5,1:(test-1)])

if(test>1&auxsas1>0){
if(AlphaSpend<sum_sa[auxsas1]){stop(c("'AlphaSpend' must be greater than or equal to ",round(sum_sa[auxsas1],6)," because it has already been used a target of ",round(sum_sa[auxsas1],6)," in the previous test."),call. =FALSE)}
                    }

sum_sa[auxsas]<- AlphaSpend

if(auxsas<N){for(ii in (auxsas+1):N){if(sum_sa[ii]<sum_sa[ii-1]){sum_sa[ii]<- sum_sa[ii-1]}}}

         }else{AlphaSpend<- c(0)}

                                                                        } # ends *
#### The adjusment for the target alpha spending terminates here.
        

total_cases<- as.numeric(cases+controls)%*%(upper.tri(matrix(0,length(cases),length(cases)),diag=T)*1)
exposed_cases<- as.numeric(cases)%*%(upper.tri(matrix(0,length(cases),length(cases)),diag=T)*1)

### Initial line to start the running inside matrix p

#if(length(total_cases)>1){base<- 1+total_cases[length(total_cases)-1]}else{base<- 1}

alpha1<- alpha
posi<- length(total_cases)+1
CIgamma<- 1-2*alpha


#-------------------------------------------------------------------------
## HERE WE CHECK IF THERE IS AT LEAST ONE ADVERSE EVENT TO PERFORM A TEST
#-------------------------------------------------------------------------

if(sum(total_cases)==0|cases[test]+controls[test]==0){
                          
if(cases[test]+controls[test]==0&sum(total_cases)>0){

result<- data.frame(matrix(0,length(total_cases)+1,12))
result[2:(length(total_cases)+1),1]<- seq(1,length(total_cases),1)
colnames(result)<- c(" "," "," ", "-----","Cumulative","--------"," "," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","Cases","Controls","Cases","Controls","E[Cases]","RR","LLR","target","actual","CV","H0 rejected")
result[2:nrow(result),2]<- as.numeric(cases) ; result[2:nrow(result),3]<- as.numeric(controls) ; result[2:nrow(result),4]<- t(as.numeric(cases)%*%(upper.tri(matrix(0,length(cases),length(cases)),diag=T)*1)) 




## Useful auxiliar variable
Exp<- rep(0,test)
for(i in 1:test){Exp[i]<- sum(1/(1+zhc[1:total_cases[i]]))  }  

result[2:(length(total_cases)+1),5]<- t(as.numeric(controls)%*%(upper.tri(matrix(0,length(controls),length(controls)),diag=T)*1)) ; result[2:nrow(result),6]<- paste(round(Exp,2))

result[nrow(result),9]<- paste(round(as.numeric(sum_sa[total_cases[(test-1)]]),4))
for(i in 2:nrow(result)){if(total_cases[i-1]>N){result[i,9]<- alpha1}else{if(total_cases[i-1]==0){result[i,9]<- 0}else{result[i,9]<- paste(round(as.numeric(sum_sa[total_cases[i-1]]),4))}}}
if(rejt>0){result[(2+rejt):nrow(result),9]<- paste("NA")}

for(i in 2:(nrow(result)-1)){if(total_cases[i-1]==0){result[i,11]<- paste("NA")}else{if(absorb[total_cases[i-1]]>total_cases[i-1]){result[i,11]<- paste("NA")}else{result[i,11]<- absorb[total_cases[i-1]]}}}
result[nrow(result),11]<- result[nrow(result)-1,11]
 

result[2:nrow(result),12]<- paste("No") ; if(rejt>0){result[(2+rejt-1):nrow(result),12]<- paste("Yes")}
result[2:nrow(result),10]<- paste(round(as.numeric(sum_sa_used),4)) ; result[nrow(result),10]<- paste(round(as.numeric(sum_sa_used[test-1]),4))
# CIgamma<- round(1-2*as.numeric(sum_sa_used[test-1]),4)
if(rejt>0){result[(2+rejt):nrow(result),10]<- paste("NA")} ; if(rejt>0){result[(2+rejt):nrow(result),11]<- paste("NA")}
result[2:nrow(result),7]<- round(as.numeric(zh*exposed_cases/(total_cases-exposed_cases)),2) 

for(i in 1:(nrow(result)-1)){if(total_cases[i]==0){result[(1+i),8]<- paste("NA") ; result[(1+i),7]<- paste("NA"); result[(1+i),9]<- paste("NA")}else{result[(1+i),8]<- round(LLR(exposed_cases[i],total_cases[i],zr),3)}}
 

if(rejt>0){ # 3
message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)      
          message(paste(c("=>    H0 was rejected on test"," ",rejt,". ","No further sequential analyses are needed.")),domain = NULL, appendLF = TRUE)  
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
          }else{
if(max(total_cases)<N){
message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
message("There is not sufficient information to perform a test at this moment.",domain = NULL, appendLF = TRUE)                                                                                                                                                                                                           
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
                       }else{
message(  " ",domain = NULL, appendLF = TRUE)
message(c("                                ",levels(title)[1]),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------", domain = NULL, appendLF = TRUE)
message(c("The upper limit on the length of surveillance has been reached."), domain = NULL, appendLF = TRUE)
if(salphaux1==1){
message("Then, the 'AlphaSpend' input is no longer used, and by default the target is alpha.")
                }                                                         
message(c("You may now end the sequential analysis without rejecting H0."), domain = NULL, appendLF = TRUE)
if(alpha1-sum_sa_used[length(total_cases)]>0){
message(c("There is still ",round(alpha1-max(sum_sa_used),6)," alpha to spend if you wish continue with more analyses."), domain = NULL, appendLF = TRUE)}else{
message(c("Additionally, the overall amount of alpha has already been spent."), domain = NULL, appendLF = TRUE)
                           }
message(" ", domain = NULL, appendLF = TRUE)  
message("-------------------------------------------------------------------------------------------", domain = NULL, appendLF = TRUE)                                                         
                                                                                                                                                              }
         } # 3 close
               
print(result,right=TRUE,row.names=FALSE)

message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: N= ",N,", alpha= ",alpha,", zp= ",z," and M= ",M,"."),domain = NULL, appendLF = TRUE)
#message(c(R.version.string,". ","Analysis performed at ",date(),"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)

if(test==1){inputSetUp[6,test]<- 0}else{inputSetUp[6,test]<- inputSetUp[6,test-1]; inputSetUp[2,test]<- inputSetUp[2,test-1]} 


                         
                                                   }else{ #1 
result<- data.frame(matrix(0,length(total_cases)+1,12))
result[2:(length(total_cases)+1),1]<- seq(1,length(total_cases),1)
colnames(result)<- c(" "," "," ", "-----","Cumulative","--------"," "," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","Cases","Controls","Cases","Controls","E[Cases]","RR","LLR","target","actual","CV","H0 rejected")
result[2:nrow(result),2]<- as.numeric(cases) ; result[2:nrow(result),3]<- as.numeric(controls) ; result[2:nrow(result),4]<- t(as.numeric(cases)%*%(upper.tri(matrix(0,length(cases),length(cases)),diag=T)*1)) 

## Useful auxiliar variable
Exp<- rep(0,test) 
for(i in 1:test){Exp[i]<- sum(1/(1+zhc[1:total_cases[i]]))  }  

result[2:(length(total_cases)+1),5]<- t(as.numeric(controls)%*%(upper.tri(matrix(0,length(controls),length(controls)),diag=T)*1)) ; result[2:nrow(result),6]<- paste(round(Exp,2))
result[2:nrow(result),7]<- paste("NA") ; result[2:nrow(result),8]<- paste("NA"); result[2:nrow(result),9]<- paste("NA") ; result[2:nrow(result),11]<- paste("NA") ; result[2:nrow(result),12]<- paste("No")


message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
message("There is no sufficient information to perform a test at this moment.",domain = NULL, appendLF = TRUE)                                                                                                                                                                                                           
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
               
print(result,right=TRUE,row.names=FALSE)

message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: N= ",N,", alpha= ",alpha,", zp= ",z," and M= ",M,"."),domain = NULL, appendLF = TRUE)
#message(c(R.version.string,". ","Analysis performed at ",date(),"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)

## Updating information for future tests.
inputSetUp[2,test]<- 0
inputSetUp[6,test]<- 0
                                                  }# close 1

                                        }else{ # 10

totalg<- total_cases
exposedg<- exposed_cases
if(length(exposed_cases)>1){exposed_cases<- exposed_cases[(sum(total_cases==0)+1):length(total_cases)]}
if(length(AlphaSpend)>1 & sum(total_cases==0)>0){alphau<- sum(AlphaSpend[1:(sum(total_cases==0)+1)]);AlphaSpend<- AlphaSpend[(sum(total_cases==0)+1):length(total_cases)];AlphaSpend[1]<- alphau}
total_cases<- total_cases[(sum(total_cases==0)+1):length(total_cases)]

##
# Verify if H0 has already been rejected in previous tests
##


      if(rejt>0){ # 11

result<- data.frame(matrix(0,length(totalg)+1,12))
result[2:nrow(result),1]<- seq(1,length(totalg),1)
colnames(result)<- c(" "," "," ", "-----","Cumulative","--------"," "," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","Cases","Controls","Cases","Controls","E[Cases]","RR","LLR","target","actual","CV","H0 rejected")
result[2:nrow(result),2]<- as.numeric(cases) ; result[2:nrow(result),3]<- as.numeric(controls) ; result[2:nrow(result),4]<- t(as.numeric(cases)%*%(upper.tri(matrix(0,length(cases),length(cases)),diag=T)*1)) 

## Useful auxiliar variable
Exp<- rep(0,test)
for(i in 1:test){Exp[i]<- sum(1/(1+zhc[1:total_cases[i]]))  }  

result[2:(length(totalg)+1),5]<- t(as.numeric(controls)%*%(upper.tri(matrix(0,length(controls),length(controls)),diag=T)*1)) ; result[2:nrow(result),6]<- paste(round(Exp,2))

for(i in 2:nrow(result)){
if(i>=rejt+2){result[i,9]<- paste("NA")}else{
if(totalg[i-1]>N){result[i,9]<- alpha1}else{
if(totalg[i-1]==0){result[i,9]<- paste("NA")}else{result[i,9]<- paste(round(as.numeric(sum_sa[totalg[i-1]]),4))}
                                                }
                                            }
                        }

for(i in 2:(nrow(result)-1)){
if(totalg[i-1]==0|i>=rejt+2){result[i,11]<- paste("NA")}else{if(absorb[totalg[i-1]]>totalg[i-1]){result[i,11]<- paste("NA")}else{result[i,11]<- absorb[totalg[i-1]]}}
                            }
result[nrow(result),11]<- paste("NA")

result[2:nrow(result),12]<- paste("No") ; result[(2+rejt-1):nrow(result),12]<- paste("Yes")

result[2:nrow(result),10]<- paste(round(as.numeric(sum_sa_used),4)) ; result[nrow(result),10]<- paste(round(as.numeric(sum_sa_used[test-1]),4))
# CIgamma<- round(1-2*as.numeric(sum_sa_used[test-1]),4)
result[(2+rejt):nrow(result),10]<- paste("NA") 

result[(2+sum(totalg==0)):nrow(result),7]<- round(as.numeric(zh[(1+sum(totalg==0)):(nrow(result)-1)]*exposed_cases/(total_cases-exposed_cases)),2)
if(sum(totalg==0)>0){result[2:(1+sum(totalg==0)),7]<- paste("NA")}

for(i in (2+sum(totalg==0)):(nrow(result))){result[i,8]<- round(LLR(exposed_cases[i-sum(totalg==0)-1],total_cases[i-sum(totalg==0)-1],zh[i-1]),3)}  
if(sum(totalg==0)>0){result[2:(1+sum(totalg==0)),8]<- paste("NA")}
                                              
                                              
                                              decision<- paste("H0 was rejected at test"," ",rejt)
                                              message(  " ",domain = NULL, appendLF = TRUE)
                                              message(c("                                ",title),domain = NULL, appendLF = TRUE)
                                              message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
                                              message(paste(c("=>    H0 was rejected on test"," ",rejt,". ","No further sequential analyses are needed.")),domain = NULL, appendLF = TRUE)
                                              message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
print(result,right=TRUE,row.names=FALSE)

message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: N= ",N,", alpha= ",alpha,", zp= ",z," and M= ",M,"."),domain = NULL, appendLF = TRUE)
#message(c(R.version.string,". ","Analysis performed at ",date(),"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)

              }else{ # related to 11

###### HERE WE HAVE AN AUXILIARY CODE TO OBTAIN P-VALUES AND CONFIDENCE INTERVALS LATER 
##### 
cumulative_prob<- function(total_cases,absorb,RRseq,zhc)
{

ps<- matrix(0,length(zhc),length(RRseq))
for(j in 1:nrow(ps)){ps[j,]<- RRseq/(RRseq+zhc[j])}
    

N<- max(total_cases)

uc<- absorb-1


p<- array(rep(0,N*(N+2)*length(RRseq)),dim=c(N,N+2,length(RRseq)))	# p[i,j] is the probability of having j-1 cases at time mu[i]

                  									# starting probabilities are all set to zero

# Calculating the first row in the p[][] matrix for which there is a chance to reject H0
# --------------------------------------------------------------------------------------

for(s in 1:absorb[1]){ p[1,s,]=dbinom(s-1,1,ps[1,])}		# Probability of having s-1 cases at time mu[1]
p[1,absorb[1]+1,]=1-pbinom(absorb[1]-1,1,ps[1,])			# Probability of rejecting H0 at time mu[1]
p2<- p

if(N>1){
for(i in 2:N){

	for(j in 1:(absorb[i])){ # This loop calculates the p[][] matrix, one column at a time, from left to right                 
			
		for(k in 1:(uc[i-1]+1)){
			p[i,j,]=p[i,j,]+p[i-1,k,]*dbinom(j-k,1,ps[i,])	# Calculates the standard p[][] cell values
                  p2[i,j,]=p[i,j,]
                                   }
                             }

 for(k in 1:(uc[i-1]+1)){
            ppaux<- 1-pbinom(absorb[i]-k,1,ps[i,]) ; ppaux2<- pbinom(absorb[i]-k,1,ps[i,]) 
		p[i,absorb[i]+1,]=p[i,absorb[i]+1,]+p[i-1,k,]*ppaux
            p2[i,absorb[i]+1,]=p2[i,absorb[i]+1,]+p2[i-1,k,]*ppaux2
                        }
                               
                 
} # end for i	
      }

prob1=0; for(i in 1:N){ prob1=prob1+p[i,absorb[i]+1,]} ; prob2<- p2[N,absorb[N]+1,]

out<- list(prob1,prob2)  
return(out)

}

###### HERE WE CLOSE THE AUXILIARY CODE TO OBTAIN P-VALUES AND CONFIDENCE INTERVALS LATER 
#####

# -------------------------------------------------------------------------
# Function produces critical value for the group binomial MaxSPRT
# -------------------------------------------------------------------------

CVGbinPar<- function(arriv,alphasa,MinCases=1,zhc,absorb,aux,t1,pref,base){

alpha<- alphasa
# T = maximum length of surveillance, defined in terms of expected counts under H0
# alpha = desired alpha level

# MinCases = The minimum number of cases for which a signal is allowed to occur, default=1

Nn<- max(arriv)  # 

ps<- 1/(1+zhc)

Perror_I<- function(nn,absor=0){

if(Nn>1){

absorb[Nn]=nn
aux[Nn]<- 1
uc<- absorb-1


p<- pref

# Calculating the first row in the p[][] matrix for which there is a chance to reject H0
# --------------------------------------------------------------------------------------

if(base==1){
for(s in 1:absorb[1]){ p[1,s]=dbinom(s-1,1,ps[1])}		# Probability of having s-1 cases at time mu[1]
p[1,absorb[1]+1]=1-pbinom(absorb[1]-1,1,ps[1])			# Probability of rejecting H0 at time mu[1]
base<- 2
           }
base2<- base
if(cont==2){base2<- Nn}
for(i in base2:Nn){

	for(j in 1:(absorb[i])){ # This loop calculates the p[][] matrix, one column at a time, from left to right                 
			
		for(k in 1:(uc[i-1]+1)){
			p[i,j]=p[i,j]+p[i-1,k]*dbinom(j-k,1,ps[i])	# Calculates the standard p[][] cell values
                                   }
                             }

 for(k in 1:(uc[i-1]+1)) 
		p[i,absorb[i]+1]=p[i,absorb[i]+1]+p[i-1,k]*(1-pbinom(absorb[i]-k,1,ps[i]))# Calculates the diagonal absorbing states where H0 is rejected

                                                
} # end for i
       	
         
# Sums up the probabilities of absorbing states when a signal occurs, to get the alpha level
# ------------------------------------------------------------------------------------------

alpha=0
for(i in 1:Nn){ alpha=alpha+p[i,absorb[i]+1]}

}else{
       absorb[1] = 1		# Contains the number of events needed at time mu[i] in order to reject H0
      
	 alpha<- 1-pbinom(0,1,ps[1])  
                           

     }## end if(Nn>1) 				
	
if(absor==0){return(list(alpha,p))}else{return(absorb[Nn])}

}# end function P_error_I
################################################################

error1<- 1
if(length(total_cases)>1){N1<- absorb[total_cases[length(total_cases)-1]]}else{N1<- 1}  ; N2<- Nn

if(N2-N1<=1){             
             result<- Perror_I(N2);error1<- result[[1]]
             if(error1>alpha){error1<- sum_sa_used[length(sum_sa_used)-1];absorb[Nn]<- Nn+1}else{
                                                                                        if(error1==alpha){absorb[Nn]<- N2}else{
                                                                                        if(N2==N1){absorb[Nn]<- N2}
                                                                                        if(N2>N1){                                                                                                 
                                                                                                 result<- Perror_I(N1);error2<- result[[1]]
                                                                                                 if(error2<=alpha){absorb[Nn]<- N1 ; error1<- error2}
                                                                                                 }
                                                                                                                              }
                                                                                              }     
            }else{
precision<- 0.0000001
lim<- log((N2-N1)/precision)/log(2)
cont<- 1
Nm<- ceiling((N2+N1)/2)	
auxNN<- 0 
while(cont<=lim&abs(alpha-error1)>precision&auxNN==0){
res<- Perror_I(Nm) 
error1<- res[[1]]
if(cont==1){cont<- cont+1; pref<- res[[2]];pref[Nn,]<- 0}
if(error1<alpha){N2<- Nm;Nm<- floor((N2+N1)/2)}else{N1<- Nm;Nm<- ceiling((N2+N1)/2)}

if(N2-N1<=1){      
             auxNN<- 1       
             result<- Perror_I(N2);error1<- result[[1]]
             if(error1>alpha){error1<- sum_sa_used[length(sum_sa_used)-1];absorb[Nn]<- Nn+1}else{
                                                                                        if(error1==alpha){absorb[Nn]<- N2}else{
                                                                                        if(N2==N1){absorb[Nn]<- N2}
                                                                                        if(N2>N1){                                                                                                 
                                                                                                 result<- Perror_I(N1);error2<- result[[1]]
                                                                                                 if(error2<=alpha){absorb[Nn]<- N1 ; error1<- error2}else{absorb[Nn]<- N2}
                                                                                                 }
                                                                                                                              }
                                                                                              }                                           
            }

                                                    }
                 }

result<- Perror_I(absorb[Nn])
p<- result[[2]]
if(absorb[Nn]>Nn){aux[Nn]<- 0}else{aux[Nn]<- 1} 
return(list(absorb[Nn],aux[Nn],error1,p))

} #end function CVGbinPar


### OBTAINING THE PROBABILITIES FOR THE OBSERVED CASES

## Here we use the alpha spending according to that specified by the user. 

                   
if(sum(AlphaSpend)>0&sum_saref==0){ # 4

                                            
                         t1<- 0
                           ps0<- 1/(1+zhc)
                            saini<- 1
                              while(saini>AlphaSpend[1]){
                                                          t1<- t1+1 
                                                          s<- 1                       
                                  while(saini>AlphaSpend[1]&s<=t1){saini<- 1-pbinom(s-1,t1,ps0[t1]);s<- s+1}                           
                                                        }

                                                         if(t1<length(absorb)){absorb[1:t1]<- seq(1,t1,1)+1;absorb[t1]<- 0}
                      
                                  }else{
t1<- 1
if(max(absorb)==0){
while( as.numeric(absorb1[t1])>t1){absorb[t1]<-  as.numeric(absorb1[t1]);t1<- t1+1}
if(t1<min(total_cases)){absorb[t1:total_cases[1]]<- seq(t1+1,total_cases[1]+1,1)}
                  }else{while( as.numeric(absorb1[t1])>t1){t1<- t1+1}}

                                       } # end 4

aux<- rep(0,m)

i<- length(total_cases)

if(rejt==0){ # 2

if(t1<=max(total_cases)){ #1 verifying if the Type I error is greater than zero for the observed sample

if(i>1){if(total_cases[i]-total_cases[i-1]>1){absorb[(total_cases[i-1]+1):(total_cases[i]-1)]<- seq(total_cases[i-1]+1,total_cases[i]-1,1)+1}}

if(total_cases[i]<N|sum(AlphaSpend)>0){result<- CVGbinPar(arriv<- total_cases[1:i],alphasa=sum_sa[total_cases[i]],MinCases,zhc,absorb,aux,t1,pref=p,base)}else{

                 result<- CVGbinPar(arriv= total_cases[1:i],alphasa=alpha1,MinCases,zhc,absorb,aux,t1,pref=p,base)
                                                                                                                                                          }

absorb[total_cases[i]]<- result[[1]]
sum_sa_used[test]<- result[[3]]
p<- result[[4]]


                     } # end #1
           }else{
                 absorb[total_cases[i]]<- total_cases[i]+1
                 sum_sa_used[test]<-sum_sa_used[test-1] 
                } # end 2

if(m>N){total_cases_aux<- c(total_cases[total_cases<=N],rep(N,sum(total_cases>N)))}else{total_cases_aux<- c(total_cases[total_cases<=N])}
if(sum(AlphaSpend)>0|t1>max(total_cases)){total_cases_aux<- total_cases}

## Checking if H0 is going to be rejected

###
#  Here we calculate the p-value and the confidence interval if the surveillance can be interrupted
###

# if(as.numeric(sum_sa_used[test])>0){CIgamma<- 1-2*as.numeric(sum_sa_used[test])}

decision1<- 0
if(length(exposed_cases)==length(total_cases)){if(sum(exposed_cases>=absorb[total_cases])>0){decision1<- 1}}else{decision1<- 1*(max(exposed_cases)>=max(absorb[total_cases]))}
if(max(total_cases)>=N){decision1<- 1}
if(rejt>0){decision1<-0}

if(decision1==1){#1
                                    if(absorb[length(absorb)]==0){absorbp<- absorb[-length(absorb)]}else{absorbp<- absorb}
                                    absorbp[length(absorbp)]<-  exposed_cases[length(exposed_cases)]
                                    for(i in max(total_cases):2){if(absorbp[i]<absorbp[i-1]){absorbp[i-1]<- absorbp[i]}}  
                                    
          if(length(exposed_cases)>1){
                                      for(i in 1:(length(exposed_cases)-1)){if(exposed_cases[i]>=absorbp[total_cases[i]]){absorbp[total_cases[i]]<- exposed_cases[i]+1}}
                                      for(i in 2:length(absorbp)){if(absorbp[i]>=absorbp[i-1]+2){absorbp[i]<- absorbp[i-1]+1}}
                                     }
                                   

                                    
                                    if(sum(exposed_cases)>0){#2

                                                             #p_value<- cumulative_prob(total_cases,absorbp,RRseq=1,zhc)
                                                             #p_value<- p_value[[1]]
                             
RR1<- 0.1
RRseq<- seq(0.01,100,0.01) ; j<- max(total_cases) ; cv<- exposed_cases[length(exposed_cases)] 
RRref<- RRseq[max(seq(1,length(RRseq),1)*(pbinom(cv-1,j,RRseq/(RRseq+max(as.numeric(zhc))))>(1-CIgamma)/2))]
RR2<- RRref

## Obtaining RR_L and RR_U
if(RR2>=1.2){RRseq<- c(seq(RR1,0.9,0.1),seq(0.91,1.1,0.01),seq(1.2,RR2,0.1))}else{if(RR2>0.9){RRseq<- c(seq(RR1,0.9,0.1),seq(0.91,RR2,0.01))}else{RRseq<- seq(RR1,RR2,0.1)}}
probs<- cumulative_prob(total_cases,absorbp,RR=RRseq,zhc)
res1<- probs[[1]]
res2<- probs[[2]]

# Obtaining RR_L
if(min(res1)<=(1-CIgamma)/2){RR_L<- max(RRseq[seq(1,length(RRseq),1)*(res1<=(1-CIgamma)/2)])}else{
                RR11<- seq(0.001,RR1-0.001,0.001)
                probs<- cumulative_prob(total_cases,absorbp,RR=RR11,zhc) ; res1<- probs[[1]]
                RR_L<- max(RR11[seq(1,length(RR11),1)*(res1<=(1-CIgamma)/2)])                                                                              
                                                                                             }

# Obtaining RR_U
if(min(res2)<=(1-CIgamma)/2){RR_U<- max(RRseq[1+seq(1,length(RRseq),1)*(res2>(1-CIgamma)/2)])}else{
                while(min(probs[[2]])>(1-CIgamma)/2){RR2<- RR2+0.1; probs<- cumulative_prob(total_cases,absorbp,RR=RR2,zhc)}
                                                                                               RR_U<- RR2
                                                                                             }

      
                                                           }else{#p_value<- 1 
                                                                RR_L<- (1-(((1+CIgamma)/2)^(1/j)))/(2-(((1+CIgamma)/2)^(1/j)))
                                                                RR_U<- (1-(((1-CIgamma)/2)^(1/j)))/(2-(((1-CIgamma)/2)^(1/j)))
                                                               }# close 2                                         
    
                   }# close 1

###_____________________________________________________________________
### Building information to show the results by using graphics and tables
###_____________________________________________________________________

auxii<- 0
i<- 1
while(i <=length(totalg)&auxii==0){if(totalg[i]>0){auxii<- 1*(exposedg[i]>=absorb[totalg[i]])};if(auxii==1){rejt<- i};i<- i+1}


if(length(exposed_cases)==length(total_cases)){if(sum(exposed_cases>=absorb[total_cases])>0){stop1<- 1}else{stop1<- 0}}else{stop1<- 0}

if(max(exposed_cases)>=max(absorb[total_cases])|stop1==1){ # 1
                                              if(sum(exposed_cases>=absorb[total_cases])==1&max(exposed_cases)>=max(absorb[total_cases])){ # 2
                                              decision<- c("Reject H0"); if(RR_L<1){RR_L<- 1}
                                              message(  " ",domain = NULL, appendLF = TRUE)
                                              message(c("                                ",title),domain = NULL, appendLF = TRUE)
                                              message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
                                              message("=>    Reject H0. No further sequential analyses are needed.",domain = NULL, appendLF = TRUE)
                                              message(c("      ",paste(round(CIgamma*100,2)),"% ","Confidence interval for the actual relative risk=(",round(RR_L,2),",",round(RR_U,2),"]"), domain = NULL, appendLF = TRUE)
                                              message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
#message(c("P-value=",round(p_value,4)), domain = NULL, appendLF = TRUE)

                                                                                                                                          } # close 2

                                                         }else{ # related to 1

decision<- c("Do not reject H0")

if(max(total_cases)>=N&sum_sa_used[length(total_cases)]<alpha1){ # 4
message(  " ",domain = NULL, appendLF = TRUE)
message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------", domain = NULL, appendLF = TRUE)
message(c("The upper limit on the length of surveillance has been reached."), domain = NULL, appendLF = TRUE) 
if(salphaux1==1){
message("Then, the 'AlphaSpend' input is no longer used, and by default the target is alpha.")
                }                                                        
message(c("You may now end the sequential analysis without rejecting H0."), domain = NULL, appendLF = TRUE)
message(c("There is still ",round(alpha1-max(sum_sa_used),6)," alpha to spend if you wish continue with more analyses."), domain = NULL, appendLF = TRUE)
message(" ", domain = NULL, appendLF = TRUE)  
message(c("  ",paste(round(CIgamma*100)),"% ","Confidence interval for the actual relative risk=[",round(RR_L,2),",",round(RR_U,2),"]"), domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------", domain = NULL, appendLF = TRUE)                                                         
#message(c("P-value=",round(p_value,4)), domain = NULL, appendLF = TRUE)

                                                               } # close 4

if(max(total_cases)>=N&sum_sa_used[length(total_cases)]==alpha1){message("The surveillance must be stopped according to the desired alpha.",domain = NULL, appendLF = TRUE)}
if(max(total_cases)<N){ # 5
message(  " ",domain = NULL, appendLF = TRUE)
message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
message("=>   Do not reject H0. Proceed to a new test as soon as you have more data.", domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
                      } # close 5
                                                              }# close 1
  
######  Building the tables

names1<- rep(0,length(total_cases))
for(i in 1:length(total_cases)){names1[i]<- paste("Test",i)}
result1<- as.data.frame(matrix(0,3,length(total_cases)))
names(result1)<- names1
result1[1,]<- as.character(total_cases)
result1[2,]<- round(sum_sa[total_cases_aux],4)
result1[3,]<- round(sum_sa_used[(sum(totalg==0)+1):length(sum_sa_used)],4)
if(posi<=length(totalg)){result1[2:3,(posi-sum(totalg==0)):length(total_cases)]<- paste("NA")}
rownames(result1)<- c("Total Cases","Target alpha spending","Actual alpha spent")

if(length(exposed_cases)==length(total_cases)){
result2<- as.data.frame(result1)
rownames(result2)<- c("Total Cases","Exposed cases","Needed to reject H0")
result2[3,]<- absorb[total_cases]
if(posi<=length(totalg)){result2[3,(posi-sum(totalg==0)):length(total_cases)]<- paste("NA")}
if(sum(absorb[total_cases]>total_cases)>0)(result2[3,absorb[total_cases]>total_cases]<- paste(rep("NA",sum(absorb[total_cases]>total_cases))))
result2[2,]<- exposed_cases
                                               }else{result2<- as.data.frame(matrix(0,3,1))
                                                     rownames(result2)<- c("Cases","Exposed cases","Needed to reject H0")
                                                     result2[1,1]<- as.character(max(total_cases))
                                                     result2[3,1]<- max(absorb[total_cases])                                                     
                                                     result2[2,1]<- max(exposed_cases)
                                                     names(result2)<- paste("Test",length(total_cases))
                                                    }

if(length(total_cases)==length(exposed_cases)){# 11

result3<- as.data.frame(matrix(0,4,length(total_cases)))
names(result3)<- names1
rownames(result3)<- c("Expected Cases","Estimated RR","LLR","CV")

## Useful auxiliar variable
Exp<- rep(0,test)
for(i in 1:test){Exp[i]<- sum(1/(1+zhc[1:total_cases[i]]))  }  

for(i in 1:length(total_cases)){
                               result3[1,i]<- paste(round(Exp[i],2))
                               result3[2,i]<- round(zh[i]*exposed_cases[i]/(total_cases[i]-exposed_cases[i]),2)
                               result3[3,i]<- round(LLR(exposed_cases[i],total_cases[i],zh[i]),3)
  if(absorb[total_cases[i]]<=total_cases[i]){result3[4,i]<- round(LLR(absorb[total_cases[i]],total_cases[i],zh[i]),3)}else{result3[4,i]<- paste("NA")} 
                               }
  if(posi<=length(total_cases)){result3[4,posi:length(total_cases)]<- paste("NA")}

                                              }else{# related to 11
                               result3<- as.data.frame(matrix(0,4,1))
                               names(result3)<- paste("Test",length(total_cases))
                               rownames(result3)<- c("Expected Cases","Estimated RR","LLR","CV")
                               result3[1,1]<- round(total_cases[length(total_cases)]/(1+zr),2)
                               result3[2,1]<- round(zh[1]*exposed_cases[length(exposed_cases)]/(total_cases[length(total_cases)]-exposed_cases[length(exposed_cases)]),2)
                               result3[3,1]<- round(LLR(exposed_cases[length(exposed_cases)],total_cases[length(total_cases)],zh[1]),3)
 if(absorb[total_cases[length(total_cases)]]<=total_cases[length(total_cases)]){result3[4,1]<- round(LLR(absorb[total_cases[length(total_cases)]],total_cases[length(total_cases)],zh[1]),6)}else{result3[4,1]<- paste("NA")} 

                                                   }# close 11

result<- data.frame(matrix(0,length(totalg)+1,12))
result[2:(length(totalg)+1),1]<- seq(1,length(totalg),1)
colnames(result)<- c(" "," "," ", "-----","Cumulative","--------"," "," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","Cases","Controls","Cases","Controls","E[Cases]","RR","LLR","target","actual","CV","H0 rejected")
zeros<- sum(totalg==0)
result[2:(length(totalg)+1),2]<- as.numeric(cases) ; result[2:(length(totalg)+1),3]<- as.numeric(controls) ; result[2:(length(totalg)+1),4]<- t(as.numeric(cases)%*%(upper.tri(matrix(0,length(cases),length(cases)),diag=T)*1)) 

## Useful auxiliar variable
Exp<- rep(0,test)
for(i in 1:test){Exp[i]<- sum(1/(1+zhc[1:total_cases[i]])) }  

result[2:(length(totalg)+1),5]<- t(as.numeric(controls)%*%(upper.tri(matrix(0,length(controls),length(controls)),diag=T)*1)) ; result[2:(length(totalg)+1),6]<- paste(round(Exp,2))
result[(zeros+2):(length(totalg)+1),7]<- t(result3[2,]); if(zeros>0){result[2:(zeros+1),7]<- paste("NA") } 
result[(zeros+2):(length(totalg)+1),8]<- t(result3[3,])
result[(zeros+2):(length(totalg)+1),9]<- t(result1[2,]); for(i in 2:nrow(result)){if(totalg[i-1]>N){result[i,9]<- alpha1}; if(totalg[i-1]==0){result[i,8]<- paste("NA") ; result[i,9]<- paste("NA")}}

result[2:nrow(result),10]<- paste(round(as.numeric(sum_sa_used),4)) 

result[(zeros+2):(length(totalg)+1),11]<- t(result2[3,]) ; if(zeros>0){result[2:(zeros+1),11]<- paste("NA")}
result[2:(length(totalg)+1),12]<- paste("No") ; if(rejt>0){result[(rejt+1):(length(totalg)+1),12]<- paste("Yes")}else{if(max(exposed_cases)>=max(absorb[total_cases])){result[(length(totalg)+1),12]<- paste("Yes")}}  

######## Graphic 1

par(mfrow=c(2,2))

## Here we define the limits to start the plot with respect to the x-axis. If we have zeros for the first entries of total_cases and rejected H0 previously, the plot will consider only the valid tests 
if(rejt==0){limx<- length(totalg)}else{limx<- rejt}


plot(seq(1,length(totalg),1),rep(max(absorb[total_cases]),length(totalg)),col="white",pch=18,xlab="Test",xaxt="n",cex.main=1.3,cex.lab=1.3,cex.main=1.5,
ylab=c("Cumulative Cases"),main=title,ylim=c(0,max(absorb[total_cases]+5,max(exposed_cases)))) 

if(length(exposed_cases)==length(total_cases)){points(seq(1,length(totalg),1),exposedg,col="blue",pch=20)
                                               lines(seq(1,length(totalg),1),exposedg,col="blue",lty=1)
                                              }else{points(length(totalg),max(exposedg),col="blue",pch=20)}

if(sum(absorb[total_cases]<=total_cases)>0){
                                                              ini<- 1
                                                              while(absorb[total_cases[ini]]>total_cases[ini]&ini<=length(total_cases[ini])){ini<- ini+1}
                                                              absorb_g<- absorb[total_cases]
                                                              xv<- seq(sum(totalg==0)+1,length(totalg),1)[absorb_g<=total_cases]
                                                              
                                                              absorb_g<- absorb_g[absorb_g<=total_cases] 
                                                              if(ini<=limx){absorb_g<- absorb_g[xv<=limx & xv>=ini] ; xv<- xv[xv<=limx & xv>=ini]}                                                                                                                  
                                                              
                                                              points(xv,absorb_g,col="red",pch=20)
                                                              lines(xv,absorb_g,col="red",lty=2)                     
                
                                           }

sequencia<- seq(1,length(totalg),1)
rotulos<- seq(1,length(totalg),1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)

points(seq(1,length(totalg),1),result[2:nrow(result),6],pch=2,col="black")
lines(seq(1,length(totalg),1),result[2:nrow(result),6],lty=1,col="black")

legend("topleft",c("Needed to reject H0 (CV)","Observed","Expected = E[Cases]"),col=c("red","blue","black"),pch=c(18,20,2),lty=c(2,1,1),bty="n")


## Graphic 2

saux<- 0
if(sum(total_cases<=N)>0){saux<- max(sum_sa[total_cases[total_cases<=N]])}
if(saux>0|max(totalg)>N){

if(sum(total_cases<=N)>0){yvt<- rep(0,sum(total_cases>N)+ sum(sum_sa[total_cases[total_cases<=N]]>0))}else{yvt<- rep(0,sum(total_cases>N))}

inix<- sum(totalg==0)+2
if(rejt==0){iniy<- nrow(result)}else{iniy<- rejt+1}
x<- rep(0,length(yvt))
cont<- 1
for(i in inix:iniy){
if(as.numeric(result[i,9])>0){
yvt[cont]<- as.numeric(result[i,9])
x[cont]<- as.numeric(result[i,1])
cont<- cont+1
                             }

                   }

if(rejt==0){yv<- as.numeric(result1[3,]) ; if(length(totalg)!=length(total_cases)){yv<- c(rep(0,length(totalg)-length(total_cases)),yv)}}else{yv<- sum_sa_used[1:rejt]}
#yvt<- as.numeric(sum_sa[total_cases_aux]) ; if(length(totalg)!=length(total_cases_aux)){yvt<- c(rep(0,length(totalg)-length(total_cases_aux)),yvt)} ; if(rejt>0){yvt<- yvt[1:rejt]}

plot(seq(1,length(totalg),1),rep(max(yvt),length(totalg)),pch=18,col="white",ylab="Alpha spending",xlab="Test",xaxt="n",cex.lab=1.3, 
ylim=c(0,max(max(yvt),0.1)),sub=" ",font.sub=1,main="Alpha Spending")

lines(x,yvt,lty=2,col="red")
points(x,yvt,pch=18,col="red")

if(max(sum_sa_used)>0){
points(seq(1,length(yv),1),yv,pch=18,col="blue")
lines(seq(1,length(yv),1),yv,lty=1,col="blue")
                      }

axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)
legend("topleft",c("Target","Actual spent"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")

                             }else{
plot(seq(1,length(totalg),1),rep(alpha1,length(totalg)),pch=18,col="white",ylab="Alpha spending",xlab="Test",xaxt="n",cex.lab=1.3, 
sub=" ",font.sub=1,main="Alpha spending"); legend("top",c("Not Applicable"),col=c("red"),pch=c(18),lty=c(2),bty="n")
                         }


######### Graphic 3

estimated<- rep(0,length(totalg))
for(i in 1:length(totalg)){if(result[i+1,7]!="NA"&result[i+1,7]!="Inf"){estimated[i]<- as.numeric(result[i+1,7])}}
if(sum(estimated)!=0){xxv<- seq(1,length(totalg),1)[estimated!=0]}else{xxv<- 0} 

plot(seq(1,length(totalg),1),rep(max(estimated),length(totalg)),pch=18,col="white",ylab="Observed relative risk",xlab="Test",xaxt="n",cex.lab=1.3, 
ylim=c(0,max(max(estimated),1)),sub=" ",font.sub=1,main="Observed Relative Risk")

if(sum(xxv)>0){
points(xxv,estimated[xxv],pch=18,col="blue")
lines(xxv,estimated[xxv],lty=1,col="blue")
              }else{legend("top",c("Not Applicable"),col=c("red"),pch=c(18),lty=c(2),bty="n")}
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)

######### Graphic 4

zeros<- sum(totalg==0)
yy<- as.numeric(result[(2+zeros):nrow(result),8])

plot(seq(1,length(totalg),1),rep(max(yy,1),length(totalg)),pch=18,col="white",ylab="Observed LLR",xlab="Test",xaxt="n",cex.lab=1.3, 
ylim=c(0,max(1,max(yy))),sub= " ",font.sub=1,main="Log-likelihood ratio")
lines(seq(1+zeros,length(totalg),1),yy,lty=1,col="blue")
points(seq(1+zeros,length(totalg),1),yy,col="blue")
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)


## Here we finish the code by offering the answers inside result

print(result,right=TRUE,row.names=FALSE)
         
message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: N= ",N,", alpha= ",alpha,", zp= ",z," and M= ",M,"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)

## Updating information for future tests.
inputSetUp[2,1:length(sum_sa)]<- as.numeric(sum_sa)
inputSetUp[6,1:test]<- sum_sa_used
inputSetUp[7,1:length(absorb)]<- absorb
inputSetUp[1,8]<- rejt

               }# closes 11

                                                                         } # closes first if for test # 10
        

## Updating information for future tests.
if(max(cases+controls)>0){inputSetUp[9:nrow(inputSetUp),]<- p;if(max(p)==0){base<- 1}else{base<- max(total_cases)+1}}
inputSetUp[1,7]<- base
write.table(inputSetUp,name)

setwd(safedir)

} ## end Analyze.Binomial function



 
