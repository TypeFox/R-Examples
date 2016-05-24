SampleSize.Poisson <-
function(alpha=0.05,power=0.9,M=1,D=0,RR=2,precision=0.000001)
{  

teste1<- 0
MinCases<- M
Late<- D
####### Tests to verify the validity of the chosen parameters
if(teste1==0){if(alpha>0.5 | alpha<(10^(-7))){teste1<- 1; out<- c("alpha must be a number in the (1e-7,0.5] interval")}}
if(RR<1 & teste1==0){teste1<- 1; out<- c("RR must be >=1.") }
if(teste1==0 & M>100){teste1<- 1; out<- c("M must be a positive integer in the range [1,100]")}
if(Late<0 & teste1==0){teste1<- 1; out<- c("Negative values for D does not make sense. Use D>=0.") }
if(M<1 & teste1==0){teste1<- 1; out<- c("M must be a positive integer in the range[1,100].") }
if(Late<0 & teste1==0){teste1<- 1; out<- c("Negative values for D does not make sense. Use 0<=D<=T.") }

#---------CODE TO CALCULATE POWER, SIGNAL TIME AND SURVEILLANCE TIME FOR GIVEN T
#######--------------------------------------------------------------------------

faux<-
function(L=30,D=0,M=1,RR=1,alpha=0.05){

# ------------------- INPUT VARIABLE ----------------------------------------------------------
# L = maximum length of surveillance, defined in terms of expected counts under H)
# RR = relative risk, RR=1 corresponds to H0
# M = The minimum number of cases for which a signal is allowed to occur
# D = Time < T for first look at the data, defined in terms of the expected counts under H0
# alpha = significance level

####### Tests to verify the validity of the chosen parameters
T<- L
teste1<- 0
MinCases<- M
Late<- D

if(T<=0){teste1<- 1; out<- c("T must be > 0")}
if(teste1==0){if(alpha>0.5 | alpha<(10^(-7))){teste1<- 1; out<- c("alpha must be a number in the (1e-7,0.5] interval")}}
if(RR<1 & teste1==0){teste1<- 1; out<- c("RR must be >=1.") }
if(teste1==0 & M>100){teste1<- 1; out<- c("M must be a positive integer.")}


if(M<1 & teste1==0){teste1<- 1; out<- c("M must be a positive integer in the range[1,100].") }

# If the parameters are incorrect in any sense, the code is interrupted and an error message is informed according to the possibilies above
#------------------------------------------------------------------------------------------------------------------------------------------
if(teste1==1){stop(out,call.=FALSE)}

## calculates the critical value

cv<- CV.Poisson(SampleSize=T,D,M,alpha)

## calculates power and signal time for relative risk equal to 2

if(length(cv)==1){
PT<- Performance.Poisson(SampleSize=T,D,M,cv,RR)
                }else{

cvc<- cv[1,1]
cvl<- cv[1,2]
resc<- Performance.Poisson(SampleSize=T,D,M,cvc,RR)
resl<- Performance.Poisson(T,D,M,cvl,RR)
powerc<- resc[1]
powerl<- resl[1]
signaltimec<- resc[2]
signaltimel<- resl[2]
SurveillanceTimec<- resc[3]
SurveillanceTimel<- resl[3]

PT<- matrix(c(powerc,powerl,signaltimec,signaltimel,SurveillanceTimec,SurveillanceTimel),ncol=2,byrow=T)
rownames(PT)<- c("Power","Signal Time","Surveillance Time")
colnames(PT)<- c("Conservative","Liberal")
                     }



# Output assigned as a vector
# ---------------------------


out=list(cv,PT)
names(out)<- c("CV","Power.SignalTime")
return(out)

}

##-------------------------------------------------------------------------------
######---------------------------------------------------------------------------

T1<- max(qexp(alpha,1),D)
T_min<- min(seq(0.001,M,0.01)[1-ppois(M-1,seq(0.001,M,0.01))>=alpha])
T1<- max(T1,T_min)
if(T1<30){T2<- 30}else{T2<- T1+30}
result<- faux(T2, D, M, RR,alpha)
pow<- result$Power.SignalTime[1]
while(pow<power){T2<- T2+30;result<- faux(T2, D, M, RR,alpha);pow<- result$Power.SignalTime[1]}
Told<- T2
Tm<- (T1+T2)/2
result<- faux(Tm, D, M, RR,alpha)
pow<- result$Power.SignalTime[1]
cont<- 0
poder<- matrix(0,31,1)
lim<- log((T2-T1)/precision)/log(2)+1
while((pow<power|power+precision<pow)&cont<lim){

                                       if(pow>power){T2<- Tm;Told<- Tm}else{T1<- Tm}
                                       Tm<- (T1+T2)/2; cont<- cont+1;result<- faux(Tm, D, M, RR,alpha); pow<- result$Power.SignalTime[1]
                                       poder[cont,1]<- pow
                                       
                                      }
L<- list(Tm)
names(L)<- c("SampleSize")
return(L)

}
