
#----------ExactBinomial.R-------------------------------------------------------------------

# Version of Feb/2015

# -------------------------------------------------------------------------------------------
# Function produces critical value for the continuous Sequential Binomial MaxSPRT
# -------------------------------------------------------------------------------------------

SampleSize.Binomial<- function(RR,alpha=0.05,power=0.9,M=1,z=1){

MinCases<- M

# alpha = desired alpha level
# MinCases = The minimum number of cases for which a signal is allowed to occur, default=1
# z = matching ratio between exposed and unexposed cases   
# RR is the relative risk
if((RR<=1|is.numeric(RR)==FALSE)){stop("'RR' must be a number greater than 1.",call. =FALSE)}
if(MinCases<1){stop("'M' must be a positive integer.",call. =FALSE)}
if(MinCases!=round(MinCases)){stop("'M' must be a positive integer.",call. =FALSE)}
if(alpha<=0|alpha>0.5|is.numeric(alpha)==FALSE){stop("'alpha' must be a number greater than zero and smaller than 0.5.",call. =FALSE)}
if(power<alpha|power>=1|is.numeric(power)==FALSE){stop("'power' must be a number greater than or equal to 'alpha' and smaller than 1.",call. =FALSE)}
if(is.numeric(MinCases)==FALSE|MinCases!=round(MinCases)|MinCases<1){stop("'M' must be a positive integer.",call. =FALSE)}
if(z<=0|is.numeric(z)==FALSE){stop("'z' must be greater than zero.",call. =FALSE)}

Nr<- 1
while(1-pbinom(Nr-1,Nr,1/(1+z))>alpha){Nr<- Nr+1}
N1<- max(Nr,M)

N2<- N1
pow<- 0
M<- as.integer(M)

# Finding a bound for N
while(pow< power){
   N2<- N2+100
   N2<- as.integer(N2)
   cv<- CV.Binomial(N=N2,alpha,M,z)[[1]] ; pow<- Performance.Binomial(N=N2,M,cv,z,RR)[[1]]   
   
                 }

# Finding the exact solution
tes<- 0
  while(N2-N1>1&tes==0){
   Nm<- round((N1+N2)/2)
   res1<- CV.Binomial(N=Nm,alpha,M,z)
   cv<- res1[[1]]
   res2<- Performance.Binomial(N=Nm,M,cv,z,RR)   
   pow<- res2[[1]]
   if(pow==power){tes<- 1;Ns<- Nm;pow1=pow}else{if(pow>power){N2<- Nm;pow1<- pow;Ns<- Nm}else{N1<- Nm}}   
                       }
   
                                 
# Ns is the solution

error1<- res1[[2]]

out<- list(Ns,cv,error1,pow1)
names(out)<- c("Required_N","cv","Type_I_Error","Actual_power")
return(out)

} #end function SampleSize.Binomial



