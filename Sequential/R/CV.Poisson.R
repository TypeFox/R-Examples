CV.Poisson<-
function(SampleSize,D=0,M=1,alpha=0.05){
alphin<- alpha
T<- SampleSize
L<- T
# ------------------- INPUT VARIABLE ----------------------------------------------------------
# T = maximum length of surveillance, defined in terms of expected counts under H0
# alpha = desired alpha level
# CVstart = CV start value, a good guess reduces computing time
# MinCases = The minimum number of cases for which a signal is allowed to occur
# Late = Time < T for first look at the data, defined in terms of the expected counts under H0

Late<- D
MinCases<- M

##### SOME PRE-CALCULATED CV's TO SAVE EXECUTION TIME, GIVEN M, T, ALPHA=0.05, AND D=0.

TableCV.PoissonM<-data.frame(matrix(c(2.85393700326128,2.36663822744395,1.7742178463146,0,0,0,0
,2.96497144429574,2.57638960908667,2.15070686609869,1.68320857573985,0,0,0
,3.0469774628286,2.68935364101625,2.34967932492501,2.00015758807598,0,0,0
,3.11041922987987,2.77748334336629,2.47487306554246,2.18732770846946,0,0,0
,3.16210581650405,2.8493269885449,2.56531970873171,2.31713928988974,1.76648463065627,0,0
,3.24500386808042,2.93740985148835,2.69918199672857,2.49889219887528,2.08947296098257,1.56463589358817,0
,3.29718279728549,3.01290868171239,2.80395477091307,2.62366847225542,2.26759506082247,1.93644717718729,0
,3.34272874946839,3.08209924178235,2.87390384946275,2.69934978606631,2.40680950082531,2.09383466767681,1.74055087831132
,3.41378224992083,3.17006191246864,2.98556046816797,2.82925856210701,2.57262664682324,2.33777083613016,2.08603188356668
,3.46795158033368,3.23800870652229,3.06424796686384,2.92156132978567,2.69058633468102,2.4848340056375,2.2814406139523
,3.51174870474631,3.29055057983813,3.12525250048283,2.99310550216978,2.7814346207046,2.58938849090246,2.41540178483103
,3.56259057516111,3.35326528420839,3.19995342291952,3.07561262580597,2.87793940574268,2.71199552491449,2.55663428340283
,3.62812301587201,3.4301407625587,3.28821566346814,3.17637032177521,2.99779165442139,2.84685765103372,2.71713745193818
,3.67631995875227,3.48796061486674,3.3566769512246,3.24963354903075,3.0810513452118,2.94727030411798,2.82771061367823
,3.71576438223103,3.53415028882636,3.40671487050728,3.3071350801791,3.14780120177428,3.01963940268373,2.91122208653942
,3.77466349620378,3.60505631658695,3.48595979980634,3.39197395348475,3.24661941308076,3.13049529955633,3.03073524477903
,3.81990269391777,3.65714243579444,3.54482617261637,3.4555214311217,3.31795545593698,3.21042821525469,3.11755331874543
,3.85575466209975,3.69888480478571,3.59056669774142,3.50521969535492,3.37419401277773,3.27148598919451,3.18419630429743
,3.91085346466642,3.7624739866046,3.65993869307721,3.58089967455607,3.45808683624598,3.36288832681964,3.28403027659113
,3.95232090599911,3.81014086931301,3.71199319899898,3.63650814981748,3.52008147061398,3.43006496056302,3.35579446678864
,3.98557719058308,3.84774767167223,3.75332945022298,3.68058354716942,3.56867928724568,3.48296604520948,3.41123480524391
,4.02533791597334,3.89271466068751,3.80241164747275,3.73238627490875,3.62614993442272,3.54430757139021,3.4766545706701
,4.07482766406364,3.94892970910266,3.86276222895003,3.79683534705039,3.69651138482233,3.61982463789605,3.55679903763937
,4.11223437028032,3.990901126178,3.90806497513453,3.84484673513135,3.74875734849489,3.67570325845673,3.61551308623033
,4.14213356396756,4.02415316719417,3.94413519839946,3.88271021939305,3.79014311617218,3.71945215752241,3.66182960350214
,4.18803098306215,4.07529664413041,3.99894986501106,3.94056304984504,3.85265830465641,3.78592974865935,3.73152355787215
,4.22263170622893,4.11369239777545,4.04002057942574,3.98377804171546,3.899238807548,3.83526454689286,3.78312575686608
,4.25030952971488,4.1443167419865,4.07263799414012,4.01808951759544,3.93617480815424,3.87418262220993,3.82390844776879
,4.29282949701916,4.19116679996068,4.12255869710559,4.07046611590226,3.9922717187982,3.93336365702725,3.88559998964833
,4.32491732176988,4.22641197548099,4.16002223662034,4.10966536848407,4.03421016321423,3.97745325368572,3.93152895644531),ncol=7,byrow=T))

###### END TABLE

##### SOME PRE-CALCULATED CV's TO SAVE EXECUTION TIME, GIVEN D, T, ALPHA=0.05, AND M=1.

TableCV.PoissonD<-data.frame(matrix(c(0,0,0,0,0,0,0,0
,2.964971,1.683209,0,0,0,0,0,0
,3.046977,2.000158,0,0,0,0,0,0
,3.110419,2.187328,1.600544,0,0,0,0,0
,3.162106,2.317139,1.766485,0,0,0,0,0
,3.245004,2.498892,2.089473,1.84232,0,0,0,0
,3.297183,2.545178,2.267595,1.936447,1.611553,0,0,0
,3.342729,2.546307,2.40681,2.093835,1.921859,0,0,0
,3.413782,2.694074,2.572627,2.337771,2.211199,1.829011,0,0
,3.467952,2.799333,2.591675,2.484834,2.298373,2.087405,1.834621,0
,3.511749,2.88072,2.683713,2.589388,2.415402,2.254018,1.96566,1.755455
,3.562591,2.970411,2.794546,2.711996,2.556634,2.347591,2.203782,2.020681
,3.628123,3.082511,2.918988,2.846635,2.717137,2.542045,2.425671,2.260811
,3.67632,3.15949,3.011001,2.886783,2.827711,2.668487,2.527763,2.432668
,3.715764,3.223172,3.080629,2.963485,2.911222,2.765594,2.634068,2.553373
,3.774663,3.313966,3.186878,3.078748,3.030735,2.903286,2.789967,2.68473
,3.819903,3.381606,3.261665,3.162197,3.117553,2.99958,2.897811,2.802863
,3.855755,3.434748,3.320749,3.226113,3.162908,3.05147,2.978063,2.890933
,3.910853,3.515052,3.407923,3.321868,3.247872,3.15182,3.090356,3.019184
,3.952321,3.574091,3.47261,3.391376,3.321971,3.232345,3.155596,3.109251
,3.985577,3.620223,3.523446,3.445695,3.379278,3.294843,3.222053,3.177847
,4.025338,3.675035,3.583195,3.509028,3.446674,3.367227,3.298671,3.238461
,4.074828,3.742844,3.655984,3.587079,3.528662,3.454679,3.391821,3.336012
,4.112234,3.792978,3.710128,3.644349,3.588871,3.518954,3.459256,3.406929
,4.142134,3.832686,3.752749,3.689355,3.636272,3.568952,3.512138,3.462111
,4.222632,3.938105,3.835265,3.808087,3.760123,3.700033,3.649189,3.605012
,4.25031,3.97371,3.874183,3.847892,3.801678,3.743656,3.694832,3.652326
,4.292829,4.028089,3.933364,3.887512,3.864597,3.809685,3.763627,3.723608
,4.324917,4.047191,3.977453,3.931529,3.911308,3.858669,3.814122,3.776275),ncol=8,byrow=T))

###### END TABLE


colnames(TableCV.PoissonD)<- c(0,1,2,3,4,6,8,10)
colnames(TableCV.PoissonM)<- c(1,2,3,4,6,8,10)

TV<- c(1,1.5,2,2.5,3,4,5,6,8,10,12,15,20,25,30,40,50,60,80,100,120,150,200,250,300,400,500,600,800,1000)
MV<- c(1,2,3,4,6,8,10)
DV<- c(0,1,2,3,4,6,8,10)

tDescont<- matrix(c(5,10,10,10,15,20,60,60,80,800,1000,1000,1,2,4,8,10,3,4,6,8,3,1,8),nrow=2,byrow=T)
# ------------------------------------------------------------

####### Tests to verify the validity of the chosen parameters

teste1<- 0

if(T<=0){teste1<- 1; out<- c("SampleSize must be > 0")}
if(teste1==0){if(alpha>0.5 | alpha<(10^(-7))){teste1<- 1; out<- c("alpha must be a number in the (1e-7,0.5] interval")}}
if(teste1==0 & M>100){teste1<- 1; out<- c("M must be a positive integer in the range [1,100]")}
#if(teste1==0 & T>1000){teste1<- 1; out<- c("Use T<=1000")}


if(teste1==0){if(1-ppois(MinCases-1,T)<alpha){
                                teste1<- 1
                                T_min<- min(seq(T,M,0.01)[1-ppois(M-1,seq(T,M,0.01))>=alpha])
                                out<- list("Does not have solution. For this M and alpha, SampleSize must be >=",T_min,".")
                                             }
             }
if(Late>T & teste1==0){teste1<- 1; out<- c("D must be <= SampleSize") }
if(Late<0 & teste1==0){teste1<- 1; out<- c("Negative values for D does not make sense. Use 0<=D<=SampleSize.") }
if(M<1 & teste1==0){teste1<- 1; out<- c("M must be a positive integer in the range[1,100].") }

if(teste1==1){stop(out,call.=FALSE)}

####### end parameters validity tests
# ------------------------------------------------------------

teste2<- teste1

####### Choosing CVstart by using the known CV table to save time.
# -------------------------------------------------------------------------

x1<- 1
x2<- 0
while(x1<=ncol(tDescont)&x2==0){if(T==tDescont[1,x1]&D==tDescont[2,x1]){x2<- 1};x1<- x1+1}

if(x2==0){
if(alpha==0.05 & T>=2 & teste1==0){
               il<- sum(TV<=T)
               ir<- sum(MV<=MinCases)
               ir2<- sum(DV<=Late)
             
 if(sum(T==TV)>0 & sum(M==MV)>0){if(Late==0){teste1<- 1;out<- TableCV.PoissonM[il,ir]}else{if(M==1 & sum(Late==DV)>0 & TableCV.PoissonD[il,ir2]>0){teste1<- 1;out<- TableCV.PoissonD[il,ir2]}else{CVstart<- 3}}}else{
           il2<- 1+sum(TV<T) 
           if(Late==0){CVstart<- (TableCV.PoissonM[il2,ir]+TableCV.PoissonM[il,ir])/2}else{if((TableCV.PoissonD[il2,ir]+TableCV.PoissonD[il,ir2])/2>0){CVstart<- (TableCV.PoissonD[il2,ir]+TableCV.PoissonD[il,ir2])/2}else{CVstart<-3}}
                                                                                                               }
                                 }else{CVstart<- 3}
         }else{CVstart<- 3}

if(T>1000){CVstart<- 4}

####### end the choice of CVstart by using the known CV table to save time.  
# -------------------------------------------------------------------------


if(teste1==0){ # calculating CV only if the parameters tests have not found a warning message and the known CV table does not have the corresponding CV solution. 

#---------------------------------------------------------------------
# Function that calculates the product log through a recursive formula
#---------------------------------------------------------------------
ProdLog <- function(z){
	x = z-z^2+1.5*z^3-(8/3)*z^4+(125/24)*z^5-(54/5)*z^6+(16807/720)*z^7
	for(i in 1:10) x = x-(x*exp(x)-z)/(exp(x)+x*exp(x))
	x
	                } # end ProdLog function 

#----------------------------------------------------------------------------------------------
# Function that calculates the probability of type I error for a given set of IMPUT parameters
#----------------------------------------------------------------------------------------------
Perror_I<- function(cv){

z = -exp(-1-cv/c)
mu = -c * ProdLog(z) 		#The expected counts under H0 that is needed to reject the null with i number of adverse events
mtemp = c(0,mu)
mmu = diff(mtemp) 		#The marginal difference of the mu[] vector


imin=MinCases
while (mu[imin] < Late) imin=imin+1
if(imin>MinCases) { 
	mu[imin-1]=Late
	mmu[imin]=mu[imin]-Late
	} #END if imin>MinCases

imax=1
while (mu[imax] < T) imax=imax+1    		# imax is the maximum number of cases that will generate a signal.            

# NOTE: If imax=1, this code will not work


if(imin<imax){

# Defining the p[][] matrix
# -------------------------

p = seq(length=(imax-1)*imax, from=0, by=0)				# p[i,j] is the probability of having j-1 cases at time mu[i]
dim(p) = c(imax-1,imax)								# i in 1:imax-1 is the rows and j in 1:imax is the column

# Calculating the first row p[imin][] in the matrix for which there is a chance to reject H0
# When MinCases=1, there is no skipping, and it is the first row in the matrix (p[1][]).
# ------------------------------------------------------------------------------------------

if(imin==MinCases) {
	for(s in 1:imin) p[imin,s] = dpois(s-1,mu[imin])			# Probability of having s-1 cases at time mu[imin], not rejectinh H0
	p[imin,imin+1] = 1-ppois(imin-1,mu[imin])				# Probability of having s+ cases at time mu[imin], rejectinh H0
	} # end if

if(imin>MinCases) {
	for(s in 1:imin) p[imin-1,s]=dpois(s-1,mu[imin-1])		# Probability of having s-1 cases at time mu[imin-1], not rejecting H0
	p[imin-1,imin+1] = 1-ppois(imin-1,mu[imin-1])				# Probability of having s+ cases at time mu[imin-1], rejecting H0
	for(s in 1:imin) 								# Probability of having s-1 cases at time mu[imin], not rejectinh H0
		for(k in 1:s) 
			p[imin,s]=p[imin,s]+p[imin-1,k]*dpois(s-k,mmu[imin])	
	for(k in 1:imin) 
		p[imin,imin+1] = p[imin,imin+1] + p[imin-1,k]*(1-ppois(imin-k,mmu[imin]))
} # end if 



# Calculating the remaining rows in the p[][] matix
# -------------------------------------------------

if(MinCases+1<=imax-1&((imin+1)<=(imax-1)))
probaux3<- 0
i<- (imin+1)
while(i <=(imax-1)&probaux3<=alphin+PRECISION) {
	for(j in 1:(i-1))							# This loop calculates the p[][] matix, one column at a time, from left to right
		for(k in 1:j) 
			p[i,j]=p[i,j]+p[i-1,k]*dpois(j-k,mmu[i])	# Calculates the standard p[][] cell values
	for(k in 1:(i-1))
		p[i,i]=p[i,i]+p[i-1,k]*dpois(i-k,mmu[i])		# Calculates the diagonal under the absorbing states, which requires a unique formula
	for(k in 1:(i-1)) 
		p[i,i+1]=p[i,i+1]+p[i-1,k]*(1-ppois(i-k,mmu[i]))# Calculates the diagonal absorbing states where H0 is rejected
 probaux3<- probaux3 + p[i,i+1]
i<- i+1
} # end for i	


pp=0
if(imax>imin)
for(k in 1:(imax-1)) pp=pp+p[imax-1,k]*(1-ppois(imax-k,T-mu[imax-1])) #Calculates the last probability to signal before time T


# Sums up the probabilities of absorbing states when a signal occurs, to get the alpha level
# ------------------------------------------------------------------------------------------

alpha=0
if(imin>MinCases) alpha=p[imin-1,imin+1]
for(i in imin:(imax-1)) alpha=alpha+p[i,i+1]					
alpha=alpha+pp


}else{alpha<- 1-ppois(imax-1,mu[imax])} # end if(imin<imax)



                      } # end Perror_I
#############################################################################################




###### Here starts the numerical procedure to find CV as a solution for a type I error equal to alpha 
#----------------------------------------------------------------------------------------------------

ALPHAFIX=alpha
PRECISION=0.00000001

LLRold=0				#Smart start values has little effect on computing time.
alphaold=1
LLR=CVstart				#Smarter start vaules reduces computing time.
alpha=0

c = 1:(2*T)
CV1<- LLRold
CV2<- LLR
teste<- 0

loop=0
while(abs(CV1-CV2)>PRECISION & abs(alpha-ALPHAFIX)>PRECISION) { #1
loop=loop+1	

if(loop>100 | alpha==alphaold | LLR<0 ){teste<- 1}								

if(teste==0){alpha<- Perror_I(LLR)} 

# Searching for the critical value that gives the exact alpha level. It tries first a linear interpolation, but, if it does not works, it is used the bisection method
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if(teste==0){ # this if() serves to verifies if the linear interpolation has failed 

LLRnew = LLR - (alpha-ALPHAFIX)*(LLR-LLRold)/(alpha-alphaold)
alphaold=alpha
LLRold=LLR
LLR=LLRnew
            }else{ # starting the bisection method in case of the linear interpolation has failed

CV1<- 0
CV2<- 20


                  while(abs(CV1-CV2)>PRECISION & abs(alpha-ALPHAFIX)>PRECISION){#2

                                                              LLR<- (CV1+CV2)/2
                                                              alpha<- Perror_I(LLR)                                                               
                                                              if(alpha<ALPHAFIX){CV2<- LLR}else{CV1<- LLR}  
                                                               

                                                                               }# end while(abs(CV1-CV2)>PRECISION & abs(alpha-ALPHAFIX)>PRECISION) #2                                                                             


                 }# end else if(teste==0)       
                                                             } #end while (abs(CV1-CV2)>PRECISION & abs(alpha-ALPHAFIX)>PRECISION) #1


## The outcome of this program will be the CV, or, in case of not being possible an exact solution, it will return a table with conservative and liberal critical values.

if(abs(alpha-ALPHAFIX)>PRECISION){
teste2<- 1
CV2<- ceiling((10^6)*CV2)/(10^6)
CV1<- floor((10^6)*CV1)/(10^6)
alpha1<- Perror_I(CV1)
alpha2<- Perror_I(CV2)

out<- matrix(c(CV2,alpha2,CV1,alpha1),ncol=2,byrow=F)
out<- as.data.frame(out)
colnames(out)<- c("Conservative","Liberal")
rownames(out)<- c("CV","alpha")
                                 }else{out<- LLR}

} # end if(if teste1==0)

if(teste2==0 & teste1==0 | teste2==0 & teste1==1){return(out)}
if(teste2==1 & teste1==0){
message("For these parameters there is no exact alpha level.")
message("Follows a conservative critical value:")
message("--------------------------------")
return(out[1,1])
                         }                         

}

#system.time(res<- CV.Poisson(SampleSize=1001,D=0,M=1,alpha=0.05))


