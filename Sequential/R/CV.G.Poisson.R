CV.G.Poisson <-
function(SampleSize,alpha=0.05,Looks,M=1){

#----------ExactPoisson.R-------------------------------------------------------------------
# Function that calculates the LLR for a given observed (c) and expected (u) number of cases
#-------------------------------------------------------------------------------------------

LLR <- function(cc,uu) {
	if(cc<=uu) x=0
	if(cc>uu) x = (uu-cc) + cc*log(cc/uu)
	x
	}



# T = maximum length of surveillance, defined in terms of expected counts under H0
# alpha = desired alpha level
# CVstart = LLR/CV start value, a good guess reduces computing time
# MinCases = The minimum number of cases for which a signal is allowed to occur, default=1
# Group = Time between looks in a group sequential trial, must be greater than 0
# Late = Time at first look at the data, defined in terms of the expected counts under H0, default=0
# NOTE: Only one of MinCases and Late can be different from the default value
CVstart=3
MinCases<- M
T<- SampleSize
Group<- T/Looks


if(alpha<=0|alpha>0.5|is.numeric(alpha)==FALSE){stop("alpha must be a number greater than zero and smaller than 0.5.",call. =FALSE)}
if(is.numeric(MinCases)==FALSE|MinCases!=round(MinCases)|MinCases<1){stop("M must be a positive integer.",call. =FALSE)}
if(T<=0|is.numeric(T)==FALSE){stop("SampleSize must be a positive number.",call. =FALSE)}
if(is.numeric(Looks)==FALSE|Looks<1|Looks!=round(Looks)){stop("Looks must be a positive integer.",call. =FALSE)}


                         if(1-ppois(0,T)<alpha|qpois(1-alpha,T)<MinCases){
                            out<- LLR(MinCases,T)                         
                            return(out)
                                      
                                               }else{ ## number 1
#--------------------------------------------------------------####


Perro_I<- function(CV){


absorb = seq(length=imax,from=0,by=0)		# Contains the number of events needed at time mu[i] in order to reject H0
for(i in 1:imax)
	while(LLR(absorb[i],mu[i])<CV) 
		absorb[i]=absorb[i]+1;



imin<- 1
count<-0
while(count==0){if(absorb[imin]<MinCases){absorb[imin]<- MinCases;imin<- imin+1};if(absorb[imin]>=MinCases | imin>imax){count<- 1}}



p = seq(length=imax*(absorb[imax]+1), from=0, by=0)	# p[i,j] is the probability of having j-1 cases at time mu[i]
									# starting probabilities are all set to zero's
dim(p) = c(imax,absorb[imax]+1)				# i in 1:imax-1 is the rows and j in 1:imax is the column

# Calculating the first row in the p[][] matrix for which there is a chance to reject H0
# --------------------------------------------------------------------------------------

for(s in 1:absorb[1]) p[1,s]=dpois(s-1,mu[1])		# Probability of having s-1 cases at time mu[1]
p[1,absorb[1]+1]=1-ppois(absorb[1]-1,mu[1])			# probability of rejecting H0 at time mu[imin]


# Calculating the remaining rows in the p[][] matix
# -------------------------------------------------
if(imax>1){
for(i in 2:imax) {
	for(j in 1:absorb[i])					# This loop calculates the p[][] matix, one column at a time, from left to right
		for(k in 1:min(j,absorb[i-1]))
			p[i,j]=p[i,j]+p[i-1,k]*dpois(j-k,mu[i]-mu[i-1])	# Calculates the standard p[][] cell values
	for(k in 1:absorb[i-1]) p[i,absorb[i]+1]=p[i,absorb[i]+1]+p[i-1,k]*(1-ppois(absorb[i]-k,mu[i]-mu[i-1]))
} # end for i	
          }

# Sums up the probabilities of absorbing states when a signal occurs, to get the alpha level
# ------------------------------------------------------------------------------------------

alpha=0
for(i in 1:imax) alpha=alpha+p[i,absorb[i]+1]					
time=0
for(i in 1:imax) time=time+i*Group*p[i,absorb[i]+1]
signaltime=time/alpha					
return(alpha)
}# end function P_erro_I
################################################################


#T=10
#alpha=0.05
#CVstart=3
#Group=1
#RR=2

ALPHAFIX=alpha
PRECISION=0.00000001

alphacont<- 0
cont<- 0
aux<- 2
CVold=0
CVnew<- CVstart			#Smart start values has little effect on computing time.
alphaold=1
CV=CVstart				#Smarter start vaules reduces computing time.
alpha=0
teste<- 0
hist<- c(CVold,alpha)

CV1<- CVold
CV2<- CV


imax = ceiling(T/Group)
					# The number of tests performed, including final time T.
mu = seq(length=imax,from=Group,by=Group)		# An array of the expected counts at each of the tests
								# mu[i] is the commulative expected count at the i'th test
mu[imax]=T							# Sets the expected count at the last test to equal T

loops=0
while(abs(CV1-CV2)>PRECISION & abs(alpha-ALPHAFIX)>PRECISION) {
loops=loops+1
						# Just a counter, not really needed
if(CV>0){alpha<- Perro_I(CV)} 

if(alpha==alphaold){teste<-1}

# Estimates a new critical value in the search of the one that gives alpha=0.05
# ---------------------------------------------------------------------------------


if(teste==0){

CVnew = CV - (alpha-ALPHAFIX)*(CV-CVold)/(alpha-alphaold)
alphaold=alpha
CVold=CV
CV=CVnew

hist<- rbind(hist,c(CVold,alpha))

            }else{if(CV>0){
if(alpha>ALPHAFIX & alphaold>ALPHAFIX){ # if 1
                                       CV1<- max(CV, CVold)
                                       CV<- CV1+abs(CV-CVold)+(loops-1)*0.1
                                       CVold<- CV1
                                       alpha<- Perro_I(CV)
                                       hist<- rbind(hist,c(CV,alpha))
                                      } # end if 1
if(alpha<ALPHAFIX & alphaold<ALPHAFIX){ # if 2
                                       CV2<- min(CV, CVold)
                                       CV<- CV2-abs(CV-CVold)-(loops-1)*0.1
                                       CVold<- CV2
                                       alpha<- Perro_I(CV)
                                       hist<- rbind(hist,c(CV,alpha))
                                      } # end if2
                          } # END IF(CV>0)


if((alpha<ALPHAFIX & alphaold>ALPHAFIX) | (alpha>ALPHAFIX & alphaold<ALPHAFIX) | CV<=0){ # if 3

CV1<- max(0,min(CV,CVold))
CV2<- max(CV, CVold)

while(abs(CV1-CV2)>PRECISION & abs(alpha-ALPHAFIX)>PRECISION){

                                                              CV<- (CV1+CV2)/2
                                                              alpha<- Perro_I(CV) 
                                                              hist<- rbind(hist,c(CV,alpha))
                                                              if(alpha<ALPHAFIX){CV2<- CV}else{CV1<- CV} 

                                                             }# end while(abs(CV1-CV2)>PRECISION & abs(alpha-ALPHAFIX)>PRECISION)
                                                                               } # end if 3


                 }# end else if(teste==0)

                

} #end while (abs(CV1-CV2)>PRECISION & abs(alpha-ALPHAFIX)>PRECISION)



out<- matrix(0,1,4)
if(teste==1){
hist<- hist[-1,]
CV2<- round(min(hist[hist[,2]<ALPHAFIX,1]),10)+0.000001
CV1<- round(max(hist[hist[,2]>ALPHAFIX,1]),10)
alpha1<- min(hist[hist[,2]>ALPHAFIX,2])
alpha2<- max(hist[hist[,2]<ALPHAFIX,2])


out[1,1:2]<- c(CV1,CV2)
out[1,3:4]<- c(alpha1,alpha2)

            }else{out[1,1:2]<- c(CV,alpha)}
if(length(out)>1){out<- out[2]}
return(out)
#return(hist[-1,])
                                            } ## close number 1
#--------------------------------------------------------------####

}


#CV.G.Poisson(SampleSize=1,alpha=0.05,Looks=1000,M=10)

