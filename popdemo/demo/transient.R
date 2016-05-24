#=================================================
#Transient dynamics
#References:
#Stott et al. 2011 Ecol. Lett., 14, 959-970
#(and citations therein)
#=================================================
#
#
#
#This matrix is for the desert tortoise
#(Gopherus Agassizii) with medium fecundity
#(Doak et al. 1994 Ecol. Appl., 4, 446-480).
Tort<-Matlab2R("[0 0 0 0 0 1.3 1.98 2.57;0.716 0.567 0 0 0 0 0 0;0 0.149 0.567 0 0 0 0 0;0 0 0.149 0.604 0 0 0 0;0 0 0 0.235 0.56 0 0 0;0 0 0 0 0.225 0.678 0 0;0 0 0 0 0 0.249 0.851 0;0 0 0 0 0 0 0.016 0.86]")
#
#
#
#-------------------------------------------------
#TRANSIENT DYNAMICS OF SPECIFIC POPULATION 
#STRUCTURES
#-------------------------------------------------
#
#We can calculate amplified and/or attenuated 
#transient indices for specified population 
#vectors. All indices standardise the matrix 
#(removing effects of asymptotic dynamics) and the 
#population vector (removing the effects of initial 
#population size). This means transient indices are 
#comparable both within and between models, no 
#matter their dimension, rate of asymptotic growth, 
#or overall population size.
#
#A population that amplifies:
Tortamp<-Matlab2R("[1;1;2;3;5;8;13;21]")
#
#Calculate reactivity 
#(amplification in the first timestep)
amp1<-reactivity(Tort, vector=Tortamp)
amp1
#
#Calculate maximum amplification 
#(largest overall amplification)
amp2<-maxamp(Tort, vector=Tortamp)
amp2
#
#Calculate inertia 
#(asymptotic amplification)
amp3<-inertia(Tort, vector=Tortamp)
amp3
#
#Indices of amplification are always >1.
#
#A population that attenuates:
Tortatt<-Matlab2R("[21;13;8;5;3;2;1;1]")
#
#Calculate first-timestep attenuation 
#(attenuation in the first timestep)
att1<-firststepatt(Tort, vector=Tortatt)
att1
#
#Calculate maximum attenuation
#(smallest overall attenuation)
att2<-maxatt(Tort, vector=Tortatt, return.t=TRUE)
att2
#
#Calculate inertia 
#(asymptotic attenuation)
att3<-inertia(Tort, vector=Tortatt)
att3
#
#Indices of attenuation are always <1.
#
#It is also possible to return the realised 
#population size (including effects of asymptotic 
#growth and initial population size):
att1.2<-firststepatt(Tort, vector=Tortatt, 
                   return.N=TRUE)
att1.2
att2.2<-maxatt(Tort, vector=Tortatt, return.t=TRUE, 
             return.N=TRUE)
att2.2
att3.2<-inertia(Tort, vector=Tortatt, 
              return.N=TRUE, t=50)
att3.2
#
#We can visualise both standardised and non-
#standardised transient dynamics on a population
#projection plot.
#Standardised dynamics:
stdpr<-project(Tort, vector=Tortatt, time=50,
               standard.A=TRUE, standard.vec=TRUE)
plot(stdpr, log="y", main="Standardised transient dynamics")
points(1, att1, col="red", pch=3)
text(1, att1, "First-timestep\nattenuation",
     col="red", adj=c(0,0))
points(att2$t, att2$maxatt, col="red", pch=3)
text(att2$t, att2$maxatt, "Maximum\nattenuation",
     col="red", adj=c(0,0.5))
points(50, att3, col="red", pch=3)
text(50, att3, "Inertia", col="red", adj=c(1,1))
#
#Non-standardised dynamics:
nstdpr<-project(Tort, vector=Tortatt, time=50)
plot(nstdpr, main="Non-standardised population dynamics")
points(1, att1.2$N, col="red", pch=3)
text(1, att1.2$N, "First-timestep\nattenuation",
     col="red", adj=c(0,0))
points(att2.2$t, att2.2$N, col="red", pch=3)
text(att2.2$t, att2.2$N, "Maximum\nattenuation",
     col="red", adj=c(0,0.5))
points(50, att3.2$N, col="red", pch=3)
text(50, att3.2$N, "Inertia", col="red", adj=c(1,1))
#
readline("Hit <Enter> for next topic, or <Esc> to cancel")
#
#
#
#3.2----------------------------------------------
#BOUNDS ON TRANSIENT DYNAMICS
#-------------------------------------------------
#
#Bounds on transient dynamics result from stage-
#biased population projections. In practise, the 
#stage-biased vectors are standard basis vectors: 
#column vectors consisting of zeroes except for a 
#1 in one of the elements.  For example, a set of 
#stage-biased vectors for a 3 by 3 matrix is 
#c(1,0,0); c(0,1,0); c(0,0,1).
#
#Transient bounds can be useful to use to inform 
#on the potential range of transient dynamics and/or 
#can be used when the structure of the population is 
#unknown. Every PPM has a full set of amplified bounds 
#and a full set of attenuated bounds, resulting from 
#biases of different stages.
#
#To calculate bounds on transient dynamics, as with 
#projecting transient dynamics, do not specify a
#population vector.
#
#Calculate the bounds on reactivity and first-timestep
#attenuation:
reac<-reactivity(Tort)
reac
firstep<-firststepatt(Tort)
firstep
#
#Calculate the bounds on maximum amplification and
#maximum attenuation:
mxamp<-maxamp(Tort, return.t=TRUE)
mxamp
mxatt<-maxatt(Tort, return.t=TRUE)
mxatt
#
#Calculate the upper and lower bounds on inertia:
upinertia<-inertia(Tort, bound="upper")
upinertia
lowinertia<-inertia(Tort, bound="lower")
lowinertia
#
#The same arguments are available to return the
#realised (non-standardised) population size
#using return.N=T.
#
#Visualising these on a population projection:
SBpr<-project(Tort, time=50, standard.A=TRUE)
plot(SBpr, log="y")
points(c(1,1), c(reac,firstep), col="red", pch=3)
text(c(1,1), c(reac,firstep), col="red",
     c("Reactivity","First-timestep\nattenuation"),
     adj=c(-0.2,0.5))
points(c(mxamp$t,mxatt$t), c(mxamp$maxamp,mxatt$maxatt),
       pch=3, col="red")
text(c(mxamp$t,mxatt$t), c(mxamp$maxamp,mxatt$maxatt),
     c("Maximum amplification","Maximum attenuation"),
     col="red", adj=c(-0.2,0.5))
points(c(35,35), c(upinertia,lowinertia), col="red", pch=3)
text(c(35,35), c(upinertia,lowinertia), col="red",
       c("Upper\ninertia","Lower\ninertia"), adj=c(-0.2,0.5))
#
#
#
#END
