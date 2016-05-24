#=================================================
#Population projections
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
#PROJECTING SPECIFIC POPULATION STRUCTURES
#-------------------------------------------------
#
#Create a population vector
Tortvec<-Matlab2R("[1;1;2;3;5;8;13;21]")
#Create a population projection of that vector
pr1<-project(Tort, vector=Tortvec, time=50)
pr1
#
#The first entry is at time=0, so the nth
#entry of the projection is the n-1th
#projection interval.
#
#We can also return the time series of population
#vectors, showing size or density of each 
#individual stage class over time:
pr1.2<-project(Tort, vector=Tortvec, time=10,
               return.vec=TRUE)
pr1.2
#
#We can standardise aspects of the projection,
#for example removing asymptotic effects by
#standardising the matrix (which scales the 
#matrix by lambda-max):
pr2<-project(Tort, vector=Tortvec, time=50,
             standard.A=TRUE)
pr2
#
#...and standardising the population structure
#by population size so overall density is 
#initially equal to 1:
pr3<-project(Tort, vector=Tortvec, time=50,
             standard.A=TRUE, standard.vec=TRUE)
pr3
#
#We can plot these projections easily
plot(pr1)
plot(pr2, main="Standardised matrix")
plot(pr3, main="Standardised matrix and population vector")
#
readline("Hit <Enter> for next topic, or <Esc> to cancel")
#
#
#
#2.2----------------------------------------------
#PROJECTING STAGE-BIASED DYNAMICS
#-------------------------------------------------
#
#Stage-biased population vectors give the
#most extreme possible transient dynamics.
#In practise, the stage-biased vectors are
#standard basis vectors: column vectors
#consisting of zeroes except for a 1 in one
#of the elements.  For example, a set of stage-
#biased vectors for a 3 by 3 matrix is c(1,0,0);
#c(0,1,0); c(0,0,1).
#
#To project stage-biased dynamics, just
#don't specify a population vector,.
SBpr1<-project(Tort, time=50)
SBpr1
#
#Notice that this is now a matrix of 
#projections: each column is a different
#stage bias.
#
#We can also return population vectors as
#before:
SBpr1.2<-project(Tort, time=10, return.vec=TRUE)
SBpr1.2
#
#This takes the form of an array, with dimensions
#[time,stage,bias].
#
#Choose the density of stage 3 in bias 1 at the
#7th time interval:
SBpr1.2$vec[8,3,1]
#
#Choose the time series of stage 2 in bias 8:
SBpr1.2$vec[,2,8]
#
#Choose the population vector of bias 7 at the 
#5th time interval:
SBpr1.2$vec[6,,7]
#
#We can also standardise the matrix as before:
SBpr2<-project(Tort, time=50, standard.A=TRUE)
#
#And this is also easy to plot:
plot(SBpr2,log="y")
#
#
#
#END
