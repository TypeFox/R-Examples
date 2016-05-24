#=================================================
#Convergence
#References:
#Caswell 2001 Matrix population models
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
#DAMPING RATIO
#-------------------------------------------------
#
#The damping ratio is the ratio of the dominant
#eigenvalue to the magnitude of the subdominant
#eigenvalue.  The larger the damping ratio, the
#quicker the model converges.
dr(Tort)
#
#Dividing the log of the damping ratio by the log
#of x gives the time (in projection intervals)
#for the the influence of the dominant eigenvalue 
#in the population projection to become
#x times larger than the subdominant.  This
#provides an estimation of time to convergence.
#
dr(Tort, x=10, return.t=TRUE)
#
#It takes 17.35 projection intervals for the
#influence of the dominant eigenvalue to become 
#10 times larger than that of the subdominant 
#eigenvalue.
#
readline("Hit <Enter> for next topic, or <Esc> to cancel")
#
#
#
#-------------------------------------------------
#DISTANCE MEASURES
#-------------------------------------------------
#
#Distance measures are a means of measuring how
#'far' from stable stage structure the population
#vector is.
#
#Keyfitz's Delta literally measures the mathematical
#distance between the two vectors in n-dimensional
#space.
Tortvec<-Matlab2R("[1;1;2;3;5;8;13;21]")
Keyfitz.delta(Tort, vector=Tortvec)
#
#Projection distance incorporates the reproductive
#value of the population:
projection.distance(Tort, vector=Tortvec)
#
#Cohen's cumulative distance desribes the length
#of the 'path' the population takes in its projection
#until stability.
Cohen.cumulative(Tort, vector=Tortvec)
#
readline("Hit <Enter> for next topic, or <Esc> to cancel")
#
#
#
#-------------------------------------------------
#SIMULATED TIME TO CONVERGENCE
#-------------------------------------------------
#
#It may be useful to know the exact time to
#convergence. popdemo contains a function to
#simulate the projection and calculate time to
#a specified accuracy of convergence.
convergence.time(Tort, vector=Tortvec,
                 accuracy=1e-2)
#
#It takes 16 years for the population to converge
#within a proportion of 0.01 (1%) of lambda-max.
#
#We can also calculate the convergence time of
#stage-biased projections if we don't specify a
#population vector.  In practise, the 
#stage-biased vectors are standard basis vectors: 
#column vectors consisting of zeroes except for a 
#1 in one of the elements.  For example, a set of 
#stage-biased vectors for a 3 by 3 matrix is 
#c(1,0,0); c(0,1,0); c(0,0,1).
convergence.time(Tort, accuracy=1e-3)
#
#These are the times it takes for each stage-biased
#projection to converge to 0.1% of lambda-max.
#
#
#
#END
