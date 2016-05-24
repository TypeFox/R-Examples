#=================================================
#DEMO INDEX:--------------------------------------
#
#1. Useful tools for working with matrices
#   1.1 USING MATLAB STYLE MATRIX NOTATION
#   1.2 PRIMITIVITY, REDUCIBILITY AND ERGODICITY
#   1.3 DEAING WITH REDUCIBLE MATRICES
#
#2. Population projection
#   2.1 PROJECTING SPECIFIC POPULATION STRUCTURES
#   2.2 PROJECTING STAGE-BIASED DYNAMICS
#
#3. Transient dynamics
#   3.1 TRANSIENT DYNAMICS OF SPECIFIC POPULATION 
#       STRUCTURES
#   3.2 BOUNDS ON TRANSIENT DYNAMICS
#
#4. Transfer function analyses
#   4.1 TRANSFER FUNCTION OF ASYMPTOTIC GROWTH
#   4.2 TRANSFER FUNCTION OF INERTIA
#   4.3 MATRIX TRANSFER FUNCTION PLOTS
#
#5. Convergence
#   5.1 DAMPING RATIO
#   5.2 DISTANCE MEASURES
#   5.3 SIMULATED TIME TO CONVERGENCE
#
#-------------------------------------------------
#=================================================
#
#Individual sections are available as seperate
#demos:
#1. demo(matrixtools)
#2. demo(projection)
#3. demo(transient)
#4. demo(transfer)
#5. demo(convergence)
#
readline("Hit <Enter> for next topic, or <Esc> to cancel")




#1.===============================================
#Useful tools for working with matrices
#References:
#Caswell 2001 Matrix population models
#(primitivity)
#Stott et al. 2010 Methods Ecol. Evol., 1, 242-252
#(reducibility and ergodicity)
#=================================================
#
#
#
#1.1----------------------------------------------
#USING MATLAB STYLE MATRIX NOTATION
#-------------------------------------------------
#
#Matlab style notation for matrices
#popdemo contains a function enabling Matlab
#style notation of matrices. As well as being
#easier to use, it crucially enables import of 
#many matrices simultaneously, using comma-
#seperated csv files imported as dataframes.
#
#This matrix is for the desert tortoise
#(Gopherus Agassizii) with medium fecundity
#(Doak et al. 1994 Ecol. Appl., 4, 446-480).
Tort<-Matlab2R("[0 0 0 0 0 1.3 1.98 2.57;0.716 0.567 0 0 0 0 0 0;0 0.149 0.567 0 0 0 0 0;0 0 0.149 0.604 0 0 0 0;0 0 0 0.235 0.56 0 0 0;0 0 0 0 0.225 0.678 0 0;0 0 0 0 0 0.249 0.851 0;0 0 0 0 0 0 0.016 0.86]")
Tort
#
readline("Hit <Enter> for next topic, or <Esc> to cancel")
#
#
#
#1.2----------------------------------------------
#PRIMITIVITY, REDUCIBILITY & ERGODICITY
#-------------------------------------------------
#
#popdemo provides a few tools for diagnosing the
#ergodic properties of the matrix.
#
#Is the matrix primitive?
is.matrix_primitive(Tort)
#
#Is the matrix reducible?
is.matrix_irreducible(Tort)
#
#Is the matrix ergodic?
is.matrix_ergodic(Tort)
#
readline("Hit <Enter> for next topic, or <Esc> to cancel")
#
#
#
#1.3----------------------------------------------
#DEALING WITH REDUCIBLE MATRICES
#-------------------------------------------------
#
#Create a reducible matrix
TortR<-Tort
TortR[7,6]<-0
is.matrix_irreducible(TortR)
#
#Block-permute the reducible matrix
blockmatrix(TortR)
TortRblock<-blockmatrix(TortR)$blockmatrix
eigen(TortR)$values
eigen(TortRblock[1:6,1:6])$values
eigen(TortRblock[7:8,7:8])$values
#
#The first diagonal block has the dominant
#eigenvalue and the second diagonal block has
#the first subdominant eigenvalue. Therefore
#the matrix should still be ergodic.
is.matrix_ergodic(TortR)
#
readline("Hit <Enter> for next topic, or <Esc> to cancel")





#2.===============================================
#Population projections
#=================================================
#
#
#
#2.1----------------------------------------------
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
readline("Hit <Enter> for next topic, or <Esc> to cancel")





#3.===============================================
#Transient dynamics
#References:
#Stott et al. 2011 Ecol. Lett., 14, 959-970
#(and citations therein)
#=================================================
#
#
#
#3.1----------------------------------------------
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
readline("Hit <Enter> for next topic, or <Esc> to cancel")





#4.===============================================
#Transfer function analyses
#Transfer function analysis is a means of exact
#perturbation analysis.  It looks at how a change
#in the vital rates of the population (the matrix
#elements or their underlying parameters)
#translates to a change in population dynamics.
#popdemo contains tools for analysing transfer
#functions of both asymptotic growth and transient
#dynamics.
#References:
#Hodgson & Townley (2004) J. Appl. Ecol. 41,
#1155-1161;
#Hodgson et al. (2006) Theor. Popul. Biol., 70,
#214-224.
#=================================================
#
#
#
#4.1----------------------------------------------
#TRANSFER FUNCTION OF ASYMPTOTIC GROWTH
#-------------------------------------------------
#
#Analysing a transfer function requires specifying
#the structure of the perturbation to the
#population.  This is done by multiplying a column
#vector and a row vector together.  The column
#vector, d, specifies the rows to be perturbed, and
#the row vector, e, specifies the columns to be
#perturbed.  The magnitude of the entries in d and
#e determines the relative magnitude of perturbation.
#
#Create a perturbation structure for the tortoise
#that targets only growth of stage 7:
d<-c(0,0,0,0,0,0,0,1)
e<-c(0,0,0,0,0,0,1,0)
d%*%t(e)
#
#Create a perturbation structure for the tortoise
#that equally targets all fecundity parameters:
d<-c(1,0,0,0,0,0,0,0)
e<-Tort[1,]
d%*%t(e)
#
#Create a perturbation structure that trades off 
#fecundity of the oldest individuals against their
#survival:
d<-c(-1,0,0,0,0,0,0,0.2)
e<-c(0,0,0,0,0,0,0,1)
d%*%t(e)
#
#...you get the idea.
#
#By specifying the PPM, the perturbation structure,
#and the range of perturbation magnitude, we can
#match perturbation to its resultant asymptotic growth.
#First, perturbing survival of old individuals:
tflam1<-tfa(Tort, d=c(0,0,0,0,0,0,0,1), e=c(0,0,0,0,0,0,1,0),
            prange=seq(-0.1,0.5,0.01))
tflam1
#
#We can easily plot this:
plot(tflam1)
#
#Perturbing fecundity:
tflam2<-tfa(Tort, d=c(1,0,0,0,0,0,0,0), e=Tort[1,],
            prange=seq(-1,2,0.1))
tflam1
#
#Plot:
plot(tflam2)
#
#Trade-off:
tflam3<-tfa(Tort, d=c(-1,0,0,0,0,0,0,0.2), e=c(0,0,0,0,0,0,0,1),
            prange=seq(-1,0.5,0.1))
tflam3
#
#Plot:
plot(tflam3)
#
#It is also possible to specify a desired lambda: for example,
#if we desire lambda=1 from the trade-off perturbation:
tfa(Tort, d=c(-1,0,0,0,0,0,0,0.2), e=c(0,0,0,0,0,0,0,1),
    lambdarange=1)
#
#This shows a perturbation magnitude of 0.6678 is required
#to achieve lambda=1.
#
readline("Hit <Enter> for next topic, or <Esc> to cancel")
#
#
#
#4.2----------------------------------------------
#THE TRANSFER FUNCTION OF INERTIA
#-------------------------------------------------
#
#The same logic applies to the transfer function
#of population inertia.  The transfer function
#requires a PPM, a perturbation structure, and a
#range of perturbation values.  Additionally,
#this formula requires either a population vector
#to work with, or a command to calculate the
#transfer function of the upper or lower bound.
#
#Transfer function that perturbs survival
#of the oldest individuals, using a specified
#population vector:
Tortvec<-Matlab2R("[1;1;2;3;5;8;13;21]")
tfin1<-inertia.tfa(Tort, vector=Tortvec, 
                   d=c(0,0,0,0,0,0,0,1), e=c(0,0,0,0,0,0,1,0),
                   prange=seq(-0.1,0.5,0.01))
tfin1
plot(tfin1)
#
#Perturb fecundity, analysing upper bound:
tfin2<-inertia.tfa(Tort, bound="upper",
                   d=c(1,0,0,0,0,0,0,0), e=Tort[1,],
                   prange=seq(-1,2,0.1))
tfin2
plot(tfin2)
#
#Perturb using trade-off, analyse lower bound:
tfin3<-inertia.tfa(Tort,bound="lower",
                   d=c(-1,0,0,0,0,0,0,0.2),e=c(0,0,0,0,0,0,0,1),
                   prange=seq(-1,0.5,0.1))
tfin3
plot(tfin3)
#
readline("Hit <Enter> for next topic, or <Esc> to cancel")
#
#
#
#4.3----------------------------------------------
#MATRIX TRANSFER FUNCTION PLOTS
#-------------------------------------------------
#
#When unsure of a perturbation structure to use,
#it can be helpful to look at how perturbing each
#element of the matrix individually is useful.
#popdemo contains functions to look at perturbation
#across the life cycle. 
#
tflammat<-tfamatrix(Tort)
plot(tflammat)
#
tfinmat1<-inertia.tfamatrix(Tort, vector=Tortvec)
plot(tfinmat1)
#
tfinmat2<-inertia.tfamatrix(Tort, bound="upper")
plot(tfinmat2)
#
tfinmat3<-inertia.tfamatrix(Tort, bound="lower")
plot(tfinmat3)
#
readline("Hit <Enter> for next topic, or <Esc> to cancel")





#5.===============================================
#Convergence
#References:
#Caswell 2001 Matrix population models
#(and citations therein)
#=================================================
#
#
#
#5.1----------------------------------------------
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
#5.2----------------------------------------------
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
#5.3----------------------------------------------
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




#=================================================
#-------------------------------------------------
#END
#-------------------------------------------------
#=================================================
