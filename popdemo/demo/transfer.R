#=================================================
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
#This matrix is for the desert tortoise
#(Gopherus Agassizii) with medium fecundity
#(Doak et al. 1994 Ecol. Appl., 4, 446-480).
Tort<-Matlab2R("[0 0 0 0 0 1.3 1.98 2.57;0.716 0.567 0 0 0 0 0 0;0 0.149 0.567 0 0 0 0 0;0 0 0.149 0.604 0 0 0 0;0 0 0 0.235 0.56 0 0 0;0 0 0 0 0.225 0.678 0 0;0 0 0 0 0 0.249 0.851 0;0 0 0 0 0 0 0.016 0.86]")
#
#
#
#-------------------------------------------------
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
#that targets only the growth of stage 7:
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
#First, perturbing growth of stage 7:
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
tflam2
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
#-------------------------------------------------
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
#Transfer function that perturbs growth
#of stage 7, using a specified
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
#-------------------------------------------------
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
#
#
#END
