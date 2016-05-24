## =============================================================================
##
## This is a stiff system of 20 non-linear Ordinary Differential Equations.
##
## It describes a chemical reaction part of the air pollution model developed at
## The Dutch National Institute of Public Health and Environmental Protection
## (RIVM),
## and consists of 25 reaction and 20 reacting compounds.
##
## The reaction rates vary from e-3 to e+12, making the model extremely stiff
## =============================================================================

#
# A FORTRAN implementation (and reference output) can be found at 
# http://pitagora.dm.uniba.it//~testset
# F. Mazzia and F. Iavernaro. Test Set for Initial Value Problem Solvers. 
# Department of Mathematics, University of Bari, August 2003. 
# Available at http://www.dm.uniba.it/~testset.

# The model is described in Verwer (1994)
# J.G. Verwer, 1994. Gauss-Seidel iteration for stiff ODEs from chemical kinetics. 
# SIAM J. Sci. Comput., 15(5):1243-1259.

# 20 chemical species are described: NO2, NO, O3P, O3, HO2, OH,
# HCHO, CO, ALD, MEO2, C2O3, CO2, PAN, CH3O, HNO3, O1D, SO2, SO4, NO3, N2O5
 
# The model describes the following reactions:
#  r1:  NO2 -> NO + O3P               
#  r2:  NO + O3 -> NO2                
#  r3:  HO2+NO -> NO2                 
#  r4:  HCHO -> 2 HO2 + CO            
#  r5:  HCHO -> CO                    
#  r6:  HCHO + OH -> HO2+CO           
#  r7:  ALD + OH -> C2O3              
#  r8:  ALD -> MEO2+HO2+C)            
#  r9:  C2O3 + NO -> NO2 + MEO2 + CO2 
#  r10: C2O3 + NO2 -> PAN             
#  r11: PAN -> C2O3 + NO2             
#  r12: MEO2 + NO -> CH3O + NO2       
#  r13: CH3O -> HCHO + HO2            
#  r14: NO2+OH -> HNO3                
#  r15: O3P -> O3                     
#  r16: O3 -> O1D                     
#  r17: O3 -> O3P                     
#  r18: O1D -> 2 OH                   
#  r19: O1D -> O3P                    
#  r20: SO2 + Oh -> SO4 + HO2         
#  r21: NO3 -> NO                     
#  r22: NO3 -> NO2 + O3P              
#  r23: NO2 + O3 -> NO3               
#  r24: NO3 + NO2 -> N2O5             
#  r25: N2O5 -> NO3 + NO2             
      
require(deTestSet)

#=======================
# the model definition
#=======================
Pollution <- function (t, y, pars) {

 r  <- vector(length = 25)
 dy <- vector(length = length(y))
      r[ 1] <- k1 *y[ 1]
      r[ 2] <- k2 *y[ 2]*y[4]
      r[ 3] <- k3 *y[ 5]*y[2]
      r[ 4] <- k4 *y[ 7]
      r[ 5] <- k5 *y[ 7]
      r[ 6] <- k6 *y[ 7]*y[6]
      r[ 7] <- k7 *y[ 9]
      r[ 8] <- k8 *y[ 9]*y[6]
      r[ 9] <- k9 *y[11]*y[2]
      r[10] <- k10*y[11]*y[1]
      r[11] <- k11*y[13]
      r[12] <- k12*y[10]*y[2]
      r[13] <- k13*y[14]
      r[14] <- k14*y[ 1]*y[6]
      r[15] <- k15*y[ 3]
      r[16] <- k16*y[ 4]
      r[17] <- k17*y[ 4]
      r[18] <- k18*y[16]
      r[19] <- k19*y[16]
      r[20] <- k20*y[17]*y[6]
      r[21] <- k21*y[19]
      r[22] <- k22*y[19]
      r[23] <- k23*y[ 1]*y[4]
      r[24] <- k24*y[19]*y[1]
      r[25] <- k25*y[20]


      dy[1]  <-  -r[1] - r[10]- r[14]- r[23]- r[24]+ r[2] + r[3]+
                  r[9] + r[11]+ r[12]+ r[22]+ r[25]
      dy[2]  <-  -r[2] - r[3] - r[9] - r[12]+ r[1] + r[21]
      dy[3]  <-  -r[15]+ r[1] + r[17]+ r[19]+ r[22]
      dy[4]  <-  -r[2] - r[16]- r[17]- r[23]+ r[15]
      dy[5]  <-  -r[3] + r[4] + r[4] + r[6] + r[7] + r[13] + r[20]
      dy[6]  <-  -r[6] - r[8] - r[14]- r[20]+ r[3] + r[18] + r[18]
      dy[7]  <-  -r[4] - r[5] - r[6] + r[13]
      dy[8]  <-  +r[4] + r[5] + r[6] + r[7]
      dy[9]  <-  -r[7] - r[8]
      dy[10] <-  -r[12]+ r[7] + r[9]
      dy[11] <-  -r[9] - r[10]+ r[8] + r[11]
      dy[12] <-  +r[9]
      dy[13] <-  -r[11]+ r[10]
      dy[14] <-  -r[13]+ r[12]
      dy[15] <-  +r[14]
      dy[16] <-  -r[18]- r[19]+ r[16]
      dy[17] <-  -r[20]
      dy[18] <-  +r[20]
      dy[19] <-  -r[21]- r[22]- r[24] + r[23] + r[25]
      dy[20] <-  -r[25]+ r[24]
      
      return(list(c(dy = dy), rate = r))
}

#=============================
# parameters, state variables
#=============================
# The species names:
spnames <- c("NO2", "NO", "O3P", "O3", "HO2",
             "OH", "HCHO", "CO", "ALD", "MEO2",
             "C2O3", "CO2", "PAN", "CH3O", "HNO3",
             "O1D", "SO2", "SO4", "NO3", "N2O5")
# State variable initial condition
   y <- rep(0,20)
   y[2]  <- 0.2
   y[4]  <- 0.04
   y[7]  <- 0.1
   y[8]  <- 0.3
   y[9]  <- 0.01
   y[17] <- 0.007
   names (y) <- spnames

# Parameters: rate coefficients
   k1  =.35    ;  k2  =.266e2 ;  k3  =.123e5 ;  k4  =.86e-3
   k5  =.82e-3 ;  k6  =.15e5  ;  k7  =.13e-3 ;  k8  =.24e5
   k9  =.165e5 ;  k10 =.9e4   ;  k11 =.22e-1 ;  k12 =.12e5
   k13 =.188e1 ;  k14 =.163e5 ;  k15 =.48e7  ;  k16 =.35e-3
   k17 =.175e-1;  k18 =.1e9   ;  k19 =.444e12;  k20 =.124e4
   k21 =.21e1  ;  k22 =.578e1 ;  k23 =.474e-1;  k24 =.178e4
   k25 =.312e1


#=============================
# application 1.   
#=============================

times <- seq(0, 10, 0.1)
# run with default tolerances, short period of time
out   <- mebdfi(y = y, times = times, func = Pollution, parms = NULL)

# increasing tolerance
out2  <- mebdfi(y = y, times = times, func = Pollution, parms = NULL,
              atol = 1e-10, rtol = 1e-10)

out2b <- gamd(y = y, times = times, func = Pollution, parms = NULL,
             atol = 1e-10, rtol = 1e-10)

# run for longer period
Times <- seq (0, 2000, 10)
out3  <- gamd(y = y, times = Times, func = Pollution, parms = NULL,
             atol = 1e-10, rtol = 1e-10)

# plotting output
mf    <- par(mfrow = c(2,2))
plot (out, out2, which = "HO2", log = "y", ylab = "log", mfrow = NULL)
legend("topright",c("tol=1e-8","tol=1e-10"),col = 1:2, lty = 1)
plot (out, out2, which = "HCHO", ylab = "conc", mfrow = NULL)

plot (out3, which = "HO2", log = "y", ylab = "log", mfrow = NULL)
plot (out3, which = "HCHO", ylab = "conc", mfrow = NULL)

mtext(side = 3, outer = TRUE, line = -1.5,
      cex = 1.5, "Pollution problem")
par (mfrow = mf)

