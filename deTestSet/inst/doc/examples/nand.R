## =============================================================================
##
##     This file is adapted from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##
##        NAND gate
##        index 0 IDE of dimension 14
##     This is revision
##     $Id: nand.F,v 1.2 2006/10/02 10:29:14 testset Exp $
##
## =============================================================================

require(deTestSet)
      
#-----------------------------------------------------------------------
#
# The network equation describing the nand gate
#             C[Y] * Y' - f[Y,t] = 0
# 
# ---------------------------------------------------------------------

Nand <- function(t,       # time point t
                 Y,       # node potentials at time point t
                 Yprime,  # rate of change of Y
                 pars) {
#-----------------------------------------------------------------------
# Voltage-dependent capacitance matrix C(Y) for the network equation
#             C(Y) * Y' - f(Y,t) = 0
#-----------------------------------------------------------------------

      CAP[1,1]   <- CGS
      CAP[1,5]   <- -CGS
      CAP[2,2]   <- CGD
      CAP[2,5]   <- -CGD
      CAP[3,3]   <- CBDBS(Y[3]-Y[5])
      CAP[3,5]   <- -CBDBS(Y[3]-Y[5])
      CAP[4,4]   <- CBDBS(Y[4]-VDD)
      CAP[5,1]   <- -CGS
      CAP[5,2]   <- -CGD
      CAP[5,3]   <- -CBDBS(Y[3]-Y[5])
      CAP[5,5]   <- CGS+CGD-CAP[5,3]+  CBDBS(Y[9]-Y[5])+C9
      CAP[5,9]   <- -CBDBS(Y[9]-Y[5])
      CAP[6,6]   <- CGS
      CAP[7,7]   <- CGD
      CAP[8,8]   <- CBDBS(Y[8]-Y[10])
      CAP[8,10]  <- -CBDBS(Y[8]-Y[10])
      CAP[9,5]   <- -CBDBS(Y[9]-Y[5])
      CAP[9,9]   <- CBDBS(Y[9]-Y[5])
      CAP[10,8]  <- -CBDBS(Y[8]-Y[10])
      CAP[10,10] <- -CAP[8,10]+CBDBS(Y[14]-Y[10])+C9
      CAP[10,14] <- -CBDBS(Y[14]-Y[10])
      CAP[11,11] <- CGS
      CAP[12,12] <- CGD
      CAP[13,13] <- CBDBS(Y[13])
      CAP[14,10] <- -CBDBS(Y[14]-Y[10])
      CAP[14,14] <- CBDBS(Y[14]-Y[10])

# ---------------------------------------------------------------------
#          PULSE: Input signal in pulse form
# ---------------------------------------------------------------------
      P1  <- PULSE(t,0.0,5.0,5.0,5.0,5.0,5.0,20.0) 
      V1  <- P1$VIN
      V1D <- P1$VIND 

      P2  <- PULSE(t,0.0,5.0,15.0,5.0,15.0,5.0,40.0)
      V2  <- P2$VIN
      V2D <- P2$VIND 

#-----------------------------------------------------------------------
# Right-hand side f[X,t] for the network equation
#             C[Y] * Y' - f[Y,t] = 0
# External reference:
#          IDS: Drain-source current
#          IBS: Nonlinear current characteristic for diode between
#               bulk and source
#          IBD: Nonlinear current characteristic for diode between
#               bulk and drain
#-----------------------------------------------------------------------

      F[1]  <- -(Y[1]-Y[5])/RGS-IDS(1,Y[2]-Y[1],Y[5]-Y[1],Y[3]-Y[5],
                Y[5]-Y[2],Y[4]-VDD)
      F[2]  <- -(Y[2]-VDD)/RGD+IDS(1,Y[2]-Y[1],Y[5]-Y[1],Y[3]-Y[5],
                Y[5]-Y[2],Y[4]-VDD)
      F[3]  <- -(Y[3]-VBB)/RBS + IBS(Y[3]-Y[5])
      F[4]  <- -(Y[4]-VBB)/RBD + IBD(Y[4]-VDD)
      F[5]  <- -(Y[5]-Y[1])/RGS-IBS(Y[3]-Y[5])-(Y[5]-Y[7])/RGD-
                IBD(Y[9]-Y[5])
      F[6]  <- CGS*V1D-(Y[6]-Y[10])/RGS-
                IDS(2,Y[7]-Y[6],V1-Y[6],Y[8]-Y[10],V1-Y[7],Y[9]-Y[5])
      F[7]  <- CGD*V1D-(Y[7]-Y[5])/RGD+
                IDS(2,Y[7]-Y[6],V1-Y[6],Y[8]-Y[10],V1-Y[7],Y[9]-Y[5])
      F[8]  <- -(Y[8]-VBB)/RBS + IBS(Y[8]-Y[10])
      F[9]  <- -(Y[9]-VBB)/RBD + IBD(Y[9]-Y[5])
      F[10] <- -(Y[10]-Y[6])/RGS-IBS(Y[8]-Y[10])-
               (Y[10]-Y[12])/RGD-IBD(Y[14]-Y[10])
      F[11] <- CGS*V2D-Y[11]/RGS-IDS(2,Y[12]-Y[11],V2-Y[11],Y[13],
                V2-Y[12],Y[14]-Y[10])
      F[12] <- CGD*V2D-(Y[12]-Y[10])/RGD+
                IDS(2,Y[12]-Y[11],V2-Y[11],Y[13],V2-Y[12],Y[14]-Y[10])
      F[13] <- -(Y[13]-VBB)/RBS + IBS(Y[13])
      F[14] <- -(Y[14]-VBB)/RBD + IBD(Y[14]-Y[10])

#             C[Y] * Y' - f[Y,t] = 0
#      Delta <- colSums(t(CAP) * Yprime) - F
      Delta <- CAP %*% Yprime - F
      return(list(Delta, pulse1 = P1$VIN, pulse2 = P2$VIN))
     }

# ---------------------------------------------------------------------------
#
# Function evaluating the drain-current due to the model of
# Shichman and Hodges
#
# ---------------------------------------------------------------------------

IDS <- function (NED,   #   NED  Integer parameter for MOSFET-type
                 VDS,   #   VDS  Voltage between drain and source
                 VGS,   #   VGS  Voltage between gate and source
                 VBS,   #   VBS  Voltage between bulk and source
                 VGD,   #   VGD  Voltage between gate and drain
                 VBD)   #   VBD  Voltage between bulk and drain

{
    if ( VDS == 0 ) return(0)
     
    if (NED == 1) { #--- Depletion-type
      VT0    =-2.43 
      CGAMMA = 0.2 
      PHI    = 1.28 
      BETA   = 5.35e-4
    } else       { # --- Enhancement-type
      VT0    = 0.2
      CGAMMA = 0.035
      PHI    = 1.01
      BETA   = 1.748e-3
    }

    if ( VDS >  0 )  # drain function for VDS>0
    {
      SQRT1 <- ifelse (PHI-VBS>0, sqrt(PHI-VBS), 0)
      VTE   <- VT0 + CGAMMA * ( SQRT1 - sqrt(PHI) )

      if ( VGS-VTE <= 0.0) IDS <- 0.  else
       if ( 0.0 < VGS-VTE & VGS-VTE <= VDS )
         IDS <- - BETA * (VGS - VTE)^ 2.0 * (1.0 + DELTA*VDS) else
      if ( 0.0 < VDS & VDS < VGS-VTE ) 
       IDS <- - BETA * VDS * (2 *(VGS - VTE) - VDS) * (1.0 + DELTA*VDS)

      }  else    {

      SQRT2 <- ifelse (PHI-VBD>0, sqrt(PHI-VBD), 0)
      VTE   <- VT0 + CGAMMA * (SQRT2 - sqrt(PHI) )

      if ( VGD-VTE <= 0.0) IDS <- 0.0  else
      if ( 0.0 < VGD-VTE & VGD-VTE <= -VDS ) 
       IDS <- BETA * (VGD - VTE)^2.0 * (1.0 - DELTA*VDS) else
      if ( 0.0 < -VDS & -VDS < VGD-VTE ) 
       IDS <- - BETA * VDS * (2 *(VGD - VTE) + VDS) *(1.0 - DELTA*VDS)
      }
    return(IDS)
}

# ---------------------------------------------------------------------------
#
# Function evaluating the current of the pn-junction between bulk and
# source due to the model of Shichman and Hodges
#
# ---------------------------------------------------------------------------

IBS <- function (VBS)     #   VBS  Voltage between bulk and source
       ifelse ( VBS <= 0.0, - CURIS * ( exp( VBS/VTH ) - 1.0) , 0.0)

# ---------------------------------------------------------------------------
#
# Function evaluating the current of the pn-junction between bulk and
# drain  due to the model of Shichman and Hodges
#
# ---------------------------------------------------------------------------
IBD <- function (VBD)     #   VBD  Voltage between bulk and drain
       ifelse ( VBD <= 0.0, - CURIS * ( exp( VBD/VTH ) - 1.0) , 0.0)

# ---------------------------------------------------------------------------
#
# Evaluating input signal at time point X
#
# ---------------------------------------------------------------------------

PULSE <- function (X,     # Time-point at which input signal is evaluated
                   LOW,   # Low-level of input signal
                   HIGH,  # High-level of input signal
                   DELAY, T1, T2, T3, PERIOD)  # to specify signal structure
# ---------------------------------------------------------------------------
# Structure of input signal:
#
#                -----------------------                       HIGH
#               /                       \
#              /                         \
#             /                           \
#            /                             \
#           /                               \
#          /                                 \
#         /                                   \
#        /                                     \
#  ------                                       ---------      LOW
#
# |DELAY|   T1  |         T2           |   T3  |
# |          P     E     R     I     O     D            |
#
# ---------------------------------------------------------------------------
{
      TIME <- X%%PERIOD
      VIN  <- LOW
      VIND <- 0.0

      if (TIME > (DELAY+T1+T2)) 
      {
        VIN  <- ((HIGH-LOW)/T3)*(DELAY+T1+T2+T3-TIME) + LOW
        VIND <- -((HIGH-LOW)/T3) } else
      if (TIME > (DELAY+T1))
      {
        VIN  <- HIGH
        VIND <- 0.0              } else
      if (TIME > DELAY) 
      {
        VIN  <- ((HIGH-LOW)/T1)*(TIME-DELAY) + LOW
        VIND <- ((HIGH-LOW)/T1)  }

      return  (list(VIN = VIN,  # Voltage of input signal at time point X
                    VIND = VIND))  # Derivative of VIN at time point X

}

# ---------------------------------------------------------------------------
#
# Function evaluating the voltage-dependent capacitance between bulk and
# drain gevalp. source  due to the model of Shichman and Hodges
#
# ---------------------------------------------------------------------------

CBDBS <- function (V)   # Voltage between bulk and drain gevalp. source
         ifelse ( V <= 0.0, CBD/sqrt(1.0-V/0.87), CBD*(1.0+V/(2.0*0.87)))


#-----------------------------------------------------------------------
# solution
# computed at Cray C90, using Cray double precision:
# Solving NAND gate using PSIDE
#
# User input:
#
# give relative error tolerance: 1d-16
# give absolute error tolerance: 1d-16
#
#
# Integration characteristics:
#
#    number of integration steps       22083
#    number of accepted steps          21506
#    number of f evaluations          308562
#    number of Jacobian evaluations      337
#    number of LU decompositions       10532
#
# CPU-time used:                         451.71 sec
#
#      y[  1] =  0.4971088699385777d+1
#      y[  2] =  0.4999752103929311d+1
#      y[  3] = -0.2499998781491227d+1
#      y[  4] = -0.2499999999999975d+1
#      y[  5] =  0.4970837023296724d+1
#      y[  6] = -0.2091214032073855d+0
#      y[  7] =  0.4970593243278363d+1
#      y[  8] = -0.2500077409198803d+1
#      y[  9] = -0.2499998781491227d+1
#      y[ 10] = -0.2090289583878100d+0
#      y[ 11] = -0.2399999999966269d-3
#      y[ 12] = -0.2091214032073855d+0
#      y[ 13] = -0.2499999999999991d+1
#      y[ 14] = -0.2500077409198803d+1
#-----------------------------------------------------------------------

RGS     <-  4
RGD     <-  4
RBS     <-  10
RBD     <-  10
CGS     <-  0.6e-4
CGD     <-  0.6e-4
CBD     <-  2.4e-5
CBS     <-  2.4e-5
C9      <-  0.5e-4
DELTA   <-  0.2e-1
CURIS   <-  1.e-14
VTH     <-  25.85
VDD     <-  5.
VBB     <-  -2.5

#-----------------------------------------------------------------------
# initialising
VBB    <- -2.5
Y      <- c(5, 5, VBB, VBB, 5, 3.62385, 5,
            VBB, VBB, 3.62385, 0, 3.62385, VBB, VBB)
Yprime <- rep(0, 14)

#-----------------------------------------------------------------------
# memory allocation
CAP   <- matrix(nrow = 14, ncol = 14, data = 0)
F     <- vector("double", 14)
times <- seq(from = 0, to = 80, by = 1) # hour

# integrate the model: low tolerances to restrict integration time
out <- mebdfi(y = Y, dy = Yprime, times, res = Nand, parms = 0,
              hini = 1e-6, rtol = 1e-6, atol = 1e-6)

# plot output

par(mar = c(4, 2, 3, 2))
plot(out, mfrow = c(4, 4), lwd = 2)

#reference solution
ref<-c(4.971088699385777,    4.999752103929311,     -2.499998781491227,
       -2.499999999999975,   4.970837023296724,     -0.2091214032073855,
       4.970593243278363,   -2.500077409198803,     -2.499998781491227,
       -0.2090289583878100, -0.2399999999966269e-3, -0.2091214032073855,
       -2.499999999999991,  -2.500077409198803)

print (t(rbind(mebdfi = out [nrow(out),2:15] ,
               reference = ref,
               delt = out [nrow(out),2:15] - ref)))
