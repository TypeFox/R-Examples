
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General 
# Public License along with this library; if not, write to the 
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
# MA  02111-1307  USA

# Copyrights (C)
# for this R-port: 
#   1999 - 2004, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:                  DESCRIPTION:
#  BesselI                    Modified Bessel Function of first kind
#  BesselK                    Modified Bessel Function of third kind
#  BesselDI                   Derivative of BesselI
#  BesselDK                   Derivative of BesselK
# INTERNAL FUNCTION:         DESCRIPTION:
#  .BesselN                   For internal use only 
#  .Bessel01                   ...
#  .Bessel.MSTA1               ...
#  .Bessel.MSTA2               ...
#  .Bessel.ENVJ                ...
################################################################################
 
 
BesselI = 
function(x, nu, expon.scaled = FALSE) 
{   # A function implemented by Diethelm Wuertz

    # Symmetry Relation:
    nu = abs(nu)
    
    # Test:
    if (nu - floor(nu) != 0) stop("nu must be an integer")
    
    # Compute:
    bessel = NULL
    for (X in x) {
        bessel = c(bessel, .BesselN(X, nu)[1])
    }
    
    # Scaling:
    if (expon.scaled) bessel = exp(-x)*bessel
    
    # Return Value:
    as.vector(bessel)
}


# ------------------------------------------------------------------------------


BesselK = 
function(x, nu, expon.scaled = FALSE) 
{   # A function implemented by Diethelm Wuertz

    # Symmetry Relation:
    nu = abs(nu)
    
    # Test:
    if (nu - floor(nu) != 0) stop("nu must be an integer")
    
    # Compute:
    bessel = NULL
    for (X in x) {
        bessel = c(bessel, .BesselN(X, nu)[2])
    }
    
    # Scaling:
    if (expon.scaled) bessel = exp(x)*bessel
    
    # Return Value:
    as.vector(bessel)  
}


# ------------------------------------------------------------------------------


BesselDI = 
function(x, nu) 
{   # A function implemented by Diethelm Wuertz

    # Symmetry Relation:
    nu = abs(nu)
    
    # Test:
    if (nu - floor(nu) != 0) stop("nu must be an integer")
    
    # Compute:
    bessel = NULL
    for (X in x) {
        bessel = c(bessel, .BesselN(X, nu)[3])
    }
    
    # Return Value:
    bessel  
}


# ------------------------------------------------------------------------------


BesselDK = 
function(x, nu) 
{   # A function implemented by Diethelm Wuertz

    # Symmetry Relation:
    nu = abs(nu)
    
    # Test:
    if (nu - floor(nu) != 0) stop("nu must be an integer")
        
    # Compute:
    bessel = NULL
    for (X in x) {
        bessel = c(bessel, .BesselN(X, nu)[4])
    }
    
    # Return Value:
    bessel
}


################################################################################
 
            
.BesselN = 
function(X, N) 
{   # A function implemented by Diethelm Wuertz
 
    # Description: 
    #   Compute modified Bessel functions In(x) and 
    #   Kn(x), and their derivatives  
    
    # Arguments:
    #   x - argument of In(x) and Kn(x)           
    #   n - order of In(x) and Kn(x)   
    
    # Value:           
    #   BI(n) --- In(x)               
    #   DI(n) --- In'(x)              
    #   BK(n) --- Kn(x)               
    #   DK(n) --- Kn'(x)              
    #   NM --- Highest order computed 

    # FUNCTION:          
                  
    # Settings:
    F = 0      
    NM = N  
        
    # Very small arguments:
    if (X == 0) {  
        if (N == 0) {
           return(c(BI = 1, BK = Inf, DI = 0, DK = -Inf)) 
        }     
       if (N >= 1) {
           return(c(BI = 0, BK = Inf, DI = 0, DK = -Inf))    
        }
    }
         
    # Start values:
    BI = BK = DI = DK = rep(NA, times = max(N+1, 2))
    Bessel = .Bessel01(X)        
    BI0 = BI[1] = Bessel["BI0"]           
    BI1 = BI[2] = Bessel["BI1"]           
    BK0 = BK[1] = Bessel["BK0"]           
    BK1 = BK[2] = Bessel["BK1"]           
    DI0 = DI[1] = Bessel["DI0"]           
    DI1 = DI[2] = Bessel["DI1"]           
    DK0 = DK[1] = Bessel["DK0"]           
    DK1 = DK[2] = Bessel["DK1"]   
           
    # Return for N=0:
    if (N == 0) return(c(BI = BI0, BK = BK0, DI = DI0, DK = DK0))
    
    # Return for N=1:
    if (N == 1) return(c(BI = BI1, BK = BK1, DI = DI1, DK = DK1)) 
        
    # Compute BI for N>1:
    if (X > 40 & N < floor(0.25*X)) {            
        H0 = BI0          
        H1 = BI1        
        for (K in 2:N) {     
            H = -2 * (K-1) / X*H1 + H0          
            BI[K+1] = H          
            H0 = H1            
            H1 = H    
        }      
    } else { 
        M = .Bessel.MSTA1(X, 200)   
        if (M < N) { 
              NM = M          
        } else {          
              M = .Bessel.MSTA2(X, N, 15)              
        }           
        F0 = 0        
        F1 = 1.0e-100      
        for (I in 0:M) {       
            # (K = M, 0, -1) 
            K = M - I
            F = 2 * (K+1) * F1/X + F0        
            if (K <= NM) BI[K+1] = F             
            F0 = F1         
            F1 = F
        }          
        S0 = BI0/F         
        for (K in 0:NM) BI[K+1] = S0 * BI[K+1]                 
    } 
    for (K in 2:NM) {      
       DI[K+1] = BI[K-1+1] - K/X * BI[K+1]             
    }
       
    # Compute BK for N>1:      
    G0 = BK0             
    G1 = BK1              
    for (K in 2:NM) {        
       G = 2 * (K-1) / X*G1 + G0           
       BK[K+1] = G          
       G0 = G1            
       G1 = G  
    }       
    for (K in 2:NM) {                 
       DK[K+1] = -BK[K-1+1] - K/X * BK[K+1]  
    }
    
    # Result:
    ans = c(BI = BI[N+1], BK = BK[N+1], DI = DI[N+1], DK = DK[N+1])
    names(ans) = NULL
    
    # Return Value:
    ans
          
}   

                     
# ------------------------------------------------------------------------------


.Bessel01 =
function(X) 
{   # A function implemented by Diethelm Wuertz
                 
    # Description: 
    #   Compute modified Bessel functions I0(x), I1(1),                
    #   K0(x) and K1(x), and their derivatives  
    
    # Arguments:    
    #   x - argument   
    
    # Values:           
    #   BI0 --- I0(x)                 
    #   DI0 --- I0'(x)                
    #   BI1 --- I1(x)                 
    #   DI1 --- I1'(x)                
    #   BK0 --- K0(x)                 
    #   DK0 --- K0'(x)                
    #   BK1 --- K1(x)                 
    #   DK1 --- K1'(x)  
    
    # FUNCTION:
                                                  
    # Compute BI and BK:
    if (X == 0) {                   
        BI0 = 1 
        BI1 = 0
        BK0 = Inf
        BK1 = Inf  
        DI0 = 0
        DI1 = 0.5
        DK0 = -Inf
        DK1 = -Inf  
        return(c(
            BI0 = BI0, BI1 = BI1, BK0 = BK0, BK1 = BK1, 
            DI0 = DI0, DI1 = DI1, DK0 = DK0, DK1 = DK1))
    }
    
    # Compute BI:
    if (X <= 18) {                                    
        bi0.fun = function(X) {
            X2 = X * X
            BI0 = 1 
            R = 1          
            for (K  in 1:50) {     
                R = 0.25 * R * X2 / (K*K)              
                BI0 = BI0 + R     
                if (abs(R/BI0) < 1.0e-15) return(BI0)           
            }
            BI0
        } 
        BI0 = bi0.fun(X) 
        bi1.fun = function(X) {
            X2 = X * X
            BI1 = 1        
            R = 1          
            for (K in 1:50) {      
                R = 0.25 * R * X2 /(K*(K+1))          
                BI1 = BI1 + R     
                if (abs(R/BI1) < 1.0e-15) return(0.5 * X * BI1)          
            }         
            0.5 * X * BI1 
        } 
        BI1 = bi1.fun(X) 
    } else {
       A = c(
             0.125,7.03125e-2,     7.32421875e-2,         1.1215209960938e-1,          
             2.2710800170898e-1,   5.7250142097473e-1,    1.7277275025845, 
             6.0740420012735,     24.380529699556,      110.01714026925,     
             551.33589612202,      3.0380905109224e03 ) 
       B = c(
            -0.375, -1.171875e-1, -1.025390625e-1,       -1.4419555664063e-1,       
            -2.7757644653320e-1,  -6.7659258842468e-1,   -1.9935317337513, 
            -6.8839142681099e0,   -2.7248827311269e01, -121.59789187654,   
            -6.0384407670507e02,  -3.3022722944809e03 )  
       K0 = 12            
       if (X >= 35) K0 = 9                 
       if (X >= 50) K0 = 7                 
       CA = exp(X) / sqrt(2 * pi * X)           
       BI0 = 1        
       XR = 1/X       
       for (K in 1:K0) BI0 = BI0 + A[K] * XR^K       
       BI0 = CA * BI0       
       BI1 = 1     
       for (K in 1:K0) BI1 = BI1 + B[K] * XR^K            
       BI1 = CA * BI1        
    }  
       
    # Compute BK:
    if (X <= 9) {                   
        bk0.fun = function(X) {
            X2 = X * X
            EL = 0.5772156649015329
            CT = -(log(X/2) + EL)              
            BK0 = 0        
            WW = BK0         
            W0 = 0         
            R = 1          
            for (K in 1:50) {     
                W0 = W0 + 1/K 
                R = 0.25 * R / (K*K) * X2           
                BK0 = BK0 + R * (W0 + CT)                
                if (abs(BK0-WW)/abs(BK0) < 1.0e-15 & abs(BK0) > 0) 
                    return(BK0 + CT)
                WW = BK0 
            }      
            BK0 + CT 
        }
        BK0 = bk0.fun(X)
    } else {  
        A1 = c(
            0.125, 0.2109375, 1.0986328125, 11.775970458984,        
            214.61706161499, 5.9511522710323e03, 2.3347645606175e05, 
            1.2312234987631e07 )    
        CB = 0.5 / X       
        XR2 = 1 / (X*X)    
        BK0 = 1        
        for (K in 1:8) BK0 = BK0 + A1[K] * XR2^K             
        BK0 = CB * BK0/BI0  
    } 
    BK1 = (1/X - BI1*BK0)/BI0    
        
    # Derivatives:          
    DI0 = BI1             
    DI1 = BI0 - BI1/X       
    DK0 = -BK1            
    DK1 = -BK0 - BK1/X                   
    
    # Return Value:
    c(
        BI0 = BI0, BI1 = BI1, BK0 = BK0, BK1 = BK1, 
        DI0 = DI0, DI1 = DI1, DK0 = DK0, DK1 = DK1)
}     

        
# ------------------------------------------------------------------------------


.Bessel.MSTA1 =
function(X, MP) 
{   # A function implemented by Diethelm Wuertz
               
    # Description: 
    #   Determine the starting point for backward recurrence such  
    #   that the magnitude of Jn(x) at that point is about 10^(-MP)      
    
    # Arguments:       
    #   x - argument of Jn(x)              
    #   MP - Value of Magnitude      

    # Value:
    #   MSTA1 - Starting point      
          
    # FUNCTION:
     
    # Settings:
    A0 = abs(X)          
    N0 = floor(1.1 * A0) + 1    
    F0 = .Bessel.ENVJ(N0, A0) - MP   
    N1 = N0 + 5             
    F1 = .Bessel.ENVJ(N1, A0) - MP   
        
    # Compute:
    for (IT in 1:20) {       
        NN = N1 - floor( (N1-N0) / (1 - F0/F1) )         
        F = .Bessel.ENVJ(NN, A0) - MP 
        if (abs(NN-N1) < 1) return(NN)        
        N0 = N1            
        F0 = F1            
        N1 = NN            
        F1 = F   
    }
    
    # Return Value:
    NN
}             


# ------------------------------------------------------------------------------


.Bessel.MSTA2 =
function(X, N, MP) 
{   # A function implemented by Diethelm Wuertz
                
    # Description: 
    #   Determine the starting point for backward recurrence such  
    #   that all Jn(x) has MP significant digits     
    
    # Arguments:       
    #   x - argument of Jn(x)      
    #   n - Order of Jn(x)         
    #   MP - Significant digit      

    # Value:
    #   MSTA2 - Starting point      
          
    # FUNCTION:
     
    # Settings:
    A0 = abs(X)         
    HMP = 0.5 * MP        
    EJN = .Bessel.ENVJ(N, A0)  
        
    if (EJN <= HMP) {                   
       OBJ = MP           
       N0 = max(floor(1.1 * A0), 1)   
    } else {  
       OBJ = HMP + EJN      
       N0 = N             
    }
     
    # Compute:
    F0 = .Bessel.ENVJ(N0, A0) - OBJ  
    N1 = N0 + 5             
    F1 = .Bessel.ENVJ(N1, A0) - OBJ
    for (IT in 1:20) {             
       NN = N1 - floor( (N1-N0) / (1 - F0/F1) )
       # print(c(N0=N0, A0=A0))  
       # print(c(F0=F0, F1=F1))
       # print(c(N=N, NN=NN, N1=N1))         
       F = .Bessel.ENVJ(NN, A0) - OBJ                 
       if (abs(NN-N1) < 1) return(NN + 10)       
       N0 = N1            
       F0 = F1            
       N1 = NN            
       F1 = F  
    }   
    
    # Return Value:
    NN + 10
}         
  
          
# ------------------------------------------------------------------------------
       
     
.Bessel.ENVJ =
function(N, X)                
{     # A function implemented by Diethelm Wuertz

      # Return Value:
      0.5 * log10(6.28*N) - N * log10(1.36*X/N)       
}             
              
  
################################################################################

              