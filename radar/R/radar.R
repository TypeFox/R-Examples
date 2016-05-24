kConstantSLP = 1013.25 # Sea-level Pressure [hPa]
kConstantP0 = 1000.0  # Reference pressure [hPa]
kConstantSpeedOfLight = 3e8 # Speed of light [m/s]
kConstantRe = 6374000.0 # Earth's radius [m]
kConstantR43 = kConstantRe*4.0/3.0 # 4/3 Approximation effective radius for standard atmosphere [m]
kConstantBoltz = 1.381e-23 # Boltzmann's constant [ m^2 kg s^-2 K^-1]


AttenuationAbsCoeff <- function(D, lam, m){
# Absorption coefficient of a spherical particle
#     From Doviak and Zrnic (1993), Eqn 3.14a or Battan (1973), Eqn 6.6
#     INPUT::
#     -----
#     D : float
#         Particle diameter [m]
#     lam : float
#         Radar wavelength [m]
#     m : float
#         Complex refractive index [unitless]
#     OUTPUT::
#     ------
#     Qa :
#         Absorption coefficient [unitless]
    Km <- (m^2 - 1) / (m^2 + 2)
    (pi^2 * D^3 / lam) * Im(-1 * Km)
    }

AttenuationScatCoeff<- function(D, lam, m){
#     Scattering coefficient of a spherical particle.
#     From Doviak and Zrnic (1993), Eqn 3.14b or Battan (1973), Eqn 6.5
#     INPUT::
#     -----
#     D : float
#         Particle diameter [m]
#     lam : float
#         Radar wavelength [m]
#     m : float
#         Complex refractive index [unitless]
#     OUTPUT::
#     ------
#     Qs : float
#         Scattering coefficient [unitless]
Km = (m^2 - 1) / (m^2 + 2)
(2 * pi^5 * D^6 / (3 * lam^4) * (abs(Km))^2)
}

AttenuationExtCoeff<- function(D,lam,m){
#     Extinction coefficient of a spherical particle.
#     From Doviak and Zrnic (1993), Eqn 3.14b or Battan (1973), Eqn 6.5
#     INPUT::
#     D : float
#         Particle diameter [m]
#     lam : float
#         Radar wavelength [m]
#     m : float
#         Complex refractive index [unitless]
#     OUTPUT::
#     ------
#     Qe : float
#         Extinction coefficient [unitless]
Qa <- AttenuationAbsCoeff(D,lam,m)
Qs <- AttenuationScatCoeff(D,lam,m)
Qa + Qs
}

ConversiondBZ2Z<- function(dBZ){
#     Conversion from dBZ (log) units to linear Z units
#     INPUT::
#     -----
#     dBZ : float
#         logarithmic reflectivity value
#     OUTPUT::
#     ------
#     Zlin : float
#         linear reflectivity units
10^(dBZ/10)
    }
    
ConversionZ2dBZ<- function(Zlin){
#     Conversion from linear Z units to dBZ (log) units
#     INPUT::
#     -----
#     Zlin : float
#         linear reflectivity units
#     OUTPUT::
#     ------
#     dBZ : float
#         logarithmic reflectivity value
10 * log10(Zlin)
    }

DopplerFreq<- function(lam, speedOfLight){
#     Frequency given wavelength.
#     INPUT::
#     -----
#     lam : float
#         Wavelength [m]
#     OUTPUT::
#     ------
#     freq : float
#         Frequency [Hz]
speedOfLight/lam
    }
    
DopplerWavelength<- function(freq,speedOfLight){
#     Wavelength given frequency.
#     INPUT::
#     -----
#     freq : float
#         Frequency [Hz]
#     OUTPUT::
#     ------
#     lam : float
#         Wavelength [m]
speedOfLight/freq
}

DopplerPulseDuration<- function(tau, speedOfLight){
#     Pulse duration from pulse length.
#     INPUT::
#     -----
#     tau : float
#         Pulse length [m]
#     OUTPUT::
#     ------
#     pDur : float
#         Pulse duration [s]
2 * tau/speedOfLight
}

DopplerPulseLength<- function(pDur, speedOfLight){
#     Pulse length from pulse duration.
#     INPUT::
#     -----
#     pDur : float
#         Pulse duration [s]
#     OUTPUT::
#     ------
#     tau : float
#         Pulse length [m]
speedOfLight * pDur/2
}

DopplerFmax<- function(PRF){
#     Maximum frequency given PRF.
#     From Rinehart (1997), Eqn 6.8
#     INPUT::
#     -----
#     PRF : float
#         Pulse repetition frequency [Hz]
#     OUTPUT::
#     ------
#     fmax : float
#         Maximum frequency [Hz]
PRF/2
}

DopplerVmax<- function(PRF, lam){
#     Nyquist velocity, or maximum unambiguous Doppler velocity (+ or -).
#     From Rinehart (1997), Eqn 6.7
#     INPUT::
#     -----
#     PRF : float
#         Radar pulse repetition frequency [Hz]
#     lam : float
#         Radar wavelength [m]
#     OUTPUT::
#     ------
#     Vmax : float
#         Nyquist velocity [m/s], +/-
PRF * lam / 4
}

DopplerRmax<- function(PRF,speedOfLight){
#     Maximum unamiguous range.
#     From Rinehart (1997), Eqn 6.11
#     INPUT::
#     -----
#     PRF : float
#         Pulse repetition frequency [Hz]
#     OUTPUT::
#     ------
#     Rmax : float
#         Maximum unambiguous range [m]
speedOfLight / (2 * PRF)
}

DopplerDilemma<- function(inFloat, lam,speedOfLight){
#     The "Doppler dilemma" is the fact that both the Nyquist velocity and 
#       unambiguous maximum range of the radar are based upon the PRF of the system.
#     However, they are inversely proportional, meaning that increasing one 
#       requires a decrease in the other.  A trade-off inherent in Doppler radar
#       systems.  This relationship allows a solution for one variable given the
#       value of the other.
#     From Rinehart (1997), Eqn 6.12
#     INPUT::
#     -----
#     inFloat : float
#         Nyquist Velocity [m/s] or Maximum unambiguous range [m]
#     lam : float
#         Radar wavelength [m]
#     OUTPUT::
#     ------
#     Out : float
#         The inFloat that is not used
(speedOfLight * lam / 8) / inFloat
}

DopplerVshift<- function(GS, psi){
#     Adjusted Doppler velocity from a mobile platform.
#     From Jorgensen (1983), Eqn 2
#     INPUT::
#     -----
#     GS : float
#         Gound speed [m/s]
#     psi : float
#         Angle between actual azimuth and fore/aft angle [deg]
#     OUTPUT::
#     ------
#     Vshift : float
#         Shift in Doppler velocity from mobile aspect [m/s]
GS * cos(pi/360*(psi))
    }
    
DopplerVmaxDual<- function(lam, PRF1, PRF2){
#     Doppler velocity from dual PRF scheme radar (+ or -).
#     From Jorgensen (1983), Eqn 2
#     INPUT::
#     -----
#     lam : float
#         Radar wavelength [m]
#     PRF1 : float
#         First Pulse repetition frequency [Hz]
#     PRF2 : float
#         Second Pulse repetition frequency [Hz]
#     OUTPUT::
#     ------
#     Vmax : float
#         Doppler velocity [m/s]
lam / (4 * ((1 / PRF1) - (1 / PRF2)))
    }

GeometryReffective<- function(dNdH=-39e-6, earthRadius){
#     Effective radius calculation.
#     Rinehart (1997), Eqn 3.9, solved for R'
#     INPUT::
#     -----
#     dNdH : float
#         Refraction [N x10^-6/km]
#     OUTPUT::
#     -----
#     R1 : float
#         Effective radius [m]
(1 / ((1/(earthRadius/1000)) + (dNdH))) * 1000
 }

GeometryHalfPowerRadius<- function(r, bwhalf){
#     Half-power radius.
#     Battan (1973), 
#     INPUT::
#     -----
#     r : float
#         Range [m]
#     bwhalf : float
#         Half-power beam width [degrees]
#     OUTPUT::
#     -----
#     Rhalf : float
#         Half-power radius [m]
(r * pi/360*(bwhalf)) / 2
}

GeometryRayHeight<- function(r, elev, H0, R1=kConstantR43){
#     Center of radar beam height calculation.
#     Rinehart (1997), Eqn 3.12, Bech et al. (2003) Eqn 3
#     INPUT::
#     -----
#     r : float
#         Range from radar to point of interest [m]
#     elev : float
#         Elevation angle of radar beam [deg]
#     H0 : float
#         Height of radar antenna [m]
#     R1 : float
#         Effective radius
#     OUTPUT::
#     -----
#     H : floa                 t
#         Radar beam height [m]
    # Convert earth's radius to km for common dN/dH values and then
    # multiply by 1000 to return radius in meters
    Term1 = sqrt(r^2 +R1^2 + 2*r*R1*sin(pi/360*(elev)))
Term1 - R1 + H0
}

GeometrySampleVolIdeal<- function(r, bwH, bwV, pLength){
#     Sample volume (idealized) assuming all power in half-power beamwidths.
#     From Rinehart (1997), Eqn 5.2
#     INPUT::
#     -----
#     r : float
#         Distance to sample volume from radar [m]
#     bwH : float
#         Horizontal beamwidth [deg]
#     bwV : float
#         Vertical beamwidth deg]
#     pLength : float
#         Pulse length [m]
#     OUTPUT::
#     -----
#     SVol : float
#         Sample Volume [m^3]
pi * (r * pi/360*(bwH)/2) * (r * pi/360*(bwV)/2) * (pLength/2)
 }
 
GeometrySampleVolGauss<- function(r, bwH, bwV, pLength){
#     Sample volume assuming transmitted energy in Gaussian beam shape.
#     From Rinehart (1997), Eqn 5.4
#     INPUT::
#     -----
#     r  : float
#         Distance to sample volume from radar [m]
#     bwH : float
#         Horizontal beamwidth [deg]
#     bwV : float
#         Vertical beamwidth deg]
#     pLength : float
#         Pulse length [m]
#     OUTPUT::
#     -----
#     SVol : float
#         Sample Volume [m]
    Numer = pi * r^2 * pi/360*(bwH) * pi/360*(bwV) * pLength
    Denom = 16 * log(2)
Numer / Denom
}

GeometryRangeCorrect<- function(r, h, E){
#     A corrected range from radar that takes into account the "loss" of 
#       ground distance because of the radar elevation angle.  This is a 
#       cumulative effect at each gate along the ray.
#     From CSU Radar Meteorology AT 741 Notes
#     INPUT::
#     -----
#     r  : float
#         Distance to sample volume from radar [m]
#     h : float
#         Height of the center of radar volume [m]
#     E : float
#         Elevation angle [deg]
#     OUTPUT::
#     -----
#     rnew : float
#         Adjusted range to sample volume [m]
    # Calculate the change in height along the ray
n<-length(h)
    dh1 = h[2:n] - h[1:(n-1)]
    # Add the 0th place in the ray at the beginning
    dh2 = c(h[1],dh1)
    # Calculate the change in distance at each gate
    a90r = pi/2  # 90 degrees in radians
    dr = dh2 / (tan(a90r - pi/360*(E)))
    # Now calculate the corrected range at each gate
cumsum(dr)
}

GeometryBeamBlockFrac<- function(Th, Bh, a){
#     Partial beam blockage fraction.
#     From Bech et al. (2003), Eqn 2 and Appendix
#     INPUT::
#     -----
#     Th : float
#         Terrain height [m]
#     Bh : float
#         Beam height [m]
#     a : float
#         Half power beam radius [m]
#     OUTPUT::
#     -----
#     PBB : float
#         Partial beam blockage fraction [unitless]
    # First find the difference between the terrain and height of
    # radar beam (Bech et al. (2003), Fig.3)
    y = Th - Bh
    Numer = (y * sqrt(a^2 - y^2)) + (a^2 * asin(y/a)) + (pi * a^2 /2)
    Denom = pi * a^2
Numer / Denom
}

SystemGainPratio<- function(P1, P2){
#     Antenna gain via power ratio.  
#     From Rinehart (1997), Eqn 2.1
#     INPUT::
#     -----
#     P1 : float
#         Power on the beam axis [W]
#     P2 : float
#         Power from an isotropic antenna [W]
#     OUPUT::
#     -----
#     G : float
#         Gain [dB]
10 * log10(P1 / P2)
}

SystemFreq<- function(lam, speedOfLight){
#     Frequency given wavelength.
#     INPUT::
#     -----
#     lam : float
#         Wavelength [m]
#     OUPUT::
#     -----
#     freq : float
#         Frequency [Hz]
speedOfLight / lam
}

Systemwavelength<- function(freq, speedOfLight){
#     Wavelength given frequency.
#     INPUT::
#     -----
#     freq : float
#         Frequency [Hz]
#     OUPUT::
#     -----
#     lam : float
 speedOfLight / freq
}

SystemRadarConst<- function(Pt, G, Tau, lam, bwH, bwV, Lm, Lr){
#     Radar constant.
#     From CSU Radar Meteorology notes, AT 741
#     INPUT::
#     -----
#     Pt : float
#         Transmitted power [W]
#     G : float
#         Antenna Gain [dB]
#     Tau : float
#         Pulse Width [s]
#     lam : float
#         Radar wavelength [m]
#     bwH : float
#         Horizontalntenna beamwidth [degrees]
#     bwV : float
#         Vertical antenna beamwidth [degrees]
#     Lm : float
#         Antenna/waveguide/coupler loss [dB]
#     Lr : float
#         Receiver loss [dB]
#     OUPUT::
#     -----
#     C : float
#         Radar constant [unitless]
    # Convert from dB to linear units
    Lmlin = 10^(Lm / 10)
    Lrlin = 10^(Lr / 10)
    Glin = 10^(G / 10)
    # Convert beamwidth to radians
    bwHr = pi/360*(bwH)
    bwVr = pi/360*(bwV)
    # Calculate the numerator
    Numer = pi^3 * c * Pt * Glin^2 * Tau * bwHr * bwVr * Lmlin * Lrlin
    # Calculate the denominator
    Denom = 1024 *log(2) * lam^2
Numer/Denom
}

SystemAntEffArea<- function(G, lam){
#     Antenna effective area.
#     From Rinehart (1997), Eqn 4.5
#     INPUT::
#     -----
#     G : float
#         Antenna Gain [dB]
#     lam : float
#         Radar wavelength [m]
#     OUPUT::
#     -----
#     Ae : float
#         Antenna effective area [unitless]
    # Convert from dB to linear units
Glin = 10^(G / 10)
Glin * lam^2 / (4 * pi)
}

SystemPowerTarget<- function(Pt, G, Asig, r){
#     Power intercepted by target.
#     From Rinehart (1997), Eqn 4.3
#     INPUT::
#     -----
#     Pt : float
#         Transmitted power [W]
#     G : float
#         Antenna gain [dB]
#     Asig : float
#         Area of target [m^2]
#     r : float
#         Distance to sample volume from radar [m]
#     OUPUT::
#     -----
#     Psig : float
#         Power intecepted by target [m]
    # Convert from dB to linear units
    Glin = 10^(G / 10)
(Pt * Glin * Asig) / (4 * pi * r^2)
}

SystemXsecBscatterSphere<- function(D, lam, K=0.93){
#     Backscatter cross-sectional area of a sphere using the Rayleigh approximation.
#     From Rinehart (1997), Eqn 4.9 and 5.7
#     INPUT::
#     -----
#     D : float
#         Diamter of targer [m]
#     lam : float
#         Radar wavelength [m]
#     K : float
#         Dielectric factor [unitless]
#     OUPUT::
#     -----
#     sig : float
#         Backscattering cross-section [m*2]
(pi^5 * K^2 * D^6) / lam^4
}

SystemNormXsecBscatterSphere<- function(D, lam, K=0.93){
#     Normalized Backscatter cross-sectional area of a sphere using the Rayleigh approximation.
#     From Rinehart (1997), Eqn 4.9 and 5.7 and Battan Ch. 4.5
#     INPUT::
#     -----
#     D : float
#         Diamter of targer [m]
#     lam : float
#         Radar wavelength [m]
#     K : float
#         Dielectric factor [unitless]
#     OUPUT::
#     -----
#     sigNorm : float
#         Normalized backscatter cross-section [unitless]
    # Calculate the cross-sectional backscatter area
    sig = SystemXsecBscatterSphere(D, lam, K)
sig/ (pi * (D/2)^2)
}

SystemSizeParam<- function(D, lam){
#     Size parameter calculation.
#     From Rinehart (1997), Eqn 4.9 and 5.7 and Battan Ch. 4.5
#     INPUT::
#     -----
#     D : float
#         Diamter of targer [m]
#     lam : float
#         Radar wavelength [m]
#     OUPUT::
#     -----
#     alpha : float
#         Size parameter [unitless]
2 * pi * D/2 / lam
}

SystemPowerReturnTarget<- function(Pt, G, lam, sig, r){
#     Power returned by target located at the center of the antenna beam pattern.
#     From Rinehart (1997), Eqn 4.7
#     INPUT::
#     -----
#     Pt : float
#         Transmitted power [W]
#     G : float
#         Antenna gain [dB]
#     lam : float
#         Radar wavelength [m]
#     sig : float
#         Backscattering cross-sectional area of target [m^2]
#     r : float
#         Distance to sample volume from radar [m]
#     OUPUT::
#     -----
#     Pr : float
#         Power returned by target [m]
    # Convert from dB to linear units
    Glin = 10^(G/10)
(Pt * Glin^2 * lam^2 * sig) / (64 * pi^3 * r^4)
}

SystemThermalNoise<- function(Bn, Units, Ts=290, k=kConstantBoltz){
#     Thermal noise power.
#     From CSU Radar Meteorology notes, AT741
#     INPUT::
#     -----
#     Bn : float
#         Receiver bandwidth [Hz]
#     Units : float
#         String of nits desired, can be 'W' or 'dBm'
#     Ts : float
#         Reciever noise temperature [K]
#     OUPUT::
#     -----
#     Nt : float
#         Thermal noise power [W or 'dBm']
     # Calculate the noise, convert if requested
    N = k * Ts * Bn
    if (toupper(Units)=='W') Nt = N else {
    if (toupper(Units)=='DBM') Nt = 10 * log10(N/10^-3) else {
        warning( "Units must be in 'W' or 'dBm'")
        Nt = NA
        }
        }
Nt
}

VariablesReflectivity<- function(Pt, G, Tau, lam, bwH, bwV, Lm, Lr, Pr, r, K=0.93){
#     Radar reflectivity.
#     From Rinehart (1993), Eqn 5.17 (See Eqn 5.14-5.16 also)
#     INPUT::
#     -----
#     Pt : float
#         Transmitted power [W]
#     G : float
#         Antenna Gain [dB]
#     Tau : float
#         Pulse Width [s]
#     lam : float
#         Radar wavelength [m]
#     bwidth : float
#         Antenna beamwidth [degrees]
#     Lm : float
#         Antenna/waveguide/coupler loss [dB]
#     Lr : float
#         Receiver loss [dB]
#     K : float
#         Dielectric factor [unitless]
#     Pr : float
#         Returned power [W]
#     r : float
#         Range to target [m]
#     OUTPUT::
#     -----
#     Ze : float
    # Call the radar constant function
    C1 = SystemRadarConst(Pt, G, Tau, lam, bwH, bwV, Lm, Lr)
Pr * r^2 / (C1 * K^2)
}

VariablesRadVel<- function(f,lam){
#     Radial velocity.
#     From Rinehart (1993), Eqn 6.6
#     INPUT::
#     -----
#     f : float
#         Frequency shift [Hz]
#     lam : float
#         Radar wavelength [m]
#     OUTPUT::
#     -----
#     Vr : float
#         Radial velocity [m/s]
 f * lam / 2
 }

VariablesCDR<- function(Zpar, Zorth){
#     Circular depolarization ratio.  
#     From Rinehart (1997), Eqn 10.2
#     INPUT::
#     -----
#     Zpar : float
#         Reflectivity in the parallel channel [mm^6/m^3]
#     Zorth : float
#         Reflectivity in the orthogonal channel [mm^6/m^3]
#     OUTPUT::
#     -----
#     CDR : float
#         Circular depolarization ratio [dB]
10* log10(Zpar/Zorth)
}

VariablesLDR<- function(Zh, Zv){
#     Linear depolarization ratio.
#     From Rinehart (1997), Eqn 10.3
#     INPUT::
#     -----
#     Zh : float
#         Horizontal reflectivity [mm^6/m^3]
#     Zv : float
#         Vertical reflectivity [mm^6/m^3]
#     OUTPUT::
#     -----
#     LDR : float
10 * log10(Zh / Zv)
}

VariablesZDR<- function(Zh, Zv){
#     Differential reflectivity.
#     From Rinehart (1997), Eqn 10.3 and Seliga and Bringi (1976)
#     INPUT::
#     -----
#     Zh : float
#         Horizontal reflectivity [mm^6/m^3]
#     Zv : float
#         Vertical reflectivity [mm^6/m^3]
#     OUTPUT::
#     -----
#     ZDR : float
#         Differential reflectivity [dB]
10 * log10(Zh / Zv)
}

VariablesZDP<- function(Zh, Zv){
#     Reflectivity difference.
#     From Rinehart (1997), Eqn 10.3
#     INPUT::
#     -----
#     Zh : float
#         Horizontal reflectivity [mm^6/m^3]
#     Zv : float
#         Vertical reflectivity [mm^6/m^3]
#     OUTPUT::
#     -----
#     ZDP : float
#         Reflectivity difference [dB]
    if (Zh > Zv){
        ZDP = 10* log10(Zh - Zv)
    } else{
        warning("Zh < Zv !")
        ZDP = NA
        }
ZDP
}

VariablesHDR<- function(dBZh, ZDR){
#     Differential reflectivity hail signature.
#     From Aydin et al. (1986), Eqns 4-5
#     INPUT::
#     -----
#     Zh : float
#         Horizontal reflectivity [dBZ]
#     ZDR : float
#         Differential reflectivity [dBZ]
#     OUTPUT::
#     -----
#     ZDP : float
#         Reflectivity difference [dB]
    # Set the f(Zdr) based upon observations
    if (ZDR <= 0) f = 27 else {
    if ((ZDR > 0) & (ZDR <= 1.74)) f = 19 * ZDR + 27 else if (ZDR > 1.74) f = 60
    }
    # Calculate HDR
dBZh - f
}
