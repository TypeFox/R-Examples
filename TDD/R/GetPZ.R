GetPZ = function(w){
  # initialize output
  PZ = list()

  ## Broadband Sensors:

  # STS-1 (360 s)
  PZ$STS1_360 = list(poles = c(-1.234e-2 + 1.234e-2i, -1.234e-2 - 1.234e-2i, -3.918e1 + 4.912e1i, -3.918e1 - 4.912e1i), np = 4, zeros = c(0, 0), nz = 2, Knorm = 3.948580E+03, Sense = 2400) # http://www.iris.edu/NRL/sensors/streckeisen/RESP.XX.NS088..BHZ.STS1.360.2400

  # Trillium 240 SN <400
  PZ$Trillium240_gen1 = list(poles = c(-0.01813 + 0.01803i, -0.01813 - 0.01803i, -124.9, -197.5 + 256.1i, -197.5 - 256.1i, -569 + 1150i, -569 - 1150i), np = 7, zeros = c(0, 0, -90, -164.2, -3203), nz = 5, Knorm = 4.532e5, Sense = 1196.6) # http://www.passcal.nmt.edu/webfm_send/1969
  # Trillium 240 SN >=400
  PZ$Trillium240_gen2 = list(poles = c(-0.0177 + 0.0176i, -0.0177 - 0.0176i, -126.7, -192 + 259.1i, -192 - 259.1i, -557.7 + 1143i, -557.7 - 1143i), np = 7, zeros = c(0, 0, -91.66, -160.1, -3207), nz = 5, Knorm = 4.517e5, Sense = 1168.2) # http://www.passcal.nmt.edu/webfm_send/1969
  
  # 3T
  PZ[['3T']] = list(poles = c(-0.037 + 0.037i, -0.037 - 0.037i, -503, -1010, -1130), np = 5, zeros = c(0, 0), nz = 2, Knorm = 2304000 * (2*pi)^3, Sense = 1500) # http://www.guralp.com/poles-and-zeroes-with-positive-normalization-factors/ (given in Hz, not rad/s)
  
  # STS-2--many generations--http://www.iris.edu/NRL/sensors/streckeisen/streckeisen_sts2_sensors.htm
  PZ$STS2_gen1 = list(poles = c(-0.037 + 0.037i, -0.037 - 0.037i, -15.99, -417.1, -187.239, -100.9 + 401.9i, -100.9 - 401.9i, -7454 + 7142i, -7454 - 7142i), np = 9, zeros = c(0, 0, -15.15, -318.6 + 401.2i, -318.6 - 401.2i), nz = 5, Knorm = 5.70624e12, Sense = 1500) # http://www.iris.edu/NRL/sensors/streckeisen/RESP.XX.NS081..BHZ.STS2_gen1.120.1500
  
  PZ$STS2_gen2 = list(poles = c(-3.7E-02 - 3.7E-02i, -3.7E-02 + 3.7E-02i, -1.095E+01, -9.844E+01 - 4.428E+02i, -9.844E+01 + 4.428E+02i, -5.568E+02 - 6.005E+01i, -5.568E+02 + 6.005E+01i, -1.391E+03, -4.936E+03 - 4.713E+03i, -4.936E+03 + 4.713E+03i, -6.227E+03, -6.909E+03 - 9.208E+03i,-6.909E+03 + 9.208E+03i, -2.550970E+02), np = 14, zeros = c(0, 0, -1.075E+01, -2.946E+02, -5.551E+02, -6.839E+02 - 1.755E+02i, -6.839E+02 + 1.755E+02i, -5.907E+03 - 3.411E+03i, -5.907E+03 + 3.411E+03i), nz = 9, Knorm = 2.35238E+17, Sense = 1500) # http://www.iris.edu/NRL/sensors/streckeisen/RESP.XX.NS083..BHZ.STS2_gen2.120.1500

  PZ$STS2_gen3 = list(poles = c(-3.700E-02 - 3.700E-02i, -3.700E-02 + 3.700E-02i, -1.564E+01, -9.734E+01 - 4.007E+02i, -9.734E+01 + 4.007E+02i, -3.748E+02, -5.203E+02, -1.053E+04 - 1.005E+04i, -1.053E+04 + 1.005E+04i, -1.330E+04, -2.550970E+02), np = 11, zeros = c(0, 0, -15.15, -176.6, -463.1 - 430.5i, -463.1 + 430.5i), nz = 6, Knorm = 3.4684E+17, Sense = 1500) # http://www.iris.edu/NRL/sensors/streckeisen/RESP.XX.NS085..BHZ.STS2_gen3.120.1500
                        
  
  # Trillium 120
  
PZ$Trillium120 = list(poles = c(-3.852e-2 + 3.658e-2i, -3.852e-2 - 3.658e-2i, -1.78e2, -1.35e2 + 1.6e2i, -1.35e2 - 1.6e2i, -6.71e2 + 1.154e3i, -6.71e2 - 1.154e3i), np = 7, zeros = c(0, 0, -0.9e2, -1.607e2, -31.08e2), nz = 5, Knorm = 3.080e5, Sense = 1201) # http://www.passcal.nmt.edu/webfm_send/1970
  
  # Compact Trillium
  PZ$TrilliumC = list(poles = c(-3.691e-2 + 3.712e-2i, -3.691e-2 - 3.712e-2i, -3.712e2, -3.739e2 + 4.755e2i, -3.739e2 - 4.755e2i, -5.884e2 + 1.508e3i, -5.884e2 - 1.508e3i), np = 7, zeros = c(0, 0, -4.341e2), nz = 3, Knorm = 8.184e11, Sense = 749.1) # http://www.iris.edu/NRL/sensors/nanometrics/RESP.XX.NS124..BHZ.TrilliumCompact.120.749
  
  ## Intermediate Sensors:

  # Trillium 40
PZ$Trillium40 = list(poles = c(-2.41e2 - 1.78e2i, -2.41e2 + 1.78e2i, -5.35e2 - 7.19e2i, -5.35e2 + 7.19e2i, -8.63e1, -1.103e-1 + 1.11e-1i, -1.103e-1 - 1.11e-1i), np = 7, zeros = c(0, 0, -6.88e1, -3.23e2, -2.53e3), nz = 5, Knorm = 1.104e5, Sense = 1553) # http://www.iris.edu/NRL/sensors/nanometrics/RESP.XX.NS120..BHZ.Trillium.40.1553
  
  # 3ESP
  PZ[['3ESP']] = list(poles = c(-1.49e-1 + 1.49e-1i, -1.49e-1 - 1.49e-1i, -503, -1010, -1130), np = 5, zeros = c(0, 0), nz = 2, Knorm = 2304000 * (2*pi)^3, Sense = 2000) # http://www.guralp.com/poles-and-zeroes-with-positive-normalization-factors/ (given in Hz, not rad/s)
  
  # 40T (30 s)
  PZ[['40T']] = list(poles = c(-1.49e-1 + 1.49e-1i, -1.49e-1 - 1.49e-1i, -503, -1010, -1130), np = 5, zeros = c(0, 0), nz = 2, Knorm = 2304000 * (2*pi)^3, Sense = 800) # http://www.guralp.com/poles-and-zeroes-with-positive-normalization-factors/ (given in Hz, not rad/s)
  
  # STS-1 (20 s)
  PZ$STS1_20 = list(poles = c(-2.221e-1 + 2.221e-1i, -2.221e-1 - 2.221e-1i, -3.918e1 + 4.912e1i, -3.918e1 - 4.912e1i), np = 4, zeros = c(0, 0), nz = 2, Knorm = 3.948580E+03, Sense = 2400) # http://www.iris.edu/NRL/sensors/streckeisen/RESP.XX.NS087..BHZ.STS1.20.2400
  
  ## Short-period Sensors:

  # 40T-1
PZ[['40T1']] = list(poles = c(-0.4443e1 + 0.4443e1i, -0.4443e1 - 0.4443e1i, -0.3920e3 + 0.8507e3i, -0.3920e3 - 0.8507e3i, -0.2199e4, -0.4712e3), np = 6, zeros = c(0, 0), nz = 2, Knorm = 585800000 * (2*pi)^4, Sense = 2000) # http://www.guralp.com/poles-and-zeroes-with-positive-normalization-factors/ (given in Hz, not rad/s)

  ## DONE!
  return(PZ[w])
}
