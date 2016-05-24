# Modified: 24 JULY 2015 SDH
 
partCO2 <-
function(
  nCO2,     # Total moles of CO2
  nCH4,     # Total moles of CH4
  mass.H2O, # Mass of water in solution in kg
  temp.c,   # Temperature in degrees C
  pH,       # pH of solution
  mu,       # Ionic strength of solution (mol/kgw)
  pres = 101325,   # Total biogas pressure (Pa)
  method = 'cont', # Method to be used, 'cont' = for continuous reactor, 'batch' = batch reactor after infinite time
  value = 'all'
  ) {

  # Activity of H+
  a.H <- 10^(-pH)

  # Temperature in K
  temp.k <- temp.c + 273.15 

  # Constants. Taken from Hafner et al. 2012 model
  # kh, k1, and k2 are from Plummer & Busenberg 1982
  # Chemical species are CO2 (aq) (includes H2CO3), HCO3-, and CO3-2
  # Henry's law constant (originally mol/kg-atm, note conversion to mol/kg-Pa as the second term)
  kh <- 10^(108.3865 + log10(1/101325) + 0.01985076*temp.k - 6919.53/temp.k - 40.45154*log10(temp.k) + 669365/temp.k^2)
  # Equilibrium constants
  k1 <- 10^(-356.3094 - 0.06091964*temp.k + 21834.37/temp.k + 126.8330*log10(temp.k) - 1684915/temp.k^2)
  k2 <- 10^(-107.8871 - 0.03252849*temp.k + 5151.79/temp.k + 38.92561*log10(temp.k) - 563713.9/temp.k^2)
  # Universal gas constant
  #R <- 0.082057  # Universal gas constant (atm L/mol-K)
  R <- 8314.4621  # Universal gas constant (Pa L/mol-K)

  # Calculate a.dh and b.dh parameters for Debye-Huckel equation. Dielectric constant, density rearranged from PHREEQC code
  # Dielectric constant
  de <- 2727.586 + 0.6224107*temp.k - 1075.112*log10(temp.k) - 52000.87/temp.k
  # Water density
  c.d <- 647.26 - temp.k
  d.H2O <- (1 + .1342489*c.d^(1/3) - 3.946263E-3*c.d)/(3.1975 - 0.3151548*c.d^(1/3) - 1.203374E-3*c.d + 7.48908E-13*c.d^4)
  # a.dh and b.dh are A and B in Debye-Huckel equation. Equations are from Truesdell & Jones (1974)
  a.dh <- 1.82483E6*d.H2O^0.5/(de*temp.k)^(3/2)
  b.dh <- 50.2916*d.H2O^0.5/(de*temp.k)^0.5

  # Activity coefficients
  gHCO3<- 10^(-a.dh*1^2*sqrt(mu)/(1+b.dh*5.4*sqrt(mu)))
  gCO3<- 10^(-a.dh*2^2*sqrt(mu)/(1+b.dh*4.5*sqrt(mu)))

  # Calculate saturated water vapor pressure in Pa
  #P.H2O <- 10^(5.08600 - 1668.21/(temp.c + 228.0))
  P.H2O <- watVap(temp.k = temp.k)

  #NTS check results more closely
  if(method=='cont') {
    # Initial guess for CO2 in biogas is total amount
    vBg <- (nCO2 + nCH4)*R*temp.k/(pres - P.H2O)  
    # Can check above expression with:
    #nH2O <- P.H2O*vBg/(R*temp.k)
    #vBg <- (nCO2 + nCH4 + nH2O)*R*temp.k/pres

    # dvol = delta vol in an iteration
    dvol <- 1
    tstart <- Sys.time()
    while(dvol > 1E-6*vBg) {
      # Molal concentrations of aqueous species at equilibrium, first equation was derived from mass action expressions and mass balance
      mCO2 <- nCO2/(1 + vBg/(kh*R*temp.k*mass.H2O) + k1/(a.H*gHCO3) + k1*k2/(a.H^2*gCO3))/mass.H2O
      mHCO3 <- 1.0*mCO2*k1/(a.H*gHCO3)
      mCO3 <- gHCO3*mHCO3*k2/(a.H*gCO3)

      # Update estimate of CO2 in biogas
      nCO2.sol <- (mCO2 + mHCO3 + mCO3)*mass.H2O
      nCO2Bg <- nCO2 - nCO2.sol

      # Update estimate of total biogas volume 
      vBg1 <- vBg
      vBg <- (nCO2Bg + nCH4)*R*temp.k/(pres - P.H2O)
      dvol <- abs(vBg - vBg1)

      # CO2 partial pressure, not returned but for checking if needed
      PCO2 <- nCO2Bg*R*temp.k/vBg
      if(as.numeric(Sys.time() - tstart, units = 'secs') > 5) {
	stop('The while loop for solution speciation seems to be stuck. Try your call without pH, conc.sub, or temp arguments. Sorry.')
      }
    } 
    # Check mass balance based on PCO2 (not needed but good to run if above equations are changed)
    #if(abs(mCO2/kh - nCO2Bg*R*temp.k/vBg)/(mCO2/kh) > 1E-6) stop('Problem with CO2 partitioning--mass balance off.')

  } else if(method=='batch') {
    # CO2 partial pressure here is at infinite time, when solution is in equlibrium with biogas as produced, not a constant or average value
    PCO2 <- nCO2/(nCO2 + nCH4)*(1 - P.H2O)
    # Concentrations of aqueous species from partial pressure of CO2
    mCO2 <- PCO2*kh
    mHCO3 <- mCO2*k1/(a.H*gHCO3)
    mCO3 <- gHCO3*mHCO3*k2/(a.H*gCO3)

    # Moles CO2 in biogas and solution
    nCO2.sol <- (mCO2 + mHCO3 + mCO3)*mass.H2O
    if(nCO2.sol>nCO2) nCO2.sol <- nCO2
    nCO2Bg <- nCO2 - nCO2.sol
  }

  # TIC in solution (mol/kgw)
  cTIC <- mCO2 + mHCO3 + mCO3

  # Calculate xCO2 and xCH4 (defined for dry biogas)
  xCO2 <- nCO2Bg/(nCO2Bg + nCH4)
  xCH4 <- nCH4/(nCO2Bg + nCH4)

  # Return xCH4 
  if(value=='xCH4') return(xCH4)

  # Complete results. as.vector just to drop any names that will be added to returned vector names
  # NTS: is there a better way than as.vector?
  return(c(nCO2Bg = as.vector(nCO2Bg), nCO2.sol = as.vector(nCO2.sol), cTIC = as.vector(cTIC), xCH4 = as.vector(xCH4)))
}
