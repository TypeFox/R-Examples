# Modified: 2 April 2016 SDH

predBg <- function(
  form = NULL,            # Character chemical formula of substrate
  mass = 1,        # Mass of substrate (g)
  mol = NULL,      # Moles of substrate 
  fs = 0,          # Fraction substrate going to cell synthesis, fs in Rittmann and McCarty
  fd = 1,          # Biodegradable fraction of substrate (g/g)
  mcomp = NULL,    # Macromolecular composition, vector with named elements NTS add details
  COD = NULL,      # Substrate total COD (g) (alternative to form or mcomp)
  conc.sub = NULL, # Substrate concentration in g per kg water (only needed for CO2 partitioning)
  pH = NULL,       # Solution pH (only needed for CO2 partitioning) 
  temp = NULL,     # Temperature in degrees C (only needed for CO2 partitioning) 
  mu = 0.1,        # Solution ionic strength (mol/kgw) (only needed for CO2 partitioning)
  shortform = NULL,       # TRUE to simplify form
  value = 'CH4'    # Output format, default is CH4 volume (mL)
  ){

  # Check arguments
  checkArgClassValue(form, c('character', 'NULL'))
  checkArgClassValue(mass, c('numeric', 'integer'), expected.range = c(0, Inf))
  checkArgClassValue(mol, c('numeric', 'integer', 'NULL'), expected.range = c(0, Inf))
  checkArgClassValue(fs, c('numeric', 'integer'), expected.range = c(0, 1))
  checkArgClassValue(fd, c('numeric', 'integer'), expected.range = c(0, 1))
  checkArgClassValue(mcomp, c('numeric', 'integer', 'NULL'))
  checkArgClassValue(COD, c('numeric', 'integer', 'NULL'), expected.range = c(0, Inf))
  checkArgClassValue(conc.sub, c('numeric', 'integer', 'NULL'), expected.range = c(0, Inf))
  checkArgClassValue(pH, c('numeric', 'integer', 'NULL'), expected.range = c(4, 10))
  checkArgClassValue(temp, c('numeric', 'integer', 'NULL'))
  checkArgClassValue(mu, c('numeric', 'integer', 'NULL'), expected.range = c(0, 0.8))
  checkArgClassValue(value, 'character')
  checkArgClassValue(tolower(value), expected.values = c('ch4', 'reactionn', 'reactionc', 'all'))

  # fe from Rittmann and McCarty (2001)
  fe <- 1 - fs
  
  # If COD is given and there is no form, just return CH4 volume in mL
  if(!missing(COD) & missing(form)) {
    #vCH4 <- fd*fe*COD*vol.mol['CH4']*0.5/(molMass('O')*2)
    vCH4 <- fd*fe*COD*vol.mol['CH4']/(molMass('O2')*2)
    vCH4 <- as.vector(vCH4)
    if(tolower(value)=='ch4') {
      return(vCH4)
    } else {
      return(data.frame(COD = COD, fd = fd, fs = fs, fe = fe, vCH4 = vCH4))
    }
  } 

  # If COD and form are both given, calculate mass and mol
  if(!missing(COD) & !missing(form)) {
    mass <- COD/calcCOD(form)
  }

  # Read chemical formula
  if(!missing(form)) {
    if(missing(shortform)) shortform <- FALSE
    if(!missing(mol) & shortform) warning('You specified quantity using the \'mol\' argument but also asked for chemical formula simplification using \'shortform\'.\nNote that the \'mol\' quantity is applied using the new (returned) formula.')
  
    # Capitalize form
    for(i in 1:length(form)) {
      if(grepl('^[a-z]', form[i])) form[i] <- toupper(form[i])
    }

  } else if(!missing(mcomp)) {
    if(missing(shortform)) shortform <- TRUE
    if(!missing(mol) & shortform) warning('You specified quantity using the \'mol\' argument but also asked for chemical formula simplification using \'shortform\'.\nNote that the \'mol\' quantity is applied using the new (returned) formula.')

    # Make sure it is a vector. Cannot provide multiple values for mcomp.
    if(!is.vector(mcomp)) stop('mcomp can only be a vector, and can only be used to specify a single composition.')

    if(abs(sum(mcomp)-1)>1E-5) {

      # Calculate mass if not given (override default)
      if(missing(mass)) {
	mass <- sum(mcomp)
	warning('mass taken as sum of mcomp: ', mass)
      }

      warning('Sum of mcomp != 1.0 so dividing all elements by the sum for calculation of formula.')
      mcomp <- mcomp/sum(mcomp)
    }

    # Code below for inert elements in mcomp has never been implemented. Would required using masscor below
    # Not sure if it will work when mol is set
    ### If there is any inert element, a mass correction factor is calculated
    ##masscor <- 1 - sum(mcomp[!tolower(names(mcomp)) %in% c('inert', 'ash', 'nondegradable', 'non-degradable')])/sum(mcomp)

    # Trim std.forms to use it in indexing operation to replace names in mcomp
    mforms <- std.forms[names(std.forms) %in% tolower(names(mcomp))]
    
    # Names must be formulas or else
    # Replace mcomp names like 'carbohydrate' with formulas
    mcnames <- names(mcomp)
    names(mcnames) <- tolower(names(mcomp))
    mcnames[names(mforms)] <- mforms
    names(mcomp) <- mcnames

    # Get molar masses and convert mass values in mcomp to moles
    mcomp <- mcomp/molMass(names(mcomp))
    mcomp <- signif(mcomp/min(mcomp), 6)

    # Create form from mcomp
    names(mcomp) <- paste0('(', names(mcomp), ')')
    form <- paste0(names(mcomp), mcomp, collapse = '')
    for(i in 1:length(form)) {
      if(grepl('^[a-z]', form[i])) form[i] <- toupper(form[i])
    }

  } else stop('Must provide one of these arguments: form, mcomp, or COD.')

  # All options except COD only continue

  # Simplify formula
  if(shortform) form <- sapply(form, readFormula, USE.NAMES = FALSE, value = 'shortform')

  # Get formula coefficients
  fc <- t(mapply(readFormula, form = form, MoreArgs = list(elements = c('C', 'H', 'N', 'O'), min.elements = c('C', 'H')), USE.NAMES = FALSE))

  # Molar mass of substrate
  mmass <- mapply(molMass, form, USE.NAMES = FALSE)

  # Moles and mass of substrate, if mol is provided, it will override default mass
  if(missing(mol)) mol <- mass/mmass else mass <- mol*mmass

  # Calculate COD'
  COD <- mapply(calcCOD, form, USE.NAMES = FALSE)*mass

  # Coefficients in overall reaction in Eq. 13.5 in Rittmann and McCarty (2001)
  d <- as.vector(4*fc[, 'C'] + fc[, 'H'] - 2*fc[, 'O'] - 3*fc[, 'N'])
  cH2O <- as.vector(2*fc[, 'C'] + fc[, 'N'] - fc[, 'O'] - 9*d*fs/20 - d*fe/4)
  cCH4 <- as.vector(d*fe/8)
  cCO2 <- as.vector(fc[, 'C'] - fc[, 'N'] - d*fs/5 - d*fe/8)
  cbio <- as.vector(d*fs/20)
  cNH4 <- as.vector(fc[, 'N'] - d*fs/20)
  cHCO3 <- as.vector(fc[, 'N'] - d*fs/20)

  # Results
  # CH4 in mL at 0C and 1 atm
  nCH4 <- cCH4*fd*mol
  vCH4 <- vol.mol['CH4']*nCH4
  vCH4 <- as.vector(vCH4)

  if(toupper(value)=='CH4') {
    return(vCH4)
  }

  if(tolower(value) %in% c('reactionn', 'reactionc')) {
    rxncoefs <- data.frame(substrate = -1, H2O = -cH2O, CH4 = cCH4, CO2 = cCO2, C5H7O2N = cbio, `NH4+` = cNH4, `HCO3-` = cHCO3)

    # If there is only one form value
    if(nrow(rxncoefs) == 1) names(rxncoefs)[1] <- form
    # Otherwise there is not a good way to organize output if there are different form values

    if(tolower(value) == 'reactionn') {
      return(rxncoefs)
    }

    reaction <- NULL

    # Otherwise make a text reaction (should move this to a separate function)
    for(i in 1:nrow(rxncoefs)) {
    
      # Change name from substrate (should switch to extraction of vector at top of loop instead)
      names(rxncoefs)[1] <- form[i]

      reactants <- - rxncoefs[i, rxncoefs[i,] < 0]
      reactants <- signif(reactants, 4)
      reactants[reactants == 1] <- ''
      reactants <- as.vector(reactants)
      reactants <- paste0(reactants, names(reactants))
      reactants <- paste(reactants, collapse = ' + ')
    
      products <- rxncoefs[i, rxncoefs[i,] > 0, drop = FALSE]
      products <- signif(products, 4)
      products[products == 1] <- ''
      products <- as.vector(products)
      products <- paste0(products, names(products))
      products <- paste(products, collapse = ' + ')
    
      rxntxt <- paste0(reactants, ' --> ', products)
      rxntxt <- gsub('NH4.', 'NH4+', rxntxt)
      rxntxt <- gsub('HCO3.', 'HCO3-', rxntxt)

      reaction <- c(reaction, rxntxt)
    
    }

    #names(reaction) <- form

    return(reaction)

  }

  # Hydrolytic water consumption, g H2O
  h <- cH2O*molMass('H2O')*fd*mol # NTS double check
  # Ammonia requirement, g N substrate (if < 0 then is produced not required)
  mNH4<- -cNH4*molMass('N')*fd*mol # NTS double check
  # Microbial biomass production
  mBio <- cbio*molMass('C5H7O2N')*fd*mol
  # CH4, CO2, and biogas in g
  mCH4 <- molMass('CH4')*nCH4
  nCO2 <- (cHCO3 + cCO2)*fd*mol
  mCO2 <- molMass('CO2')*nCO2
  # Biogas composition (per production, but does not consider TIC in solution)
  fCH4 <- nCH4/(nCH4 + nCO2)

  # Partition CO2 between headspace and solution if H2O mass, pH, and temperature are available
  if(!missing(conc.sub)[1] & !missing(pH)[1] & !missing(temp)[1] & all(nCO2 > 1E-8)) {
    # Note that method is hard-wired to 'cont' for now. 
    out <- t(mapply(partCO2, nCO2 = nCO2, nCH4 = nCH4, mass.H2O = mass/conc.sub, temp.c = temp, pH = pH, mu = mu, pres = 101325, method = 'cont', value = 'all', USE.NAMES = FALSE))
    nCO2Bg <- out[, 'nCO2Bg']
    nCO2.sol <- out[, 'nCO2.sol']
    cTIC <- out[, 'cTIC']
    xCH4 <- out[, 'xCH4']

    # Volume of CO2 left in biogas
    vCO2 <- vol.mol['CO2']*nCO2Bg
    # Dry biogas volume
    vBg <- vCH4 + vCO2

    # Masses of CO2 in solution and biogas
    mCO2Bg <- molMass('CO2')*nCO2Bg
    mCO2.sol <- molMass('CO2')*nCO2.sol

    if(all(fs==0)) {
      results <- data.frame(form = form, mass = mass, mol.mass = as.vector(mmass), moles = as.vector(mol), COD = COD, fd = fd, conc.sub = conc.sub, temp = temp, pH = pH, hydro = h, fCH4 = fCH4, xCH4 = xCH4, vCH4 = vCH4, vCO2 = vCO2, vBg = vBg, mCH4 = mCH4, mCO2 = mCO2, mCO2Bg = mCO2Bg, mCO2.sol = mCO2.sol, cTIC = cTIC)
    } else {
      results <- data.frame(form = form, mass = mass, mol.mass = as.vector(mmass), moles = as.vector(mol), COD = COD, fs = fs, fe = fe, fd = fd, conc.sub = conc.sub, temp = temp, pH = pH, hydro = h, fCH4 = fCH4, xCH4 = xCH4, vCH4 = vCH4, vCO2 = vCO2, vBg = vBg, mCH4 = mCH4, mCO2 = mCO2, mCO2Bg = mCO2Bg, mCO2.sol = mCO2.sol, cTIC = cTIC, m.bio = mBio, N.req = mNH4)
    }
    rownames(results) <- 1:nrow(results)
    return(results)
  }

  # Not enough info for CO2 partitioning
  if(all(fs==0)) {
    results <- data.frame(form = form, mass = mass, mol.mass = as.vector(mmass), moles = as.vector(mol), COD = COD, hydro = h, fCH4 = fCH4, vCH4 = vCH4, mCH4 = mCH4, mCO2 = mCO2)
  } else {
    results <- data.frame(form = form, mass = as.vector(mass), mol.mass = as.vector(mmass), moles = as.vector(mol), COD = as.vector(COD), fs = fs, fe = fe, fd = fd, hydro = h, fCH4 = fCH4, vCH4 = vCH4, mCH4 = mCH4, mCO2 = mCO2, m.bio = mBio, N.req = mNH4)
  }
  rownames(results) <- 1:nrow(results)
  return(results)
}
