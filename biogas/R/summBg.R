# Modified: 14 April 2016 SDH

summBg <- function(
  vol,
  setup,
  id.name = 'id',
  time.name = 'time',
  descrip.name = 'descrip',
  inoc.name = NULL,
  norm.name = NULL,
  inoc.m.name = 'minoc',
  vol.name = 'cvCH4',
  imethod = 'linear',
  extrap = FALSE,
  when = 30,
  show.obs = FALSE, 
  sort = TRUE) 
{

  # Argument checks~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  checkArgClassValue(vol, 'data.frame')
  checkArgClassValue(setup, 'data.frame')
  checkArgClassValue(id.name, 'character')
  checkArgClassValue(time.name, c('character', 'NULL'))
  checkArgClassValue(descrip.name, c('character', 'NULL'))
  checkArgClassValue(inoc.name, c('character', 'NULL'))
  checkArgClassValue(norm.name, c('character', 'NULL'))
  checkArgClassValue(inoc.m.name, c('character', 'NULL'))
  checkArgClassValue(vol.name, 'character')
  # Skip imethod, since it is checked in interp()
  checkArgClassValue(extrap, 'logical')
  checkArgClassValue(when, c('numeric', 'integer', 'character', 'NULL'))
  checkArgClassValue(show.obs, 'logical')
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Echo response variable
  message('Response variable (volume) is ', deparse(substitute(vol)), '$', vol.name, '.')

  # Check for missing columns in vol
  if(class(when) %in% c('numeric', 'integer')) {
    if(any(missing.col <- !c(id.name, time.name, vol.name) %in% names(vol))){
      stop('Specified columns in vol data frame (', deparse(substitute(vol)), ') not found: ', c(id.name, time.name, vol.name)[missing.col], '.')
    } 
  } else { # when is 'end' or 'meas'
    if(any(missing.col <- !c(id.name, vol.name) %in% names(vol))){
      stop('Specified columns in vol data frame (', deparse(substitute(vol)), ') not found: ', c(id.name, vol.name)[missing.col], '.')
    } 
  }

  # Check for missing columns in setup
  if(any(missing.col <- !c(id.name, descrip.name) %in% names(setup))){
    stop('Specified columns in setup data frame (', deparse(substitute(setup)), ') not found: ', c(id.name, descrip.name)[missing.col], '.')
  } 

  # Check that inoc.name and norm.name can be found in setup data frame
  if(!is.null(inoc.name) && !inoc.name %in% setup[, descrip.name]) {
    stop('inoc.name ', deparse(substitute(inoc.name)), ' not found in ', deparse(substitute(setup)), '$', descrip.name, '.')
  }

  if(!is.null(norm.name) && !norm.name %in% names(setup)) {
    stop('norm.name ', deparse(substitute(norm.name)), ' not found in the column names of ', deparse(substitute(setup)), '.')
  }

  # Add other checks here

  # Trim setup based on ids and check again for inoc.name and norm.name~~~~~~~~~~~~~~~~~~~
  # Find reactor/bottle IDs present in both vol and setup
  ids <- intersect(setup[, id.name], vol[, id.name])

  setup <- setup[setup[, id.name] %in% ids, ]

  if(!is.null(inoc.name) && !inoc.name %in% setup[, descrip.name]) {
    stop('inoc.name ', deparse(substitute(inoc.name)), ' no longer in setup after trimming--are reactors present in setup missing in vol?')
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Remove inoc ids
  if(!is.null(inoc.name)) {
    ids <- setup[setup[, descrip.name]!=inoc.name, id.name]
    ids.inoc <- setup[setup[, descrip.name]==inoc.name, id.name]
  }

  # Check for duplicates in setup and vol~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(any(duplicated(setup[, id.name]))) {
    stop('Duplicated reactor IDs (', id.name, ' column) in setup dataframe! This must be an error.')
  }

  if(any(duplicated(vol[, c(id.name, time.name)]))) {
    stop('Duplicated ID (', id.name, ' column) x time (', time.name, ' column) in vol dataframe! This must be an error.')
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Interpolate cvCH4 to common time for each reactor~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Or select values for analysis (when = 'end' or 'meas')

  if(class(when) %in% c('numeric', 'integer')) {
    summ1 <- expand.grid(id = ids, time = when)
    names(summ1) <- c(id.name, time.name)

    # Then interpolate
    for(i in ids) {
      dc <- vol[vol[, id.name]==i, ]
      # Interpolate if more than one value is present

      if(nrow(dc)>1) {
        summ1[summ1[, id.name]==i, vol.name] <- interp(dc[, time.name], dc[, vol.name], time.out = when, method = imethod, extrap = extrap)
      } else {
	if(dc[, time.name]==when) { # `when` argument matches the single time present
          summ1[summ1[, id.name]==i, vol.name] <- dc[, vol.name]
	} else {
          summ1[summ1[, id.name]==i, vol.name] <- NA
      	  warning('There is only a single ', vol.name, ' value for reactor ', i,', and it does not match the specified when (', when, '). Interpolation is not possible.')
	}
      }

    }

  } else if(length(when) == 1 && tolower(when) == 'end') { # User just wants to use latest values of volume

    summ1 <- data.frame(id = ids, time = NA, vol = NA)
    names(summ1) <- c(id.name, time.name, vol.name)

    # Sort, in order to find latest values
    vol <- vol[order(vol[, id.name], vol[, vol.name]), ]

    for(i in ids) {
      dc <- vol[vol[, id.name]==i, ]
      # Select the last row from sorted data frame
      summ1[summ1[, id.name]==i, c(time.name, vol.name)] <- dc[nrow(dc), c(time.name, vol.name)]
    }

  } else if(length(when) == 1 && when == 'meas') { # Return values for all measurement times, which may differ among reactors

    summ1 <- vol[vol[, id.name] %in% ids, c(id.name, time.name, vol.name)]

  } else {

    stop('when argument not recognized. Options are numeric or integer vector, \"end\" or \"meas\".')

  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Get all unique times
  times.summ <- unique(summ1[, time.name])

  # Work with inoculum data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Now interpolate inoculum-only reactors to all unique times
  if(!is.null(inoc.name)) {

    summ.inoc <- expand.grid(id = ids.inoc, time = times.summ)

    # Then interpolate inoculum production (each inoc reactor) to each time
    for(i in ids.inoc) {

      dc <- vol[vol[, id.name]==i, ]

      # Interpolate if more than one value is present
      if(nrow(dc)>1) {
        summ.inoc[summ.inoc$id==i, vol.name] <- interp(dc[, time.name], dc[, vol.name], time.out = times.summ, method = imethod, extrap = extrap)
      } else {

	if(dc[, time.name]==times.summ) { # `when` argument matches the single time present
          summ.inoc[summ.inoc$id==i, vol.name] <- dc[, vol.name]
	} else {
          summ.inoc[summ.inoc$id==i, vol.name] <- NA
      	  warning('There is only a single ', vol.name, ' value for reactor ', i,', and it does not match the specified when (', when, '). Interpolation is not possible.')
	}

      }

    }

    # Merge to add mass inoculum and VS in substrate
    summ.inoc <- merge(setup, summ.inoc, by.x = id.name, by.y = 'id')

    # Volume contribution per unit inoculum mass
    summ.inoc$vol.mi <- summ.inoc[, vol.name]/summ.inoc[, inoc.m.name]

    # Mean and sd volume contribution per unit inoc mass
    inoc.vol <- data.frame()

    for(i in times.summ) {
      vol.mi <- summ.inoc[summ.inoc$time == i, 'vol.mi']
      inoc.vol <- rbind(inoc.vol, c(time = i, mn = mean(vol.mi), s = sd(vol.mi)))
    }

    names(inoc.vol) <- c(time.name, 'vol.mi.mn', 'vol.mi.sd')
    # inoc.vol has mean and sd vol per unit mass inoc for all times

  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Samples
  # Add mass of inoculum and VS in substrate
  summ1 <- merge(setup, summ1, by = id.name)

  if(!is.null(inoc.name)) {

    # Merge inoculum normalized volumes with sample data
    summ1 <- merge(summ1, inoc.vol, by = time.name)

    # Calculate and substract inoc contribution
    # Next three lines only for returning additional info when show.obs = TRUE
    summ1[, paste0(vol.name, '.tot')] <- summ1[, vol.name]
    summ1[, paste0(vol.name, '.inoc')] <- summ1$vol.mi.mn*summ1[, inoc.m.name]
    summ1[, 'fv.inoc'] <- summ1[, paste0(vol.name, '.inoc')]/summ1[, paste0(vol.name, '.tot')]

    # Correct vol for inoculum
    summ1[, vol.name] <- summ1[, vol.name] - summ1$vol.mi.mn*summ1[, inoc.m.name]

    # Add sd in volume produced by inoculum for use below in error propagation
    summ1[, 'sd.inoc'] <- summ1$vol.mi.sd*summ1[, inoc.m.name]

  } else {

    # NTS: How did I handle this before 10 Feb 2016?
    summ1[, 'sd.inoc'] <- 0

  }

  # Normalization
  if(!is.null(norm.name)) { 
    # Normalize remaining vol by norm.name (typically by substrate VS)
    summ1[, vol.name] <- summ1[, vol.name]/summ1[, norm.name]

    # Normalize sd contribution from inoc by the same value
    summ1[, 'sd.inoc'] <- summ1[, 'sd.inoc']/summ1[, norm.name]

    # Next two lines only for returning additional info when show.obs = TRUE
    # Only have the .tot and .inoc columns when inoc is subtracted out
    if(!is.null(inoc.name) && inoc.name %in% setup[, descrip.name]) { 
      summ1[, paste0(vol.name, '.tot')] <- summ1[, paste0(vol.name, '.tot')]/summ1[, norm.name]
      summ1[, paste0(vol.name, '.inoc')] <- summ1[, paste0(vol.name, '.inoc')]/summ1[, norm.name]
    }
  } 

    # Calculate means and sd for a summary
  if(!show.obs) {
    # Summarize by description
    summ2 <- unique(summ1[, c(time.name, descrip.name)]) # NTS: may want to put time second

    for(i in unique(summ1[, descrip.name])){
      dd <- summ1[summ1[, descrip.name]==i, ]
      for(j in unique(dd[, time.name])) {
        ddd <- dd[dd[, time.name]==j, ]
        summ2[summ2[, descrip.name]==i & summ2[, time.name]==j, 'mean'] <- mean(na.omit(ddd[, vol.name]))
        #summ2[summ2[, descrip.name]==i, 'sd'] <- signif(sd(na.omit(ddd[, vol.name])), 5) # Original code with no contribution of inoc to sd
        summ2[summ2[, descrip.name]==i & summ2[, time.name]==j, 'sd'] <- sqrt(sd(na.omit(ddd[, vol.name]))^2 + mean(ddd[, 'sd.inoc'])^2) # NTS: Probably there is a better way to do this
        summ2[summ2[, descrip.name]==i & summ2[, time.name]==j, 'n'] <- sum(!is.na(ddd[, vol.name]))  
      }
    }
  } else { # If show.obs = TRUE, just return individual observations
    #summ1$sd.inoc <- NULL
    summ2 <- summ1[order(summ1[, descrip.name], summ1[, id.name], summ1[, time.name]), ]
  }


  # More messages~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Messages about inoculum 
  if(!is.null(inoc.name) && inoc.name %in% setup[, descrip.name]) { # Inoculum contribution subtracted
    #message('Inoculum contribution subtracted based on ', deparse(substitute(setup.orig)), '$', inoc.m.name, '.') 
    message('Inoculum contribution subtracted based on setup$', inoc.m.name, '.') 
  } else {
      message('Inoculum contribution not subtracted.') 
  }

  # Message about normalization
  if(!is.null(norm.name)) { 
    #message('Response normalized by ', deparse(substitute(setup)), '$', norm.name, '.')
    message('Response normalized by setup$', norm.name, '.')
  } else {
    message('No normalization by substrate mass.')
  }

  # Select columns
  if(!show.obs) {
    summ2 <- summ2[ , c(descrip.name, time.name, 'mean', 'sd', 'n')]
  } 

  # Sort result
  if(sort) {
    summ2 <- summ2[order(summ2[, descrip.name], summ2[, time.name]), ]
  } else {
    # Get original reactor order from setup
    descrip.order <- 1:length(unique(setup[, descrip.name]))
    names(descrip.order) <- setup[!duplicated(setup[, descrip.name]), descrip.name]

    # Sort
    summ2 <- summ2[order(descrip.order[as.character(summ2[, descrip.name])], summ2[, time.name]), ]
  }

  # Row names
  rownames(summ2) <- 1:nrow(summ2)

  return(summ2)

}
