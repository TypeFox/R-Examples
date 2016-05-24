ecoval.river.create <- function(phys         = list(msk.morphol.1998.create,
                                                    msk.physapp.2007.create),
                                physagg      = "utility.aggregate.addmin",
                                physpar      = numeric(0),
                                chem         = list(msk.nutrients.2010.create,
                                                    val.micropoll.create,
                                                    val.heavymetals.create),
                                chemagg      = "utility.aggregate.addmin",
                                chempar      = numeric(0),
                                biol         = list(msk.diatoms.2007.create,
                                                    val.invertebrates.create,
                                                    msk.fish.2004.create),
                                biolagg      = "utility.aggregate.addmin",
                                biolpar      = numeric(0),
                                ecolagg      = "utility.aggregate.addmin",
                                ecolpar      = numeric(0),
                                language     = "English",
                                dictionaries = NA,
                                col          = "black")
{
  dict <- ecoval.dict(language,dictionaries)

  phys.nodes  <- list()
  phys.branch <- NA
  if ( length(phys) > 0 )
  {
    for ( i in 1:length(phys) )
    {
      if ( is.function(phys[[i]]) )
      {
        phys.nodes[[i]] <- phys[[i]](language=language,dictionaries=dictionaries,col=col)
      }
      else
      {
        phys.nodes[[i]] <- phys[[i]]
      }
    }
    par <- c(rep(1,length(phys.nodes)),0.5)
    if ( length(physpar) > 0 ) par <- physpar
    phys.branch <- utility.aggregation.create(name.node = ecoval.translate("N_phys",dict),
                                              nodes     = phys.nodes,
                                              name.fun  = physagg,
                                              par       = par,
                                              col       = col,
                                              required  = TRUE)
  }

  chem.nodes  <- list()
  chem.branch <- NA
  if ( length(chem) > 0 )
  {
    for ( i in 1:length(chem) )
    {
      if ( is.function(chem[[i]]) )
      {
        chem.nodes[[i]] <- chem[[i]](language=language,dictionaries=dictionaries,col=col)
      }
      else
      {
        chem.nodes[[i]] <- chem[[i]]
      }
    }
    par <- c(rep(1,length(chem.nodes)),0.5)
    if ( length(chempar) > 0 ) par <- chempar
    chem.branch <- utility.aggregation.create(name.node = ecoval.translate("N_chem",dict),
                                              nodes     = chem.nodes,
                                              name.fun  = chemagg,
                                              par       = par,
                                              col       = col,
                                              required  = TRUE)
  }
  
  biol.nodes  <- list()
  biol.branch <- NA
  if ( length(biol) > 0 )
  {
    for ( i in 1:length(biol) )
    {
      if ( is.function(biol[[i]]) )
      {
        biol.nodes[[i]] <- biol[[i]](language=language,dictionaries=dictionaries,col=col)
      }
      else
      {
        biol.nodes[[i]] <- biol[[i]]
      }
    }
    par <- c(rep(1,length(biol.nodes)),0.5)
    if ( length(biolpar) > 0 ) par <- biolpar
    biol.branch <- utility.aggregation.create(name.node = ecoval.translate("N_biol",dict),
                                              nodes     = biol.nodes,
                                              name.fun  = biolagg,
                                              par       = par,
                                              col       = col,
                                              required  = TRUE)
  }
    
  ecol.nodes <- list()
  i <- 1
  if ( !is.na(phys.branch[[1]]) ) { ecol.nodes[[i]] <- phys.branch; i <- i+1 }
  if ( !is.na(chem.branch[[1]]) ) { ecol.nodes[[i]] <- chem.branch; i <- i+1 }
  if ( !is.na(biol.branch[[1]]) ) { ecol.nodes[[i]] <- biol.branch; i <- i+1 }
  if ( i == 1 ) return(NA)
  par <- c(rep(1,length(ecol.nodes)),0.5)
  if ( length(ecolpar) > 0 ) par <- ecolpar
  ecol <- utility.aggregation.create(name.node = ecoval.translate("N_ecol",dict),
                                     nodes     = ecol.nodes,
                                     name.fun  = ecolagg,
                                     par       = par,
                                     col       = col)
   return(ecol)
}


