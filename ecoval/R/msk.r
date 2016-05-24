msk.create <- function(language         = "English",
                       dictionaries     = NA,
                       col              = "black",
                       modify.nutrients = F)
{
  dict <- ecoval.dict(language,dictionaries)
  
  morphol <- msk.morphol.1998.create(language=language,dictionaries=dictionaries,col=col)
  hydrol  <- msk.hydrol.2011.create (language=language,dictionaries=dictionaries,col=col)
  physapp <- msk.physapp.2007.create(language=language,dictionaries=dictionaries,col=col)
  phys.branch <- utility.aggregation.create(name.node = ecoval.translate("N_phys",dict),
                                            nodes     = list(morphol,hydrol,physapp),
                                            name.fun  = "utility.aggregate.addmin",
                                            par       = c(0.4,0.4,0.2,0.5),
                                            names.par = c("w_morphol","w_hydrol","w_physapp",
                                                          "w_add_phys"),
                                            col       = col,
                                            required  = TRUE)
  
  nutrients <- msk.nutrients.2010.create(language=language,dictionaries=dictionaries,col=col,
                                         modify=modify.nutrients)
  chem.branch <- utility.aggregation.create(name.node = ecoval.translate("N_chem",dict),
                                            nodes     = list(nutrients),
                                            name.fun  = "utility.aggregate.addmin",
                                            par       = c(1,0.5),
                                            names.par = c("w_nutrients","w_add_chem"),
                                            col       = col,
                                            required  = TRUE)
  
  diatoms       <- msk.diatoms.2007.create      (language=language,dictionaries=dictionaries,col=col)
  invertebrates <- msk.invertebrates.2010.create(language=language,dictionaries=dictionaries,col=col)
  fish          <- msk.fish.2004.create         (language=language,dictionaries=dictionaries,col=col)
  biol.branch <- utility.aggregation.create(name.node = ecoval.translate("N_biol",dict),
                                            nodes     = list(diatoms,invertebrates,fish),
                                            name.fun  = "utility.aggregate.addmin",
                                            par       = c(1,1,1,0.5),
                                            names.par = c("w_diatoms","w_invertebrates","w_fish",
                                                          "w_add_biol"),
                                            col       = col,
                                            required  = TRUE)
  
  ecol <- utility.aggregation.create(name.node = ecoval.translate("N_ecol",dict),
                                     nodes     = list(phys.branch,chem.branch,biol.branch),
                                     name.fun  = "utility.aggregate.addmin",
                                     par       = c(1,1,1,0.5),
                                     names.par = c("w_phys","w_chem","w_biol","w_add_ecol"),
                                     col       = col)
  return(ecol)
}


