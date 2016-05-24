rbeauti  <- function(s, file, taxonset){
	
  # prepare alignment
  # ------------
  ss  <-  as.list(s)
  ss <- lapply(ss, function(x)  paste(as.character(x), collapse = ""))
  ss <- lapply(ss, gsub, pattern = "n", replacement = "?")
  ss <- data.frame(id = paste("seq", names(ss), sep = "_"),
                   taxon = names(ss),
                   value = unlist(ss),
                   stringsAsFactors = FALSE)
  id <- tolower(substr(ss$taxon[1], 1, 3))
  sequence <- apply(ss, 1, function(x) xmlNode("sequence", attrs = c(x[1], x[2], x[3])))
  data <- xmlNode("data", 
                  attrs = c(id = id, name = "alignment"),
                  .children = sequence)
  
 
  ## assemble node <state>
  ## ---------------------
  N.taxonset <- xmlNode("taxonset", 
                      attrs = c(id = paste("TaxonSet.", id, sep = ""), spec = "TaxonSet"),
                      .children = list(xmlNode("data",
                                               attrs = c(idref = id,
                                                         name = "alignment"))))
  tree <- xmlNode("tree", 
                  attrs = c(id = paste("Tree.t:", id, sep = ""),
                                    name = "stateNode"),
                  .children = list(N.taxonset))
  state <- xmlNode("state", 
                   attrs = c(id = "state",
                            storeEvery = "5000"),
                   .children = list(tree,
                                    xmlNode("parameter", "1.0",
                                            attrs = c(id = paste("birthRate.t:", id, sep = ""),
                                                      name = "stateNode"))))
  
  ## assemble node <init>
  ## ---------------------
  parameter <- xmlNode("parameter", "1.0", 
                       attrs = c(id = paste("randomPopSize.t:", id, sep = ""),
                                 name = "popSize"))
  populationModel <- xmlNode("populationModel",
                            attrs = c(id = paste("ConstantPopulation0.t:", id, sep = ""),
                                      spec = "ConstantPopulation"),
                            .children = list(parameter))
  init <- xmlNode("init", 
                  attrs = c(estimate = "false",
                                    id = paste("RandomTree.t:", id, sep = ""),
                                    initial = paste("@Tree.t:", id, sep = ""),
                                    spec = "beast.evolution.tree.RandomTree",
                                    taxa = paste("@", id, sep = "")),
                  .children = list(populationModel))
  
  ## assemble node <distribution>
  ## ----------------------------
  distribution <- xmlNode("distribution", 
                          attrs = c(birthDiffRate = paste("@birthRate.t:", id, sep = ""),
                                    id = paste("YuleModel.t:", id, sep = ""),
                                    spec = "beast.evolution.speciation.YuleModel",
                                    tree = paste("@Tree.t:", id, sep = "")))
  prior <- xmlNode("prior",
                   attrs = c(id = paste("YuleBirthRatePrior.t:", id, sep = ""),
                             name = "distribution",
                             x = paste("@birthRate.t:", id, sep = "")),
                   .children = list(xmlNode("Uniform",
                                            attrs = c(id = "Uniform.0",
                                                      name ="distr",
                                                      upper = "Infinity"))))
  prior.distribution <- xmlNode("distribution", distribution, prior,
                                attrs = c(id = "prior",
                                          spec = "util.CompoundDistribution"))
  if ( !missing(taxonset) ) {
    ts <- lapply(taxonset, rbeauti.taxonset, id = id)
    prior.distribution <- addChildren(prior.distribution, kids = ts)
  }
  
  siteModel <- xmlNode("siteModel", 
                       attrs = c(id = paste("SiteModel.s:", id, sep = ""),
                                 spec = "SiteModel"),
                       .children = list(xmlNode("parameter",  "1.0",
                                                attrs = c(estimate = "false",
                                                          id = paste("mutationRate.s:", id, sep = ""),
                                                          name = "mutationRate")),
                                        xmlNode("parameter", "1.0",
                                                attrs = c(estimate = "false",
                                                          id = paste("gammaShape.s:", id, sep = ""),
                                                          name = "shape")),
                                        xmlNode("parameter", "0.0", 
                                                attrs = c(estimate = "false",
                                                          id = paste("proportionInvariant.s:", id, sep = ""),
                                                          lower = "0.0",
                                                          name = "proportionInvariant",
                                                          upper = "1.0")),
                                        xmlNode("substModel", attrs = c(id = paste("JC69.s:", id, sep = ""),
                                                                        spec="JukesCantor"))))
  branchRateModel <- xmlNode("branchRateModel",
                             attrs = c(id = paste("StrictClock.c:", id, sep = ""),
                                       spec = "beast.evolution.branchratemodel.StrictClockModel"),
                             .children = list(xmlNode("parameter", "1.0",
                                                      attrs = c(estimate = "false",
                                                                id = paste("clockRate.c:", id, sep = ""),
                                                                name = "clock.rate"))))
  distribution <- xmlNode("distribution",
                          attrs = c(data = paste("@", id, sep = ""),
                                    id = paste("treeLikelihood.", id, sep = ""),
                                    spec = "TreeLikelihood",
                                    tree = paste("@Tree.t:", id, sep = "")),
                          .children = list(siteModel, branchRateModel))
  likelihood.distribution <- xmlNode("distribution", distribution,
                                     attrs = c(id = "likelihood",
                                               spec = "util.CompoundDistribution"))
  distribution <- xmlNode("distribution", 
                          attrs = c(id = "posterior",
                                    spec = "util.CompoundDistribution"),
                          .children = list(prior.distribution, likelihood.distribution))
  
  ## assemble node <logger> 1
  ## ------------------------
  logger1 <- xmlNode("logger", attrs = c(fileName = paste(id, ".log", sep = ""),
                                  id = "tracelog",
                                  logEvery = "1000",
                                  model = "@posterior",
                                  sanitiseHeaders = "true",
                                  sort="smart"))
  logger1 <- addChildren(logger1, kids = list(
    xmlNode("log", attrs = c(idref = "posterior")),
    xmlNode("log", attrs = c(idref = "likelihood")),
    xmlNode("log", attrs = c(idref = "prior")),
    xmlNode("log", attrs = c(idref = paste("treeLikelihood.", id, sep = ""))),
    xmlNode("log", attrs = c(id = paste("TreeHeight.t:", id, sep = ""),
                             spec = "beast.evolution.tree.TreeHeightLogger",
                             tree = paste("@Tree.t:", id, sep = ""))),
    xmlNode("log", attrs = c(idref = paste("YuleModel.t:", id, sep = ""))),
    xmlNode("parameter", attrs = c(idref = paste("birthRate.t:", id, sep = ""), 
                                   name = "log"))
    ))
  if ( !missing(taxonset) ){
    tslog <- lapply(taxonset, function(x) xmlNode("log", attrs = c(idref = paste(x$id, ".prior", sep = ""))))
    logger1 <- addChildren(logger1, kids = tslog)
  }
  
  ## assemble node <logger> 2
  ## ------------------------
  logger2 <- xmlNode("logger", attrs = c(id = "screenlog", logEvery = "1000"))
  logger2 <- addChildren(logger2, kids = list(
    xmlNode("log", attrs = c(idref = "posterior")),
    xmlNode("log", attrs = c(arg = "@posterior",
                             id = "ESS.0",
                             spec = "util.ESS")),
    xmlNode("log", attrs = c(idref = "likelihood")),
    xmlNode("log", attrs = c(idref = "prior"))
    ))
  
  ## assemble node <logger> 3
  ## ------------------------
  logger3 <- xmlNode("logger", attrs = c(fileName = "$(tree).trees",
                                  id = paste("treelog.t:", id, sep = ""),
                                  logEvery = "1000",
                                  mode = "tree"))
  logger3 <- addChildren(logger3, kids = list(
    xmlNode("log", attrs = c(id = paste("TreeWithMetaDataLogger.t:", id, sep = ""),
                           spec = "beast.evolution.tree.TreeWithMetaDataLogger",
                           tree = paste("@Tree.t:", id, sep = "")))
  ))
  
  ## assemble node <run> 
  ## ---------------------
  run <- xmlNode("run", attrs = c(chainLength = "10000000", 
                                  id = "mcmc",
                                  spec = "MCMC"))
  run <- addChildren(run, kids = list(state,
                                      init,
                                      distribution,
                                      xmlNode("operator", attrs = c(id = paste("YuleBirthRateScaler.t:", id, sep = ""),
                                                                    parameter = paste("@birthRate.t:", id, sep = ""),
                                                                    scaleFactor = "0.75",
                                                                    spec = "ScaleOperator",
                                                                    weight = "3.0")),
                                      xmlNode("operator", attrs = c(id = paste("treeScaler.t:", id, sep = ""),
                                                                    scaleFactor = "0.5",
                                                                    spec = "ScaleOperator",
                                                                    tree = paste("@Tree.t:", id, sep = ""),
                                                                    weight = "3.0")),
                                      xmlNode("operator", attrs = c(id = paste("treeRootScaler.t:", id, sep = ""),
                                                                    rootOnly = "true",
                                                                    scaleFactor = "0.5",
                                                                    spec = "ScaleOperator",
                                                                    tree = paste("@Tree.t:", id, sep = ""),
                                                                    weight = "3.0")),
                                      xmlNode("operator", attrs = c(id = paste("UniformOperator.t:", id, sep = ""), 
                                                                    spec = "Uniform",
                                                                    tree = paste("@Tree.t:", id, sep = ""),
                                                                    weight = "30.0")),
                                      xmlNode("operator", attrs = c(id = paste("SubtreeSlide.t:", id, sep = ""),
                                                                    spec = "SubtreeSlide",
                                                                    tree = paste("@Tree.t:", id, sep = ""),
                                                                    weight = "15.0")),
                                      xmlNode("operator", attrs = c(id = paste("narrow.t:", id, sep = ""),
                                                                    spec = "Exchange",
                                                                    tree = paste("@Tree.t:", id, sep = ""),
                                                                    weight = "15.0")),
                                      xmlNode("operator", attrs = c(id = paste("wide.t:", id, sep = ""),
                                                                    isNarrow = "false",
                                                                    spec = "Exchange",
                                                                    tree = paste("@Tree.t:", id, sep = ""),
                                                                    weight = "3.0")),
                                      xmlNode("operator", attrs = c(id = paste("WilsonBalding.t:", id, sep = ""),
                                                                    spec = "WilsonBalding",
                                                                    tree = paste("@Tree.t:", id, sep = ""),
                                                                    weight = "3.0")),
                                      logger1,
                                      logger2,
                                      logger3))
  
  ## assemble nodes <map> 
  ## --------------------
  m <- mm <- c("Beta", "Exponential", "InverseGamma", "LogNormalDistributionModel", "Gamma", 
            "Uniform", "Prior", "LaplaceDistribution", "OneOnX", "Normal")
  m[m == "Prior"] <- "prior"; m[m == "LogNormalDistributionModel"] <- "LogNormal"
  mm <- paste("beast.math.distributions", mm, sep = ".")
  m <- cbind(name = m, mm)
  maps <- apply(m, 1, function(x) xmlNode("map", x[2], 
                                     attrs = x[1]))
  
  
  ## assemble node <beast> 
  ## ---------------------
  namespace <- c("beast.core", "beast.evolution.alignment", "beast.evolution.tree.coalescent", 
                 "beast.core.util", "beast.evolution.nuc", "beast.evolution.operators", 
                 "beast.evolution.sitemodel", "beast.evolution.substitutionmodel", "beast.evolution.likelihood")
  namespace <- paste(namespace, collapse = ":")
  beast <- xmlNode("beast", attrs = c(beautitemplate = "Standard",
                                      beautistatus = "",
                                      namespace = namespace,
                                      version = "2.0"))
  beast <- addChildren(beast, kids = list(data))
  beast <- addChildren(beast, kids = maps)
  beast <- addChildren(beast, kids = list(run))
  
  if ( missing(file) ) return(beast)
  else {
    saveXML(beast, file = file, encoding = "UTF-8",
            prefix = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>')
  }
}