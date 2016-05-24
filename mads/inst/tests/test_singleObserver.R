library(mrds)
library(testthat)

context("Single Observer Analyses")

test_that("Test Analyses", {


  #datasetup
  ex.filename<-system.file("testData/input_checks/ddf_dat.robj", package="mads")
  load(ex.filename)
  ex.filename<-system.file("testData/input_checks/obs_table.robj", package="mads")
  load(ex.filename)
  ex.filename<-system.file("testData/input_checks/region_table.robj", package="mads")
  load(ex.filename)
  ex.filename<-system.file("testData/input_checks/sample_table.robj", package="mads")
  load(ex.filename)

  #run ddf analyses
  ddf.1 <- ddf(dsmodel = ~mcds(key = "hn", formula = ~ size), method='ds', data=ddf.dat,meta.data=list(width=4))
  ddf.2 <- ddf(dsmodel = ~mcds(key = "hr", formula = ~ size), method='ds', data=ddf.dat,meta.data=list(width=4))
  ddf.3 <- ddf(dsmodel = ~mcds(key = "hn", formula = ~ 1, adj.series = "cos", adj.order = c(2)), method='ds', data=ddf.dat,meta.data=list(width=4, mono=TRUE))
  #think this should have been fixed in mrds
  ddf.1$data$detected <- rep(1, nrow(ddf.1$data))
  ddf.2$data$detected <- rep(1, nrow(ddf.2$data))
  ddf.3$data$detected <- rep(1, nrow(ddf.3$data))

  #Multi-analysis options
  model.names              <- list("CD"=c("ddf.1","ddf.2","ddf.3"), "WD"=c("ddf.1","ddf.2","ddf.3"), "UnidDol"=c("ddf.1","ddf.2","ddf.3"))
  ddf.models               <- list("ddf.1" = ddf.1, "ddf.2" = ddf.2, "ddf.3" = ddf.3)
  species.code.definitions <- list("UnidDol" = c("CD","WD"))
  species.presence         <- list("A" = c("CD","WD"))
  covariate.uncertainty    <- NULL
  ddf.model.options        <- list(criterion="AIC")
  ddf.model.options$distance.naming.conv <- TRUE
  bootstrap                <- TRUE
  bootstrap.options        <- list(resample="samples", n=2, quantile.type = 7)
  dht.options              <- list(convert.units = 1)

  set.seed(747)
  results.to.compare <- execute.multi.analysis(
                              species.code = names(model.names),
                              unidentified.sightings = species.code.definitions,
                              species.presence = species.presence,
                              covariate.uncertainty = covariate.uncertainty,
                              models.by.species.code = model.names,
                              ddf.model.objects = ddf.models,
                              ddf.model.options = ddf.model.options,
                              region.table = region.table,
                              sample.table = sample.table,
                              obs.table = obs.table,
                              bootstrap = bootstrap,
                              bootstrap.option = bootstrap.options,
                              silent = FALSE)

  set.seed(747)
  MAE.warnings    <- NULL
  species.code    <- names(model.names)
  ddf.model.info  <- check.ddf.models(model.names, ddf.models)
  clusters        <- ddf.model.info$clusters
  double.observer <- ddf.model.info$double.observer
  # If the user has not specified the criteria set it
  if(is.null(ddf.model.options$criterion)){
    ddf.model.options$criterion <- "AIC"
  }
  # If the user has not specified the species field name set it
  if(is.null(ddf.model.options$species.field.name)){
    ddf.model.options$species.field.name <- "species"
  }

  ##################################
  expect_that(clusters, is_true())
  expect_that(double.observer, is_false())
  ##################################

  species.code.definitions <- check.species.code.definitions(species.code.definitions, species.code)
  unidentified.species     <- species.code.definitions$unidentified
  species.code.definitions <- species.code.definitions$species.code.definitions

  ##################################
  expect_that(unidentified.species, is_true())
  ##################################

  species.presence         <- check.species.presence(species.presence, species.code, strata.name = as.character(region.table$Region.Label))

  ##################################
  expect_that(names(species.presence), is_identical_to(c("A")))
  expect_that(species.presence[[1]], is_identical_to(c("CD","WD")))
  ##################################

  species.presence.compare <- species.presence
  species.presence         <- NULL
  species.presence         <- check.species.presence(species.presence, species.code, strata.name = as.character(region.table$Region.Label))

  ##################################
  #expect_that(species.presence, is_identical_to(species.presence.compare))
  rm(species.presence.compare)
  ##################################

  covariate.uncertainty    <- check.covar.uncertainty(covariate.uncertainty)
  check.bootstrap.options(bootstrap, bootstrap.options$resample, bootstrap.options$n, sample.table)
  bootstrap.options$n <- ifelse(bootstrap, bootstrap.options$n, 1)
  #Make master copies of all the datasets
  ddf.dat.master      <- get.datasets(model.names, ddf.models)
  unique.model.names  <- ddf.dat.master$unique.model.names
  model.index         <- ddf.dat.master$model.index
  ddf.dat.master      <- ddf.dat.master$ddf.dat.master

  ##################################
  expect_that(unique.model.names, is_identical_to(list("CD" = c("ddf.1", "ddf.2", "ddf.3"))))
  test <- c("CD","CD","CD")
  names(test) <- c("CD","WD","UnidDol")
  expect_that(model.index, is_identical_to(test))
  rm(test)
  expect_that(length(ddf.dat.master), equals(1))
  expect_that(nrow(ddf.dat.master[[1]]), equals(nrow(ddf.1$data)))
  ##################################

  obs.table.master    <- obs.table
  sample.table.master <- sample.table
  #Create storage for results (only for the species codes not the unidentified codes)
  bootstrap.results <- create.result.arrays(species.code, species.code.definitions, region.table, clusters, bootstrap.options$n)
  bootstrap.ddf.statistics <- create.param.arrays(unique.model.names, ddf.models, bootstrap.options$n, ddf.model.options$criterion)

  ##################################
  expect_that(names(bootstrap.ddf.statistics), matches("CD"))
  expect_that(dimnames(bootstrap.results$individual.summary)[[4]], is_identical_to(c("CD","WD")))
  ##################################

  n=1
  #Resample Data
  bootstrap = TRUE
  if(bootstrap){
    ddf.dat.working <- resample.data(resample=bootstrap.options$resample, obs.table.master, sample.table.master, ddf.dat.master, double.observer)
    obs.table       <- ddf.dat.working$obs.table
    sample.table    <- ddf.dat.working$sample.table
    ddf.dat.working <- ddf.dat.working$ddf.dat.working
  }else{
    ddf.dat.working <- ddf.dat.master
  }

  ##################################
  expect_that(length(unique(sample.table$Sample.Label)), equals(length(unique(sample.table.master$Sample.Label))))
  expect_that(table(sample.table$Region), is_identical_to(table(sample.table.master$Region)))
  expect_that(nrow(ddf.dat.working[[1]]), equals(nrow(obs.table)))
  expect_that(length(which(ddf.dat.working[[1]]$object%in%obs.table$object)), equals(nrow(obs.table)))
  expect_that(ddf.dat.working[["CD"]]$distance[ddf.dat.working[["CD"]]$object == 16], equals(ddf.dat.master[["CD"]]$distance[ddf.dat.master[["CD"]]$object == 16]))
  ##################################

  #ddf.dat.working.check <- ddf.dat.working
  if(!is.null(covariate.uncertainty)){
    ddf.dat.working <- resample.covariates(ddf.dat.working, covariate.uncertainty, MAE.warnings)
    MAE.warnings <- ddf.dat.working$MAE.warnings
    ddf.dat.working <- ddf.dat.working$ddf.dat.working
  }

  ##################################
  #expect_that(ddf.dat.working[["10"]]$object, is_identical_to(ddf.dat.working.check[["10"]]$object))
  #expect_that(ddf.dat.working[["10"]]$scaledtotsize[1] == ddf.dat.working.check[["10"]]$scaledtotsize[1], is_false())
  #expect_that(ddf.dat.working[["10"]]$distance[ddf.dat.working[["10"]]$object == 106], equals(ddf.dat.master[["10"]]$distance[ddf.dat.master[["10"]]$object == 106]))
  ##################################

  #Fit ddf models to all species codes
  ddf.results <- fit.ddf.models(ddf.dat.working, unique.model.names, ddf.models, ddf.model.options$criterion, bootstrap.ddf.statistics, n, MAE.warnings)
  if(class(ddf.results) == "list"){
    bootstrap.ddf.statistics <- ddf.results$bootstrap.ddf.statistics
    ddf.results <- ddf.results$ddf.results
  }else{
    #If the ddf results are not valid for all species move to next bootstrap iteration
    MAE.warnings <- ddf.results
    next
  }

  ##################################
  expect_that(as.numeric(bootstrap.ddf.statistics[["CD"]][["ddf.2"]]$ds.param[n,1:2]), equals(as.numeric(ddf.results[[1]]$ds$aux$ddfobj$scale$parameters)))
  expect_that(bootstrap.ddf.statistics[["CD"]][["ddf.2"]]$AIC[n] < bootstrap.ddf.statistics[["CD"]][["ddf.1"]]$AIC[n], is_true())
  expect_that(ddf.results[[1]]$criterion, equals(bootstrap.ddf.statistics[["CD"]][["ddf.2"]]$AIC[n]))
  ##################################

  dht.results <- calculate.dht(species.code, ddf.model.options$species.field.name, model.index, ddf.results, region.table, sample.table, obs.table, dht.options)

  ##################################
  expect_that(names(dht.results), is_identical_to(c("CD","WD","UnidDol")))
  expect_that(dht.results[[1]]$clusters$summary$n[1]+dht.results[[2]]$clusters$summary$n[1]+dht.results[[3]]$clusters$summary$n[1], equals(nrow(obs.table)))
  ##################################

  if(unidentified.species){
    formatted.dht.results <- prorate.unidentified(dht.results, species.code.definitions, species.presence, clusters)
  }else{
    formatted.dht.results <- format.dht.results(dht.results, species.code, clusters)
  }

  ##################################
  expect_that(length(formatted.dht.results), equals(2))
  expect_that(names(formatted.dht.results), is_identical_to(c("CD","WD")))
  expect_that(dht.results[[1]]$clusters$N$Estimate[1]+dht.results[[2]]$clusters$N$Estimate[1]+dht.results[[3]]$clusters$N$Estimate[1], equals(formatted.dht.results[[1]]$clusters$N$Estimate[1]+formatted.dht.results[[2]]$clusters$N$Estimate[1]))
  expect_that(as.numeric(((formatted.dht.results[["CD"]]$clusters$N$Estimate[1]-dht.results[["CD"]]$clusters$N$Estimate[1])/formatted.dht.results[["CD"]]$clusters$N$Estimate[1])*100), equals(formatted.dht.results[["CD"]]$clusters$N$PercentUnidentified[1], tolerance = 0.0001))
  ##################################

  bootstrap.results <- accumulate.results(n, bootstrap.results, formatted.dht.results, clusters)

  ##################################
  expect_that(bootstrap.results$clusters.N["Total","PercentUnidentified",1,"CD"], equals(bootstrap.results$clusters.N["Total","PercentUnidentified",1,"WD"]))
  expect_that(bootstrap.results$clusters.N["Total","Estimate",1,"WD"], equals(as.numeric(formatted.dht.results[["WD"]]$clusters$N$Estimate[1])))
  expect_that(bootstrap.results$individual.N["Total","PercentUnidentified",1,"CD"], equals(as.numeric(((bootstrap.results$individual.N["Total","Estimate",1,"CD"]- dht.results[["CD"]]$individual$N$Estimate[1])/bootstrap.results$individual.N["Total","Estimate",1,"CD"])*100), tolerance = 0.001))
  ##################################

  n=2
  #Resample Data
  bootstrap = TRUE
  if(bootstrap){
    ddf.dat.working <- resample.data(resample=bootstrap.options$resample, obs.table.master, sample.table.master, ddf.dat.master, double.observer)
    obs.table       <- ddf.dat.working$obs.table
    sample.table    <- ddf.dat.working$sample.table
    ddf.dat.working <- ddf.dat.working$ddf.dat.working
  }else{
    ddf.dat.working <- ddf.dat.master
  }

  ##################################
  expect_that(length(unique(sample.table$Sample.Label)), equals(length(unique(sample.table.master$Sample.Label))))
  expect_that(table(sample.table$Region), is_identical_to(table(sample.table.master$Region)))
  expect_that(nrow(ddf.dat.working[[1]]), equals(nrow(obs.table)))
  expect_that(length(which(ddf.dat.working[[1]]$object%in%obs.table$object)), equals(nrow(obs.table)))
  expect_that(ddf.dat.working[["CD"]]$distance[ddf.dat.working[["CD"]]$object == 11], equals(ddf.dat.master[["CD"]]$distance[ddf.dat.master[["CD"]]$object == 11]))
  ##################################

  ddf.dat.working.check <- ddf.dat.working
  if(!is.null(covariate.uncertainty)){
    ddf.dat.working <- resample.covariates(ddf.dat.working, covariate.uncertainty, MAE.warnings)
    MAE.warnings <- ddf.dat.working$MAE.warnings
    ddf.dat.working <- ddf.dat.working$ddf.dat.working
  }

  ##################################
  expect_that(ddf.dat.working[["CD"]], is_identical_to(ddf.dat.working.check[["CD"]]))
  rm(ddf.dat.working.check)
  ##################################

  #Fit ddf models to all species codes
  ddf.results <- fit.ddf.models(ddf.dat.working, unique.model.names, ddf.models, ddf.model.options$criterion, bootstrap.ddf.statistics, n, MAE.warnings)
  if(class(ddf.results) == "list"){
    bootstrap.ddf.statistics <- ddf.results$bootstrap.ddf.statistics
    ddf.results <- ddf.results$ddf.results
  }else{
    #If the ddf results are not valid for all species move to next bootstrap iteration
    MAE.warnings <- ddf.results
    next
  }

  ##################################
  expect_that(as.numeric(bootstrap.ddf.statistics[["CD"]][["ddf.1"]]$ds.param[n,1:2]), equals(as.numeric(ddf.results[[1]]$ds$aux$ddfobj$scale$parameters)))
  expect_that(bootstrap.ddf.statistics[["CD"]][["ddf.2"]]$AIC[n] > bootstrap.ddf.statistics[["CD"]][["ddf.1"]]$AIC[n], is_true())
  expect_that(ddf.results[[1]]$criterion, equals(bootstrap.ddf.statistics[["CD"]][["ddf.1"]]$AIC[n]))
  ##################################

  dht.results <- calculate.dht(species.code, ddf.model.options$species.field.name, model.index, ddf.results, region.table, sample.table, obs.table, dht.options)
  if(unidentified.species){
    formatted.dht.results <- prorate.unidentified(dht.results, species.code.definitions, species.presence, clusters)
  }else{
    formatted.dht.results <- format.dht.results(dht.results, species.code, clusters)
  }
  if(unidentified.species){
    formatted.dht.results <- prorate.unidentified(dht.results, species.code.definitions, species.presence, clusters)
  }else{
    formatted.dht.results <- format.dht.results(dht.results, species.code, clusters)
  }

  ##################################
  expect_that(length(formatted.dht.results), equals(2))
  expect_that(names(formatted.dht.results), is_identical_to(c("CD","WD")))
  expect_that(dht.results[[1]]$clusters$N$Estimate[1]+dht.results[[2]]$clusters$N$Estimate[1]+dht.results[[3]]$clusters$N$Estimate[1], equals(formatted.dht.results[[1]]$clusters$N$Estimate[1]+formatted.dht.results[[2]]$clusters$N$Estimate[1]))
  expect_that(as.numeric(((formatted.dht.results[["CD"]]$clusters$N$Estimate[1]-dht.results[["CD"]]$clusters$N$Estimate[1])/formatted.dht.results[["CD"]]$clusters$N$Estimate[1])*100), equals(formatted.dht.results[["CD"]]$clusters$N$PercentUnidentified[1], tolerance = 0.0001))
  ##################################

  bootstrap.results <- accumulate.results(n, bootstrap.results, formatted.dht.results, clusters)

  ##################################
  expect_that(bootstrap.results$clusters.N["Total","PercentUnidentified",2,"CD"], equals(bootstrap.results$clusters.N["Total","PercentUnidentified",2,"WD"]))
  expect_that(bootstrap.results$individual.N["Total","PercentUnidentified",2,"WD"], equals(as.numeric(((bootstrap.results$individual.N["Total","Estimate",2,"WD"]- dht.results[["WD"]]$individual$N$Estimate[1])/bootstrap.results$individual.N["Total","Estimate",2,"WD"])*100), tolerance = 0.001))
  expect_that(bootstrap.results$Expected.S["Total","new.Expected.S",2,"CD"], equals(as.numeric(formatted.dht.results[["CD"]]$individual$N$Estimate[1]/formatted.dht.results[["CD"]]$clusters$N$Estimate[1])))
  ##################################

  #process results
  results <- process.bootstrap.results(bootstrap.results, model.index, clusters, bootstrap.ddf.statistics, bootstrap.options$quantile.type, analysis.options = list(bootstrap = bootstrap, n = bootstrap.options$n, covariate.uncertainty = covariate.uncertainty, clusters = clusters, double.observer = double.observer, unidentified.species = unidentified.species, species.code.definitions = species.code.definitions, model.names = model.names))

  class(results) <- "ma"
  class(results$analysis.options) <- "ma.analysis"
  class(results$species) <- "ma.allspecies"
  for(sp in seq(along = results$species)){
    class(results$species[[sp]]) <- "ma.species"
  }
  if(!is.null(results$unidentified)){
    class(results$unidentified) <- "ma.allunid"
    for(sp in seq(along = results$unidentified)){
      class(results$unidentified[[sp]]) <- "ma.unid"
    }
  }

  ##################################
  expect_that(results, is_identical_to(results.to.compare))
  ##################################

  #rm(.Random.seed)
})








