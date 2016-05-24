library(mrds)
library(testthat)

context("Data Input")

test_that("Test data input checks", {


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
  #ddf.3 <- ddf(dsmodel = ~mcds(key = "unif", formula = ~ 1, adj.series = "cos", adj.order = c(2)), method='ds', data=ddf.dat,meta.data=list(width=4))
  ddf.3 <- ddf(dsmodel = ~mcds(key = "hn", formula = ~ 1, adj.series = "cos", adj.order = c(2)), method='ds', data=ddf.dat,meta.data=list(width=4, mono=TRUE))
  ddf.4 <- ddf(dsmodel = ~mcds(key = "hn", formula = ~ size), method='trial',mrmodel=~glm(link='logit',formula=~distance + size + sex + exposure), data=ddf.dat,meta.data=list(width=4))
  ddf.5 <- ddf(dsmodel = ~mcds(key = "hn", formula = ~ size), method='trial.fi',mrmodel=~glm(link='logit',formula=~distance + size + sex + exposure), data=ddf.dat,meta.data=list(width=4))
  #ddf.6 <- ddf(dsmodel = ~mcds(key = "hn", formula = ~ size), method='io',mrmodel=~glm(link='logit',formula=~distance + size + sex + exposure), data=ddf.dat,meta.data=list(width=4))
  ddf.7 <- ddf(dsmodel = ~mcds(key = "hn", formula = ~ size), method='io.fi',mrmodel=~glm(link='logit',formula=~distance + size + sex + exposure), data=ddf.dat,meta.data=list(width=4))
  ddf.8 <- ddf.7
  class(ddf.8) <- c("rem", "ddf")
  ddf.9 <- ddf.7
  class(ddf.9) <- c("glm")
  ddf.10 <- ddf.3
  ddf.10$data  <- ddf.10$data[,-3]

  #set up information for mads
  #Multi-analysis options
  model.names              <- list("CD"=c("ddf.1","ddf.2"), "WD"=c("ddf.1","ddf.2"), "UnidDol"=c("ddf.1","ddf.2"))
  ddf.models               <- list("ddf.1" = ddf.1, "ddf.2" = ddf.2)
  species.code.definitions <- list("UnidDol" = c("CD","WD"))
  species.presence         <- list("1. A" = c("CD","WD"))
  covariate.uncertainty    <- NULL
  ddf.model.options        <- list(criterion="AIC")
  ddf.model.options$distance.naming.conv <- TRUE
  bootstrap                <- TRUE
  bootstrap.options        <- list(resample="samples", n=1)
  seed.array               <- NULL

  #~~~~~~~~~~~~~~~~~~~~~~~~ TEST check.ddf.models(...) ~~~~~~~~~~~~~~~~~~~~~~~~~

  # check that there is an error when the same model is selected multiple times for a given species
  model.names              <- list("CD"=c("ddf.1","ddf.1"), "WD"=c("ddf.1","ddf.2"), "UnidDol"=c("ddf.1","ddf.2"))
  ddf.models               <- list("ddf.1" = ddf.1, "ddf.2" = ddf.2)
  expect_that(results <- execute.multi.analysis(
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
                              silent = FALSE),
           throws_error("The model names are not unique for species CD."))

  # check that there is an error when ds models are mixed with io/trial
  model.names              <- list("CD"=c("ddf.1","ddf.4"), "WD"=c("ddf.1","ddf.2"), "UnidDol"=c("ddf.1","ddf.2"))
  ddf.models               <- list("ddf.1" = ddf.1, "ddf.2" = ddf.2, "ddf.4" = ddf.4)
  expect_that(results <- execute.multi.analysis(
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
                              silent = FALSE),
              throws_error("Unsupported model types have been selected: trial"))
  # check that there is an error when io models are mixed with trial
  model.names              <- list("CD"=c("ddf.4","ddf.5"), "WD"=c("ddf.4","ddf.5"), "UnidDol"=c("ddf.4","ddf.7"))
  ddf.models               <- list("ddf.4" = ddf.4, "ddf.5" = ddf.5, "ddf.7" = ddf.7)
  expect_that(results <- execute.multi.analysis(
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
                              silent = FALSE),
           throws_error("Unsupported model types have been selected: trial, trial.fi, trial, trial.fi, trial, io.fi"))
  # check that there is an error when models don't exists
  model.names              <- list("CD"=c("ddf.40","ddf.5"), "WD"=c("ddf.4","ddf.5"), "UnidDol"=c("ddf.4","ddf.5"))
  ddf.models               <- list("ddf.4" = ddf.4, "ddf.5" = ddf.5)
  expect_that(results <- execute.multi.analysis(
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
                              silent = FALSE),
           throws_error("ddf object 1, analysis name ddf.40, for species code CD has not been provided."))
  # check that there is an error when some models contain data with cluster size and some contain data withough cluster size
  model.names              <- list("CD"=c("ddf.2","ddf.3"), "WD"=c("ddf.2","ddf.3"), "UnidDol"=c("ddf.2","ddf.10"))
  ddf.models               <- list("ddf.2" = ddf.2, "ddf.3" = ddf.3, "ddf.10" = ddf.10)
  expect_that(results <- execute.multi.analysis(
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
                              silent = FALSE),
           throws_error("Cluster size must be present in all datasets within the ddf models or none."))

  #~~~~~~~~~~~~~~~~ TEST check.species.code.definitions(...) ~~~~~~~~~~~~~~~~~~~

  model.names              <- list("CD"=c("ddf.1","ddf.2"), "WD"=c("ddf.1","ddf.2"), "UnidDol"=c("ddf.1","ddf.2"))
  ddf.models               <- list("ddf.1" = ddf.1, "ddf.2" = ddf.2)
  species.code.definitions <- list("UnidDol"=c("CD","WD"), "CD"=NULL)
  expect_that(results <- execute.multi.analysis(
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
                              silent = FALSE),
           throws_error("No species codes specified for CD in the species code definitions list."))
  species.code.definitions <- list("UnidDol"=c("CD","WD"), "CP"="CP")
  expect_that(results <- execute.multi.analysis(
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
                              silent = FALSE),
           throws_error("Species mismatch in ddf models and species code definitions. Models not suppled for all species or models supplied for species not included in species code definitions."))
  species.code.definitions <- list("UnidDol"=c("CD","WD"), "CD"="WD")
  expect_that(results <- execute.multi.analysis(
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
                              silent = FALSE),
           throws_error("Incorrect species code definition for species CD. If only a single code is entered it must match the name of the list element."))
  species.code.definitions <- list("UnidDol"=c("UnidDol","CD","WD"))
  expect_that(results <- execute.multi.analysis(
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
                              silent = FALSE),
           throws_error("Incorrect species code definition for species UnidDol. Unidentified code cannot be prorated to itself."))
  species.code.definitions <- list("UnidDol"=c("CD","WD"), "UnidDel" = c("CD", "WD", "UnidDol"))
  model.names              <- list("CD"=c("ddf.1","ddf.2"), "WD"=c("ddf.1","ddf.2"), "UnidDol"=c("ddf.1","ddf.2"), "UnidDel"=c("ddf.1","ddf.2"))
  ddf.models               <- list("ddf.1" = ddf.1, "ddf.2" = ddf.2)
  expect_that(results <- execute.multi.analysis(
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
                              silent = FALSE),
           throws_error("Incorrect species code definition for species UnidDel. An Unidentified code cannot be prorated to another unidentified code."))
  species.code.definitions <- list("UnidDol"=c("CD","WD"), "UnidDol"=c("CD","HP"))
  model.names              <- list("CD"=c("ddf.1","ddf.2"), "WD"=c("ddf.1","ddf.2"), "UnidDol"=c("ddf.1","ddf.2"))
  ddf.models               <- list("ddf.1" = ddf.1, "ddf.2" = ddf.2)
  expect_that(results <- execute.multi.analysis(
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
                              silent = FALSE),
           throws_error("Multiple species code entries in the species code definitions list."))
  species.code.definitions <- list("UnidDol"=c("CD","WD"))
  model.names              <- list("CD"=c("ddf.1","ddf.2"), "WD"=c("ddf.1","ddf.2"), "UnidDol"=c("ddf.1","ddf.2"), "HP"=c("ddf.1","ddf.2"))
  ddf.models               <- list("ddf.1" = ddf.1, "ddf.2" = ddf.2)
  expect_that(results <- execute.multi.analysis(
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
                              silent = FALSE),
           throws_error("Species mismatch in ddf models and species code definitions. Models not suppled for all species or models supplied for species not included in species code definitions."))

  #~~~~~~~~~~~~~~~~~~~ TEST check.covar.uncertainty(...) ~~~~~~~~~~~~~~~~~~~~~~~

  model.names              <- list("CD"=c("ddf.1","ddf.2"), "WD"=c("ddf.1","ddf.2"), "UnidDol"=c("ddf.1","ddf.2"))
  ddf.models               <- list("ddf.1" = ddf.1, "ddf.2" = ddf.2)
  covariate.uncertainty = data.frame(variable.layer = c("observation"), variable.name = c("scaledtotsize"), cor.factor.layer = "numeric", cor.factor.name = 1, uncertainty.layer = c("observation"), uncertainty.name = c("totsizecv"), uncertainty.measure = c("CV"), sampling.distribution = c("Norm"))
  expect_that(results <- execute.multi.analysis(
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
                              silent = FALSE),
           throws_error("An unsupported sampling distribution has been chosen for covariate uncertainty. Only one of the following may be specified: Normal, Normal.Absolute, Lognormal.BC, Poisson, TruncPoisson.BC"))

  covariate.uncertainty = data.frame(variable.layer = c("observation"), variable.name = c("scaledtotsize"), cor.factor.layer = "numeric", cor.factor.name = 1, uncertainty.layer = c("observation"), uncertainty.name = c("totsizecv"), uncertainty.measure = c("CV"), sampling.distribution = c("Normal"))
  expect_that(results <- execute.multi.analysis(
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
                              silent = FALSE),
           throws_error("Invalid names for the covariates or associated uncertainty have been specified in the covariate uncertainty dataframe."))
  
  
  #Test that the species.field.name argument works
  ddf.1 <- ddf(dsmodel = ~mcds(key = "hn", formula = ~ size), method='ds', data=ddf.dat,meta.data=list(width=4))
  ddf.2 <- ddf(dsmodel = ~mcds(key = "hr", formula = ~ size), method='ds', data=ddf.dat,meta.data=list(width=4))
  model.names              <- list("CD"=c("ddf.1","ddf.2"), "WD"=c("ddf.1","ddf.2"), "UnidDol"=c("ddf.1","ddf.2"))
  ddf.models               <- list("ddf.1" = ddf.1, "ddf.2" = ddf.2)
  species.code.definitions <- list("UnidDol" = c("CD","WD"))
  species.presence         <- list("1. A" = c("CD","WD"))
  covariate.uncertainty    <- NULL
  ddf.model.options        <- list(criterion="AIC")
  ddf.model.options$distance.naming.conv <- TRUE
  bootstrap                <- FALSE
  bootstrap.options        <- list(resample="samples", n=1)
  seed.array               <- NULL
  
  set.seed(444)
  species.field.name <- execute.multi.analysis(
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
  
  #add a new column of a different name the same as species
  ddf.dat$test <- ddf.dat$species
  #Remove the species column
  new.dat <- ddf.dat[,-13]
  
  #re-run ddf analyses so they have the updated data
  ddf.1 <- ddf(dsmodel = ~mcds(key = "hn", formula = ~ size), method='ds', data=new.dat,meta.data=list(width=4))
  ddf.2 <- ddf(dsmodel = ~mcds(key = "hr", formula = ~ size), method='ds', data=new.dat,meta.data=list(width=4))
  #Use new models!
  ddf.models               <- list("ddf.1" = ddf.1, "ddf.2" = ddf.2)
  
  #provide the new species.field.name
  ddf.model.options <- list(criterion="AIC", species.field.name = "test")
  set.seed(444)
  alternative.field.name <- execute.multi.analysis(
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
  
  expect_identical(species.field.name, alternative.field.name)
  
})
