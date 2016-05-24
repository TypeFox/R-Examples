library(mrds)
library(testthat)

context("Data Resample Checks")

test_that("Resamples the data correctly", {     

  #datasetup
  ex.filename<-system.file("testData/input_checks/ddf_dat.robj", package="mads")
  load(ex.filename)
  ex.filename<-system.file("testData/input_checks/obs_table.robj", package="mads")
  load(ex.filename)
  ex.filename<-system.file("testData/input_checks/region_table.robj", package="mads")
  load(ex.filename)
  ex.filename<-system.file("testData/input_checks/sample_table.robj", package="mads")
  load(ex.filename)
  
  #fit models to full dataset
  ddf.1 <- ddf(dsmodel = ~mcds(key = "hn", formula = ~ size), method='ds', data=ddf.dat,meta.data=list(width=4)) 
  ddf.1$data$detected <- rep(1, nrow(ddf.1$data))
  ddf.4 <- ddf(dsmodel = ~mcds(key = "hn", formula = ~ size), method='trial',mrmodel=~glm(link='logit',formula=~distance + size + sex + exposure), data=ddf.dat,meta.data=list(width=4)) 
  
  obs.table.master <- obs.table
  sample.table.master <- sample.table
  
  ddf.dat.master <- list("CD" = ddf.1$data) 

  #Resample data - resampling samples (ds analysis)                                                         
  ddf.dat.working <- resample.data(resample="samples", obs.table.master, sample.table.master, ddf.dat.master, double.observer = FALSE) 
  obs.table       <- ddf.dat.working$obs.table 
  sample.table    <- ddf.dat.working$sample.table
  ddf.dat.working <- ddf.dat.working$ddf.dat.working  

  ##################################
  expect_that(length(unique(sample.table$Sample.Label)), equals(length(unique(sample.table.master$Sample.Label))))
  expect_that(table(sample.table$Region), is_identical_to(table(sample.table.master$Region)))
  expect_that(nrow(ddf.dat.working[[1]]), equals(nrow(obs.table))) #note this is specific to the case where there is only one dataset
  expect_that(length(which(ddf.dat.working[[1]]$object%in%obs.table$object)), equals(nrow(obs.table)))
  
  #Expect that the distance associated with an observation has not changed
  #SOMETIMES FAILS DUE TO THE RENUMBERING. THE NEW ID NOT ALWAYS FOUND
  #object.id <- ddf.dat.working[["CD"]]$object[5]
  #expect_that(ddf.dat.working[["CD"]]$distance[ddf.dat.working[["CD"]]$object == object.id], equals(ddf.dat.master[["CD"]]$distance[ddf.dat.master[["CD"]]$object == object.id]))
  ##################################
  
  #Resample data - resampling individual observations (ds analysis)                                                           
  ddf.dat.working <- resample.data(resample="observations", obs.table.master, sample.table.master, ddf.dat.master, double.observer = FALSE) 
  obs.table       <- ddf.dat.working$obs.table 
  sample.table    <- ddf.dat.working$sample.table
  ddf.dat.working <- ddf.dat.working$ddf.dat.working  

  ##################################
  expect_that(length(unique(sample.table$Sample.Label)), equals(length(unique(sample.table.master$Sample.Label))))
  expect_that(table(sample.table$Region), is_identical_to(table(sample.table.master$Region)))
  expect_that(table(ddf.dat.working[[1]]$species), is_identical_to(table(ddf.dat.master[[1]]$species)))     
  expect_that(length(which(ddf.dat.working[[1]]$object%in%obs.table$object)), equals(nrow(ddf.dat.working[[1]])))
  
  #Expect that the distance associated with an observation has not changed
  #SOMETIMES FAILS DUE TO THE RENUMBERING. THE NEW ID NOT ALWAYS FOUND
  #object.id <- ddf.dat.working[["CD"]]$object[5]
  #expect_that(ddf.dat.working[["CD"]]$distance[ddf.dat.working[["CD"]]$object == object.id], equals(ddf.dat.master[["CD"]]$distance[ddf.dat.master[["CD"]]$object == object.id]))
  ##################################
  
  ddf.dat.master <- list("CD" = ddf.4$data) 
  
  #Resample data - resampling samples (mrds analysis)                                                            
  ddf.dat.working <- resample.data(resample="samples", obs.table.master, sample.table.master, ddf.dat.master, double.observer = TRUE) 
  obs.table       <- ddf.dat.working$obs.table 
  sample.table    <- ddf.dat.working$sample.table
  ddf.dat.working <- ddf.dat.working$ddf.dat.working
  
  #n.obs <- NULL
  #for(i in 1:999){
  #  ddf.dat.working <- resample.data(resample="samples", obs.table.master, sample.table.master, ddf.dat.master, double.observer = TRUE)$ddf.dat.working 
  #  n.obs <- c(n.obs, length(table(ddf.dat.working[[1]]$object)))
  #}
  #length(table(ddf.dat.master[[1]]$object))
  #[1] 162
  #mean(n.obs)
  #[1] 162.1632

  
  ##################################
  #Expect that there are still the same number of samplers
  expect_that(length(unique(sample.table$Sample.Label)), equals(length(unique(sample.table.master$Sample.Label))))
  #Expect that there are still the same number of samplers in each strata
  expect_that(table(sample.table$Region), is_identical_to(table(sample.table.master$Region)))
  #Expect that there are 2 rows for each object (double observer)
  expect_that(length(which(table(ddf.dat.master[[1]]$object) != 2)), equals(0))
  #Expect that there is an entry in the obs table for all objects in the ddf data
  expect_that(length(which(ddf.dat.working[[1]]$object%in%obs.table$object)), equals(nrow(ddf.dat.working[[1]])))
  #Expect that the distance associated with an observation has not changed
  #SOMETIMES FAILS DUE TO THE RENUMBERING. THE NEW ID NOT ALWAYS FOUND
  #object.id <- ddf.dat.working[["CD"]]$object[5]
  #expect_that(ddf.dat.working[["CD"]]$distance[ddf.dat.working[["CD"]]$object == object.id], equals(ddf.dat.master[["CD"]]$distance[ddf.dat.master[["CD"]]$object == object.id]))
  ##################################
  
  #Resample data - resampling individual observations (mrds analysis)                                                            
  ddf.dat.working <- resample.data(resample="observations", obs.table.master, sample.table.master, ddf.dat.master, double.observer = TRUE) 
  obs.table       <- ddf.dat.working$obs.table 
  sample.table    <- ddf.dat.working$sample.table
  ddf.dat.working <- ddf.dat.working$ddf.dat.working    
  
  ##################################
  #Expect that there are still the same number of samplers
  expect_that(length(unique(sample.table$Sample.Label)), equals(length(unique(sample.table.master$Sample.Label))))
  #Expect that there are still the same number of samplers in each strata
  expect_that(table(sample.table$Region), is_identical_to(table(sample.table.master$Region)))
  #Expect that there are the same number of observations
  expect_that(length(unique(ddf.dat.master[[1]]$object)), equals(length(unique(ddf.dat.working[[1]]$object))))
  #Expect that there are 2 rows for each object (double observer)
  expect_that(length(which(table(ddf.dat.master[[1]]$object) != 2)), equals(0))
  #Expect that there is an entry in the obs table for all objects in the ddf data
  expect_that(length(which(ddf.dat.working[[1]]$object%in%obs.table$object)), equals(nrow(ddf.dat.working[[1]])))
  
  #Expect that the distance associated with an observation has not changed
  #SOMETIMES FAILS DUE TO THE RENUMBERING. THE NEW ID NOT ALWAYS FOUND
  #object.id <- ddf.dat.working[["CD"]]$object[5]
  #expect_that(ddf.dat.working[["CD"]]$distance[ddf.dat.working[["CD"]]$object == object.id], equals(ddf.dat.master[["CD"]]$distance[ddf.dat.master[["CD"]]$object == object.id]))
  
  #Expect that there are the same number of sightings of each species
  expect_that(table(ddf.dat.working[[1]]$species), is_identical_to(table(ddf.dat.master[[1]]$species)))  
  #Expect that the labels have not changed
  expect_that(table(ddf.dat.working[[1]]$label), is_identical_to(table(ddf.dat.master[[1]]$label)))     
  ##################################

})  
  


