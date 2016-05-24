# Test uploading a glmnet model
# 
# This test is only run if the causata_environment.R file can be sourced.
# 
# Author: davidb
###############################################################################
library(testthat)
library(RCurl)
library(rjson)
library(glmnet)
library(XML)
library(stringr)
library(Causata)

source("utils.R")
has.causata.connection <- has.causata.connection.function()

context("UploadModel")

test1 <- function() {
  test_that("We can build a glmnet model, upload it, get scores from it, and that the match the scores from predict()", 
          
  # Connect to causata
  #
  with.local.connection(function(conn) {
    
    # Upload a known set of variables
    #
    causata.variables <- list(
      "customer-id"="SELECT PROPERTY customer-id",
      "page-view-count"="COUNT page-view",
      "purchase-count"="COUNT purchase")
    
    with.primary.variables(conn, causata.config, causata.variables, {
      
      # Get some data to use in training a model. (We will predict number of purchases from number of page views)
      #
      query <- Query() + WithVariables('page-view-count', 'purchase-count') + Limit(10000)
      causata.data <- GetCausataData(conn, query)
      causata.data <- CleanNaFromContinuousP(causata.data)
      formula.obj <- formula("`purchase.count_All.Past` ~ `page.view.count_All.Past`")
      customers.matrix <- model.matrix(formula.obj, data=causata.data$df)    
      model <- cv.glmnet(x=customers.matrix, y=causata.data$df[["purchase.count_All.Past"]], family="gaussian", alpha=0.9)
      
      variable.definition <- VariableDefinition("model")
      model.definition <- ModelDefinition(model, causata.data, formula.obj, model$lambda.min)
      
      upload.succeeded <- UploadModel(causata.config, model.definition, variable.definition)
      if (!upload.succeeded) {
        print("model upload failed")
        stop("model upload failed")
      }
      print("Created variable model")
      
      delete.variable.at.end(causata.config, "model", {
        # Now, check that the model matches the predictions from R
        #
        wait.for.variable(conn, "model")
        test.causata.data <- GetCausataData(conn, query + WithVariables("model") + Limit(50))
        test.df = test.causata.data$df
        
        # Make the same missing value substitutions, then convert into a matrix that predict.cv.glmnet can accept
        #
        purchase.missing.value <- causata.data$variableList[['purchase.count_All.Past']]$missingReplacement
        test.df[["purchase.count_All.Past"]][is.na(test.df[["purchase.count_All.Past"]])] <- purchase.missing.value
        test.matrix <- model.matrix(formula.obj, data=test.df)
        
        predictions <- predict(model, newx=test.matrix, s="lambda.min")
        errors <- test.df[["model_Time.Independent"]] - as.vector(predictions)
        large.errors <- errors[errors > 0.00000000000001]
        stopifnot(length(large.errors) == 0)
      })
    })
  })
)
}

# This is a copy-paste of the other test.  The only change is to use UploadModelWithValidation instead
# of UploadModel.  I've copied the code instead of making function that accepts an upload function
# because the output when a test fails is already bad for debugging purposes, and further indirection
# will not improve the situation
#
test2 <- function() {
  test_that("We can build a glmnet model, upload it, get scores from it, and that the match the scores from predict()", 
  
  # Connect to causata
  #
  with.local.connection(function(conn) {
    
    # Upload a known set of variables
    #
    causata.variables <- list(
      "customer-id"="SELECT PROPERTY customer-id",
      "page-view-count"="COUNT page-view",
      "purchase-count"="COUNT purchase")
    
    with.primary.variables(conn, causata.config, causata.variables, {
      
      # Get some data to use in training a model. (We will predict number of purchases from number of page views)
      #
      query <- Query() + WithVariables('page-view-count', 'purchase-count') + Limit(10000)
      causata.data <- GetCausataData(conn, query)
      causata.data <- CleanNaFromContinuousP(causata.data)
      formula.obj <- formula("`purchase.count_All.Past` ~ `page.view.count_All.Past`")
      customers.matrix <- model.matrix(formula.obj, data=causata.data$df)    
      model <- cv.glmnet(x=customers.matrix, y=causata.data$df[["purchase.count_All.Past"]], family="gaussian", alpha=0.9)
      
      variable.definition <- VariableDefinition("model")
      model.definition <- ModelDefinition(model, causata.data, formula.obj, model$lambda.min)
      
      upload.succeeded <- UploadModelWithValidation(
        causata.config,
        model.definition,
        variable.definition,
        conn,
        QueryFunction, # FIXME: add a query function
        0.0000000000001
      )
      
      if (!is.logical(upload.succeeded) || upload.succeeded != TRUE) {
        print("model upload failed")
        stop("model upload failed")
      }
      print("Created variable model")
      
      delete.variable.at.end(causata.config, "model", {
        # Now, check that the model matches the predictions from R
        #
        wait.for.variable(conn, "model")
        test.causata.data <- GetCausataData(conn, query + WithVariables("model") + Limit(50))
        test.df = test.causata.data$df
        
        # Make the same missing value substitutions, then convert into a matrix that predict.cv.glmnet can accept
        #
        purchase.missing.value <- causata.data$variableList[['purchase.count_All.Past']]$missingReplacement
        test.df[["purchase.count_All.Past"]][is.na(test.df[["purchase.count_All.Past"]])] <- purchase.missing.value
        test.matrix <- model.matrix(formula.obj, data=test.df)
        
        predictions <- predict(model, newx=test.matrix, s="lambda.min")
        errors <- test.df[["model_Time.Independent"]] - as.vector(predictions)
        large.errors <- errors[errors > 0.0000000000001]
        stopifnot(length(large.errors) == 0)
      })
    })
  })
)
}

if (has.causata.connection){
  for (i in 1:100) {
    test1()
    test2()
  }
}