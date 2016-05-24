# Tests that variables can be reliably uploaded and then detected
# Author David Barker <support@causata.com>

library(testthat)
library(RCurl)
library(rjson)
library(XML)
library(stringr)
library(Causata)

source("utils.R")

context("ReliableModelCreation")

has.causata.connection <- has.causata.connection.function()

if (has.causata.connection){
  # we need a causata connection to run the tests below, test for a connection

  test_that("model upload is reliable", {
    randomName <- function() paste(sample(letters[1:26], 20, replace=TRUE), collapse="")
            
    createVariables <- function(variables) {
      for (name in names(variables)) {
        Config.CreatePrimaryVariable(causata.config, name, variable.expression=variables[[name]])
      }
    }
            
    detectVariables <- function(conn, variables) {
      for (name in names(variables)) {
        wait.for.variable(conn, name)
      }
    }
            
    deleteVariables <- function(variables) {
      for (name in names(variables)) {
        Config.DeleteVariable(causata.config, name)
      }
    }
    
    conn <- NULL
    tryCatch({
      
      conn <- local.connection()
        
      causata.variables <- list(
        "customer-id"="SELECT PROPERTY customer-id",
        "page-view-count"="COUNT page-view",
        "purchase-count"="COUNT purchase")
      
      createVariables(causata.variables)
      detectVariables(conn, causata.variables)
      
      query <- Query() + WithVariables('page-view-count', 'purchase-count') + Limit(10000)
      causata.data <- GetCausataData(conn, query)
      causata.data <- CleanNaFromContinuousP(causata.data)
      formula.obj <- formula("`purchase.count_All.Past` ~ `page.view.count_All.Past`")
      customers.matrix <- model.matrix(formula.obj, data=causata.data$df)    
      model <- cv.glmnet(x=customers.matrix, y=causata.data$df[["purchase.count_All.Past"]], family="gaussian", alpha=0.9)
      model.definition <- ModelDefinition(model, causata.data, formula.obj, model$lambda.min)
      
      
      for (i in 1:100) {
        vars <- list()
        tryCatch({
          variable.definition <- VariableDefinition(randomName())
          vars[[variable.definition$name]] <- "THERE IS NO EXPRESSION, IT'S A MODEL"
          
          upload.succeeded <- UploadModel(causata.config, model.definition, variable.definition)
          
          detectVariables(conn, vars)
        }, finally={
          deleteVariables(vars)
        })
      }
      
    }, finally={
      deleteVariables(causata.variables)
      close.connection(conn)
    })
  })

} # end of if for causata connection
