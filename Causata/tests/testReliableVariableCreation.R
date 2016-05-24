# Tests that variables can be reliably uploaded and then detected
# This test is only executed if the causata_environment.R file can be sourced
# Author David Barker <support@causata.com>

library(testthat)
library(RCurl)
library(rjson)
library(XML)
library(stringr)
library(Causata)

source("utils.R")
has.causata.connection <- has.causata.connection.function()

context("ReliableVariableCreation")

if (has.causata.connection){
  test_that("variable upload is reliable", {
    randomName <- function() paste(sample(letters[1:26], 20, replace=TRUE), collapse="")
            
    createVariables <- function(variables) {
      for (name in names(variables)) {
        Config.CreatePrimaryVariable(causata.config, variable.name=name, variable.expression=variables[[name]])
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
        
      for (i in 1:100) {
        variables <- list()
        variables[[randomName()]] <- "SELECT PROPERTY customer-id"
        
        tryCatch({
          createVariables(variables)
          detectVariables(conn, variables)
        }, finally={
          deleteVariables(variables)
        })
      }
      
    }, finally={
      close.connection(conn)
    })
  })
} # end of if has.causata.connection