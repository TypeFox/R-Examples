#
# file: tests/testthat/test_mcSimulation.R
#
# This file is part of the R-package decisionSupport
# 
# Authors: 
#   Lutz Göhring <lutz.goehring@gmx.de>
#   Eike Luedeling (ICRAF) <eike@eikeluedeling.com>
#
# Copyright (C) 2015 World Agroforestry Centre (ICRAF) 
#	http://www.worldagroforestry.org
# 
# The R-package decisionSupport is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# The R-package decisionSupport is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with the R-package decisionSupport.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################################
context("Testing mcSimulation()")

set.seed(100)
# Number of model runs to perform the Monte Carlo simulation:
n= 1000000
tolerance=3/sqrt(n)

test_that("Difference of two uncorrelated normally distributed variables is normally distributed
          (randomMethod=\"calculate\", functionSyntax=\"data.frameNames\") (1).",{
            # Define the model for the profit:
            profitModel <- function(x){
              x$revenue-x$costs
            }
            # Read the estimate for revenue and costs:
            profitEstimate<-estimate_read_csv("profit-1.csv")
            # Calculate means from 95%-confidence intervalls:
            meanRevenue <- mean(c(profitEstimate$marginal["revenue","lower"], profitEstimate$marginal["revenue","upper"]) ) 
            meanCosts <- mean(c(profitEstimate$marginal["costs","lower"], profitEstimate$marginal["costs","upper"]) ) 
            # Calculate standard deviations from 95%-confidence intervalls:
            sdRevenue <- 0.5 * (profitEstimate$marginal["revenue","upper"] - profitEstimate$marginal["revenue","lower"]) / qnorm(0.95)
            sdCosts <- 0.5 * (profitEstimate$marginal["costs","upper"] - profitEstimate$marginal["costs","lower"]) / qnorm(0.95)
            # Calculate expected moments for profit = revenue - costs:
            meanProfitExpected <- meanRevenue - meanCosts
            sdProfitExpected <- sqrt( sdRevenue^2 + sdCosts^2)
            # Run the Monte Carlo Simulation:
            profitSimulation <- mcSimulation(estimate=profitEstimate, 
                                             model_function=profitModel, 
                                             numberOfModelRuns=n,
                                             randomMethod="calculate",
                                             functionSyntax="data.frameNames")
            expect_equal(mean(profitSimulation$y$output_1), meanProfitExpected, tolerance=tolerance)
            expect_equal(sd(profitSimulation$y$output_1), sdProfitExpected, tolerance=tolerance)
          })
test_that("Difference of two uncorrelated normally distributed variables is normally distributed
          (randomMethod=\"calculate\", functionSyntax=\"data.frameNames\") (2).",{
            # Define the model for the profit:
            profitModel <- function(x){
              x[["revenue"]]-x[["costs"]]
            }
            # Read the estimate for revenue and costs:
            profitEstimate<-estimate_read_csv("profit-1.csv")
            # Calculate means from 95%-confidence intervalls:
            meanRevenue <- mean(c(profitEstimate$marginal["revenue","lower"], profitEstimate$marginal["revenue","upper"]) ) 
            meanCosts <- mean(c(profitEstimate$marginal["costs","lower"], profitEstimate$marginal["costs","upper"]) ) 
            # Calculate standard deviations from 95%-confidence intervalls:
            sdRevenue <- 0.5 * (profitEstimate$marginal["revenue","upper"] - profitEstimate$marginal["revenue","lower"]) / qnorm(0.95)
            sdCosts <- 0.5 * (profitEstimate$marginal["costs","upper"] - profitEstimate$marginal["costs","lower"]) / qnorm(0.95)
            # Calculate expected moments for profit = revenue - costs:
            meanProfitExpected <- meanRevenue - meanCosts
            sdProfitExpected <- sqrt( sdRevenue^2 + sdCosts^2)
            # Run the Monte Carlo Simulation:
            profitSimulation <- mcSimulation(estimate=profitEstimate, 
                                             model_function=profitModel, 
                                             numberOfModelRuns=n,
                                             randomMethod="calculate",
                                             functionSyntax="data.frameNames")
            expect_equal(mean(profitSimulation$y$output_1), meanProfitExpected, tolerance=tolerance)
            expect_equal(sd(profitSimulation$y$output_1), sdProfitExpected, tolerance=tolerance)
          })
test_that("Difference of two uncorrelated normally distributed variables is normally distributed
          (randomMethod=\"calculate\", functionSyntax=\"matrixNames\").",{
            # Define the model for the profit:
            profitModel <- function(x){
              x[,"revenue"]-x[,"costs"]
            }
            # Read the estimate for revenue and costs:
            profitEstimate<-estimate_read_csv("profit-1.csv")
            # Calculate means from 95%-confidence intervalls:
            meanRevenue <- mean(c(profitEstimate$marginal["revenue","lower"], profitEstimate$marginal["revenue","upper"]) ) 
            meanCosts <- mean(c(profitEstimate$marginal["costs","lower"], profitEstimate$marginal["costs","upper"]) ) 
            # Calculate standard deviations from 95%-confidence intervalls:
            sdRevenue <- 0.5 * (profitEstimate$marginal["revenue","upper"] - profitEstimate$marginal["revenue","lower"]) / qnorm(0.95)
            sdCosts <- 0.5 * (profitEstimate$marginal["costs","upper"] - profitEstimate$marginal["costs","lower"]) / qnorm(0.95)
            # Calculate expected moments for profit = revenue - costs:
            meanProfitExpected <- meanRevenue - meanCosts
            sdProfitExpected <- sqrt( sdRevenue^2 + sdCosts^2)
            # Run the Monte Carlo Simulation:
            profitSimulation <- mcSimulation(estimate=profitEstimate, 
                                             model_function=profitModel, 
                                             numberOfModelRuns=n,
                                             randomMethod="calculate",
                                             functionSyntax="matrixNames")
            expect_equal(mean(profitSimulation$y$output_1), meanProfitExpected, tolerance=tolerance)
            expect_equal(sd(profitSimulation$y$output_1), sdProfitExpected, tolerance=tolerance)
          })
test_that("Difference of two uncorrelated normally distributed variables is normally distributed
          (randomMethod=\"calculate\", functionSyntax=\"plainNames\") (1).",{
            # Reset number of model runs to perform the Monte Carlo simulation 
            #  because "plainNames" is slower than "data.frameNames" and "matrixNames":
            n = 10000
            tolerance=3/sqrt(n)
            # Define the model for the profit:
            profitModel <- function() revenue - costs
            # Read the estimate for revenue and costs:
            profitEstimate<-estimate_read_csv("profit-1.csv")
            # Calculate means from 95%-confidence intervalls:
            meanRevenue <- mean(c(profitEstimate$marginal["revenue","lower"], profitEstimate$marginal["revenue","upper"]) ) 
            meanCosts <- mean(c(profitEstimate$marginal["costs","lower"], profitEstimate$marginal["costs","upper"]) ) 
            # Calculate standard deviations from 95%-confidence intervalls:
            sdRevenue <- 0.5 * (profitEstimate$marginal["revenue","upper"] - profitEstimate$marginal["revenue","lower"]) / qnorm(0.95)
            sdCosts <- 0.5 * (profitEstimate$marginal["costs","upper"] - profitEstimate$marginal["costs","lower"]) / qnorm(0.95)
            # Calculate expected moments for profit = revenue - costs:
            meanProfitExpected <- meanRevenue - meanCosts
            sdProfitExpected <- sqrt( sdRevenue^2 + sdCosts^2)
            # Run the Monte Carlo Simulation:
            profitSimulation <- mcSimulation(estimate=profitEstimate, 
                                             model_function=profitModel, 
                                             numberOfModelRuns=n,
                                             randomMethod="calculate",
                                             functionSyntax="plainNames")
            expect_equal(mean(profitSimulation$y$output_1), meanProfitExpected, tolerance=tolerance)
            expect_equal(sd(profitSimulation$y$output_1), sdProfitExpected, tolerance=tolerance)
          })
test_that("Difference of two correlated normally distributed variables is normally distributed
          (randomMethod=\"calculate\", functionSyntax=\"data.frameNames\") (1).",{
            # Define the model for the profit:
            profitModel <- function(x){
              x$revenue-x$costs
            }
            # Read the estimate for revenue and costs:
            profitEstimate<-estimate_read_csv("profit-2.csv")
            # Calculate means from 95%-confidence intervalls:
            meanRevenue <- mean(c(profitEstimate$marginal["revenue","lower"], profitEstimate$marginal["revenue","upper"]) ) 
            meanCosts <- mean(c(profitEstimate$marginal["costs","lower"], profitEstimate$marginal["costs","upper"]) ) 
            # Calculate standard deviations from 95%-confidence intervalls:
            sdRevenue <- 0.5 * (profitEstimate$marginal["revenue","upper"] - profitEstimate$marginal["revenue","lower"]) / qnorm(0.95)
            sdCosts <- 0.5 * (profitEstimate$marginal["costs","upper"] - profitEstimate$marginal["costs","lower"]) / qnorm(0.95)
            covRevenueCosts <- sdRevenue * sdCosts * profitEstimate$correlation_matrix["revenue","costs"]
            # Calculate expected moments for profit = revenue - costs:
            meanProfitExpected <- meanRevenue - meanCosts
            sdProfitExpected <- sqrt( sdRevenue^2 - 2*covRevenueCosts + sdCosts^2)
            # Run the Monte Carlo Simulation:
            profitSimulation <- mcSimulation(estimate=profitEstimate, 
                                             model_function=profitModel, 
                                             numberOfModelRuns=n,
                                             randomMethod="calculate",
                                             functionSyntax="data.frameNames")
            expect_equal(mean(profitSimulation$y$output_1), meanProfitExpected, tolerance=tolerance)
            expect_equal(sd(profitSimulation$y$output_1), sdProfitExpected, tolerance=tolerance)
          })
test_that("Difference of two correlated normally distributed variables is normally distributed
          (randomMethod=\"calculate\", functionSyntax=\"plainNames\") (1).",{
            # Reset number of model runs to perform the Monte Carlo simulation 
            #  because "plainNames" is slower than "data.frameNames" and "matrixNames":
            n = 10000
            tolerance=3/sqrt(n)
            # Define the model for the profit:
            profitModel <- function() revenue - costs
            # Read the estimate for revenue and costs:
            profitEstimate<-estimate_read_csv("profit-2.csv")
            # Calculate means from 95%-confidence intervalls:
            meanRevenue <- mean(c(profitEstimate$marginal["revenue","lower"], profitEstimate$marginal["revenue","upper"]) ) 
            meanCosts <- mean(c(profitEstimate$marginal["costs","lower"], profitEstimate$marginal["costs","upper"]) ) 
            # Calculate standard deviations from 95%-confidence intervalls:
            sdRevenue <- 0.5 * (profitEstimate$marginal["revenue","upper"] - profitEstimate$marginal["revenue","lower"]) / qnorm(0.95)
            sdCosts <- 0.5 * (profitEstimate$marginal["costs","upper"] - profitEstimate$marginal["costs","lower"]) / qnorm(0.95)
            covRevenueCosts <- sdRevenue * sdCosts * profitEstimate$correlation_matrix["revenue","costs"]
            # Calculate expected moments for profit = revenue - costs:
            meanProfitExpected <- meanRevenue - meanCosts
            sdProfitExpected <- sqrt( sdRevenue^2 - 2*covRevenueCosts + sdCosts^2)
            # Run the Monte Carlo Simulation:
            profitSimulation <- mcSimulation(estimate=profitEstimate, 
                                             model_function=profitModel, 
                                             numberOfModelRuns=n,
                                             randomMethod="calculate",
                                             functionSyntax="plainNames")
            expect_equal(mean(profitSimulation$y$output_1), meanProfitExpected, tolerance=tolerance)
            expect_equal(sd(profitSimulation$y$output_1), sdProfitExpected, tolerance=tolerance)
          })
test_that("5 dimensional estimate and 2 dimensional named model function are simulated:
          (randomMethod=\"calculate\", functionSyntax=\"plainNames\") (1).",{
            # Number of simulations (only a small number needed for the testing that its running):
            n=10
            # Define some model for the profit and schnecke:
            estimate5dModel2d<-function(x){
              for(i in names(x)) assign(i, as.numeric(x[i]))
              
              a<-(b+c)^d+e+f-g
              aa<-b+c+d-e+f*g
              return(list(profit=a,schnecke=aa))
            }
            # Read the estimate for a, b, c, d, e and f from file:
            estimate5d<-estimate_read_csv("estimate5d.csv")
            
            # Run the Monte Carlo Simulation:
            profitSimulation <- mcSimulation(estimate=estimate5d, 
                                             model_function=estimate5dModel2d, 
                                             numberOfModelRuns=n,
                                             randomMethod="calculate",
                                             functionSyntax="plainNames")
          })
test_that("4 dimensional estimate and 2 dimensional unnamed model function are simulated:
          (randomMethod=\"calculate\", functionSyntax=\"plainNames\") (1).",{
            # Number of simulations:
            n=10
            # Create the current estimate from text:
            estimateText<-"variable,  distribution, lower, upper
               revenue1,  posnorm,      100,   1000
               revenue2,  posnorm,      50,    2000
               costs1,    posnorm,      50,    2000
               costs2,    posnorm,      100,   1000"
            estimate4d<-as.estimate(read.csv(header=TRUE, text=estimateText,
                                             strip.white=TRUE, stringsAsFactors=FALSE))
            # The welfare function:
            estimate4dModel2d <- function(x){
              # Assign the variable names to the function environement:
              for(i in names(x)) assign(i, as.numeric(x[i]))
              
              list(revenue1 + revenue2 - costs1 - costs2, revenue1)
            }
            # Calculate the Individual EVPI:
            
            # Run the Monte Carlo Simulation:
            profitSimulation <- mcSimulation(estimate=estimate4d, 
                                             model_function=estimate4dModel2d, 
                                             numberOfModelRuns=n,
                                             randomMethod="calculate",
                                             functionSyntax="plainNames",
                                             verbosity=0)
            expect_false(is.null(profitSimulation$y$output_1))
            expect_false(is.null(profitSimulation$y$output_2))
            expect_true( is.null(profitSimulation$y$output_3))
          })
test_that("2 dimensional estimate and 1 dimensional  model function returning an unnamed list are simulated:
          (randomMethod=\"calculate\", functionSyntax=\"plainNames\") (1).",{
            #########################################################
            # Create the estimate object:
            variable=c("revenue","costs")
            distribution=c("norm","norm")
            lower=c(10000,  5000)
            upper=c(100000, 50000)
            costBenefitEstimate<-as.estimate(variable, distribution, lower, upper)
            # Number of simulations:
            n=10
            profit1<-function(x){
              # Assign the variable names to the function environment:
              for(i in names(x)) assign(i, as.numeric(x[i]))
              #list(revenue-costs,revenue+costs)
              list(revenue-costs)
            }
            # Perform the Monte Carlo simulation:
            predictionProfit1<-mcSimulation( estimate=costBenefitEstimate,
                                             model_function=profit1,
                                             numberOfModelRuns=n,
                                             functionSyntax="plainNames")
            expect_false(is.null(predictionProfit1$y$output_1))
            expect_true(is.null(predictionProfit1$y$output_2))
            expect_true(is.null(predictionProfit1$y$output_3))
          })
test_that("2 dimensional estimate and 2 dimensional  model function returning an unnamed list are simulated:
          (randomMethod=\"calculate\", functionSyntax=\"data.frameNames\") (1).",{
            #########################################################
            # Create the estimate object:
            variable=c("revenue","costs")
            distribution=c("norm","norm")
            lower=c(10000,  5000)
            upper=c(100000, 50000)
            costBenefitEstimate<-as.estimate(variable, distribution, lower, upper)
            # Number of simulations:
            n=10
            # (a) Define the model function without name for the return value:
            profit1<-function(x){
              list(x[["revenue"]]-x[["costs"]],x[["revenue"]]+x[["costs"]])
            }
            # Perform the Monte Carlo simulation:
            profitSimulation<-mcSimulation( estimate=costBenefitEstimate,
                                             model_function=profit1,
                                             numberOfModelRuns=n,
                                             functionSyntax="data.frameNames",
                                             verbosity=0)
            expect_false(is.null(profitSimulation$y$output_1))
            expect_false(is.null(profitSimulation$y$output_2))
            expect_true( is.null(profitSimulation$y$output_3))
          })
test_that("5 dimensional estimate and 2 dimensional named model function are simulated:
          (randomMethod=\"calculate\", functionSyntax=\"plainNames\") (1).",{
            # Number of simulations (only a small number needed for the testing that its running):
            n=10
            # Define some model for the profit and schnecke:
            estimate5dModel2d<-function(){
              a<-(b+c)^d+e+f-g
              aa<-b+c+d-e+f*g
              return(list(profit=a,schnecke=aa))
            }
            # Read the estimate for a, b, c, d, e and f from file:
            estimate5d<-estimate_read_csv("estimate5d.csv")
            
            # Run the Monte Carlo Simulation:
            profitSimulation <- mcSimulation(estimate=estimate5d, 
                                             model_function=estimate5dModel2d, 
                                             numberOfModelRuns=n,
                                             randomMethod="calculate",
                                             functionSyntax="plainNames")
          })
test_that("4 dimensional estimate and 2 dimensional unnamed model function are simulated:
          (randomMethod=\"calculate\", functionSyntax=\"plainNames\") (1).",{
            # Number of simulations:
            n=10
            # Create the current estimate from text:
            estimateText<-"variable,  distribution, lower, upper
            revenue1,  posnorm,      100,   1000
            revenue2,  posnorm,      50,    2000
            costs1,    posnorm,      50,    2000
            costs2,    posnorm,      100,   1000"
            estimate4d<-as.estimate(read.csv(header=TRUE, text=estimateText,
                                             strip.white=TRUE, stringsAsFactors=FALSE))
            # The welfare function:
            estimate4dModel2d <- function(){
              list(revenue1 + revenue2 - costs1 - costs2, revenue1)
            }
            # Calculate the Individual EVPI:
            
            # Run the Monte Carlo Simulation:
            profitSimulation <- mcSimulation(estimate=estimate4d, 
                                             model_function=estimate4dModel2d, 
                                             numberOfModelRuns=n,
                                             randomMethod="calculate",
                                             functionSyntax="plainNames",
                                             verbosity=0)
            expect_false(is.null(profitSimulation$y$output_1))
            expect_false(is.null(profitSimulation$y$output_2))
            expect_true( is.null(profitSimulation$y$output_3))
          })
test_that("2 dimensional estimate and 1 dimensional  model function returning an unnamed list are simulated:
          (randomMethod=\"calculate\", functionSyntax=\"plainNames\") (1).",{
            #########################################################
            # Create the estimate object:
            variable=c("revenue","costs")
            distribution=c("norm","norm")
            lower=c(10000,  5000)
            upper=c(100000, 50000)
            costBenefitEstimate<-as.estimate(variable, distribution, lower, upper)
            # Number of simulations:
            n=10
            profit1<-function(){
              list(revenue-costs)
            }
            # Perform the Monte Carlo simulation:
            predictionProfit1<-mcSimulation( estimate=costBenefitEstimate,
                                             model_function=profit1,
                                             numberOfModelRuns=n,
                                             functionSyntax="plainNames")
            expect_false(is.null(predictionProfit1$y$output_1))
            expect_true(is.null(predictionProfit1$y$output_2))
            expect_true(is.null(predictionProfit1$y$output_3))
          })
