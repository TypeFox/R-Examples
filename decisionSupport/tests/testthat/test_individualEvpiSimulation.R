#
# file: tests/testthat/test_individualEvpiSimulation.R
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
context("Testing individualEvpiSimulation()")

set.seed(100)
# Number of model runs for the underlying Monte Carlo simulation:
n= 1000
tolerance=3/sqrt(n)

test_that("Individual EVPI is calculated for 2-D correltated current estimate.",{
  # Define the model for the profit:
  profitModel <- function(x){
    x$revenue-x$costs
  }
  # Read the estimate for revenue and costs:
  profitEstimate<-estimate_read_csv("profit-2.csv")
  # Run the EVPI Simulation:
  evpiSimulation <- individualEvpiSimulation(welfare=profitModel, 
                                             currentEstimate=profitEstimate, 
                                             numberOfModelRuns=n,
                                             randomMethod="calculate",
                                             functionSyntax="data.frameNames")
})
test_that("Individual EVPI is calculated for 3-D correlated current estimate.",{
  # Define the model for the profit:
  profitModel <- function(x){
    x$revenue1 + x$revenue2 -x$costs1
  }
  # Create an estimate from text (with correlated components):
  estimateTextMarg<-"variable,  distribution, lower,  upper
                     revenue1,          norm,   500,   1000
                     revenue2,          norm,    50,   2000
                       costs1,          norm,    50,   1100"
  estimateTextCor<-",         revenue1, revenue2, costs1
                  revenue1,          1,     -0.3,    0.5
                  revenue2,       -0.3,        1,    0.2
                  costs1,          0.5,      0.2,      1"
  profitEstimate<-as.estimate(read.csv(header=TRUE, text=estimateTextMarg,
                                       strip.white=TRUE, stringsAsFactors=FALSE),
                              correlation_matrix=data.matrix(read.csv(text=estimateTextCor,
                                                                      row.names=1,
                                                                      strip.white=TRUE)))
  # Run the EVPI Simulation:
  evpiSimulation <- individualEvpiSimulation(welfare=profitModel, 
                                             currentEstimate=profitEstimate, 
                                             numberOfModelRuns=n,
                                             randomMethod="calculate",
                                             functionSyntax="data.frameNames")
})
test_that("Individual EVPI is calculated for 4-D partly correlated current estimate.",{
  # Define the model for the profit:
  profitModel <- function(x){
    x$revenue1 + x$revenue2 - (x$costs1 + x$costs2)
  }
  # Create an estimate from text (with correlated components):
  estimateTextMarg<-"variable,  distribution,    lower,   upper
                    revenue1,           norm,      100,    1000
                    revenue2,           norm,       50,    2000
                      costs1,           norm,       50,    2000
                      costs2,           norm,      100,    1000"
  estimateTextCor<-",         revenue1, costs2
                  revenue1,        1,   -0.3
                  costs2,       -0.3,      1"
  profitEstimate<-as.estimate(read.csv(header=TRUE, text=estimateTextMarg,
                                       strip.white=TRUE, stringsAsFactors=FALSE),
                              correlation_matrix=data.matrix(read.csv(text=estimateTextCor,
                                                                      row.names=1,
                                                                      strip.white=TRUE)))
  # Run the EVPI Simulation:
  evpiSimulation <- individualEvpiSimulation(welfare=profitModel, 
                                             currentEstimate=profitEstimate, 
                                             numberOfModelRuns=n,
                                             randomMethod="calculate",
                                             functionSyntax="data.frameNames")
})
test_that("Individual EVPI is calculated for 4-D uncorrelated current estimate
           and 2 dimensional unnamed model function
          (randomMethod=\"calculate\", functionSyntax=\"plainNames\") (1).",{
            # Number of simulations:
            n=10
            # Create the current estimate from text:
            estimateText<-"variable,  distribution, lower, upper
            revenue1,  posnorm,      100,   1000
            revenue2,  posnorm,      50,    2000
            costs1,    posnorm,      50,    2000
            costs2,    posnorm,      100,   1000"
            currentEstimate<-as.estimate(read.csv(header=TRUE, text=estimateText,
                                                  strip.white=TRUE, stringsAsFactors=FALSE))
            # The welfare function:
#            profitModel <- function(x,varnames){
            profitModel <- function(x){
              # Assign the variable names to the function environement:
#               tt<-table(varnames)
#               for (t in 1:length(tt))
#                 assign(names(tt[t]),as.numeric(x[which(varnames==names(tt[t]))]))
              for(i in names(x)) assign(i, as.numeric(x[i]))
              
              list(revenue1 + revenue2 - costs1 - costs2, revenue1)
            }
            # Calculate the Individual EVPI:
            individualEvpiResult<-individualEvpiSimulation(welfare=profitModel,
                                                           currentEstimate=currentEstimate,
                                                           numberOfModelRuns=n,
                                                           functionSyntax="plainNames",
                                                           verbosity=0)
          })
test_that("Example from Hubbard (2014), ch. 7, The value of information for ranges
           is reproduced correctly.",{
             # Number of simulations:
             n=100000
             # Relative comparison tolerance:
             tolerance=4.5/sqrt(n)
             # Example from Hubbard (2014), ch. 7, The value of information for ranges:
             ## welfare function: w_{PA}(x) = px - c with 
             ##    x: units sold,
             ##    P_x: Normal distribution with 90\%-confidence interval: $[c_l,c_u]:=[1.5 \cdot 10^5, 3.0 \cdot 10^5]$
             ##    p: unit price (=25$)
             ##    c: campaign costs
             sales<-estimate("norm", 1.5e+05, 3.0e+05, variable="sales")
             p<-25
             c<-5e+06
             profitModel<-function(x) {
               list(Profit = p*x$sales - c)
             }
             # Run IndividualEVPI simulation:
             individualEvpiResult<-individualEvpiSimulation(welfare=profitModel,
                                                            currentEstimate=sales,
                                                            numberOfModelRuns=n,
                                                            functionSyntax="data.frameNames")
             ## Parameters of the distribution of sales: 
             ##  mu = (c_l + c_u)/2
             ##  sigma = (mu - c_l)/(Phi^{-1}(0.95))
             mu<-(sales$marginal["sales","lower"] + sales$marginal["sales","upper"])/2
             sigma<-(mu - sales$marginal["sales","lower"])/qnorm(0.95)
             ## Calculation of theoretical expectation values:
             ###  EL_{PA}(current)=( c - p*mu) Phi( 1/sigma(c/p - mu) ) 
             ###                    + p*sigma^2 /(sqrt(2*pi)*sigma) e^(-1/(2*sigma^2) (c/p - mu)^2 )
             elPaCurrent_calc<-(c - p*mu)*pnorm( (c/p - mu)/sigma ) + p*sigma^2 * dnorm(x=c/p,mean=mu,sd=sigma)
             ### ENB_PA(current) = p*mu - c
             enbPaCurrent_calc <- p*mu - c
             ###  EL_{SQ}(current)=( p*mu - c ) Phi( 1/sigma(mu - c/p ) ) 
             ###                    + p*sigma^2 /(sqrt(2*pi)*sigma) e^(-1/(2*sigma^2) (c/p - mu)^2 )
             elSqCurrent_calc<-(p*mu - c)*pnorm( (mu - c/p)/sigma ) + p*sigma^2 * dnorm(x=c/p,mean=mu,sd=sigma)
             ### EOL(current)=EL_PA, c<=p*mu; =EL_SQ, otherwise
             eolCurrent_calc<-ifelse(c<=p*mu, yes=elPaCurrent_calc, no=elSqCurrent_calc)
             ### EL_PA(prospective)=0, c<=p*mu; =c-p*mu, otherwise
             elPaProspective_calc<-ifelse(c<=p*mu, yes=0, no=c-p*mu)
             ### ENB_PA(prospective) = ENB_PA(current) 
             enbPaProspective_calc <- enbPaCurrent_calc
             ### EL_SQ(prospective)=p*mu-c, c<=p*mu; =0, otherwise
             elSqProspective_calc<-ifelse(c<=p*mu, yes=p*mu-c, no=0)
             ### EOL(prospective)=0
             eolProspective_calc<-0
             ### IndividualEVPI[x=mu]=EOL(current)
             IndividualEvpiXMu_calc<-eolCurrent_calc
             ## Safe the corresponding simulated values:
             elPaCurrent_sim <-individualEvpiResult$current$elPa[["Profit"]]
             enbPaCurrent_sim<-individualEvpiResult$current$enbPa[["Profit"]]
             elSqCurrent_sim <-individualEvpiResult$current$elSq[["Profit"]]
             eolCurrent_sim  <-individualEvpiResult$current$eol[["Profit"]]
             elPaProspective_sim <-individualEvpiResult$prospective$sales$elPa[["Profit"]]
             enbPaProspective_sim<-individualEvpiResult$prospective$sales$enbPa[["Profit"]]
             elSqProspective_sim <-individualEvpiResult$prospective$sales$elSq[["Profit"]]
             eolProspective_sim  <-individualEvpiResult$prospective$sales$eol[["Profit"]]
             IndividualEvpiXMu_sim<-individualEvpiResult$evi["Profit","sales"]
             ## Compare theoretical with simulated values:
             expect_equal(elPaCurrent_sim,  elPaCurrent_calc,  tolerance=tolerance, scale=elPaCurrent_calc,  use.names=FALSE)
             expect_equal(enbPaCurrent_sim, enbPaCurrent_calc, tolerance=tolerance, scale=enbPaCurrent_calc, use.names=FALSE)
             expect_equal(elSqCurrent_sim,  elSqCurrent_calc,  tolerance=tolerance, scale=elSqCurrent_calc,  use.names=FALSE)
             expect_equal(eolCurrent_sim,   eolCurrent_calc,   tolerance=tolerance, scale=eolCurrent_calc,   use.names=FALSE)
             expect_equal(elPaProspective_sim,  elPaProspective_calc,  tolerance=tolerance, scale=elPaProspective_calc,  use.names=FALSE)
             expect_equal(enbPaProspective_sim, enbPaProspective_calc, tolerance=tolerance, scale=enbPaProspective_calc, use.names=FALSE)
             expect_equal(elSqProspective_sim,  elSqProspective_calc,  tolerance=tolerance, scale=elSqProspective_calc,  use.names=FALSE)
             expect_equal(eolProspective_sim,  eolProspective_calc,  tolerance=tolerance, scale=eolProspective_calc,  use.names=FALSE)
             expect_equal(IndividualEvpiXMu_sim,  IndividualEvpiXMu_calc,  tolerance=tolerance, scale=IndividualEvpiXMu_calc,  use.names=FALSE)
           })
