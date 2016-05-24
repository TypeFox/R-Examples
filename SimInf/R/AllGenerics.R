## siminf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015  Stefan Engblom
## Copyright (C) 2015  Stefan Widgren
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

##' Run siminf stochastic simulation algorithm
##'
##' @rdname run-methods
##' @docType methods
##' @param model The siminf model to run.
##' @param threads Number of threads. Default is NULL, i.e. to use the
##' number of available processors.
##' @param seed Random number seed. Default is NULL, i.e. to use a
##' time-seed.
##' @return \code{siminf_model} with result from simulation.
##' @examples
##' ## Create a 'SISe' demo model with 1 node and
##' ## initialize it to run over 1000 days.
##' model <- demo_model(nodes = 1, days = 1000, model = "SISe")
##' run(model)
setGeneric("run",
           signature = "model",
           function(model,
                    threads  = NULL,
                    seed     = NULL) standardGeneric("run"))

##' Susceptible
##'
##' Extracts the number of susceptible
##' @rdname susceptible-methods
##' @docType methods
##' @param model The \code{model} to extract the susceptible from
##' @param ... Additional arguments affecting the measure
##' @param age For models with age categories, the age category to
##' extract.
##' @param i Indices specifying the nodes to include when extracting
##' the number of susceptible. Default is NULL, which includes all
##' nodes.
##' @param by The number to increment the sequence of time points
##' starting from 1. Default is 1, which gives the number of
##' susceptible at every time point.
##' @keywords methods
##' @export
##' @examples
##' ## Create a 'SISe' demo model with 5 nodes and initialize
##' ## it to run over 10 days.
##' model <- demo_model(nodes = 5, days = 10, model = "SISe")
##'
##' ## Run the model and save the result
##' result <- run(model)
##'
##' ## Extract the number of susceptible individuals in each
##' ## node after each time step in the simulation
##' susceptible(result)
##'
##' ## Extract the number of susceptible individuals in the
##' ## first node after each time step in the simulation
##' susceptible(result, i = 1)
##'
##' ## Extract the number of susceptible individuals in the
##' ## first and third node after each time step in the simulation
##' susceptible(result, i = c(1, 3))
##'
##' ## Extract the number of susceptible individuals in the first
##' ## and third node after every other time step in the simulation
##' susceptible(result, i = c(1, 3), by = 2)
##'
##' ## Create a 'SISe3' demo model with 5 nodes and initialize
##' ## it to run over 10 days.
##' model <- demo_model(nodes = 5, days = 10, model = "SISe3")
##'
##' ## Run the model and save the result
##' result <- run(model)
##'
##' ## Extract the sum all of susceptible individuals in all age
##' ## categories in each node after each time step in the simulation
##' susceptible(result)
##'
##' ## Extract the number of susceptible individuals in the first age
##' ## category in each node after each time step in the simulation
##' susceptible(result, age = 1)
##'
##' ## Extract the sum of susceptible individuals in the first and
##' ## second age category in each node after each time step in
##' ## the simulation
##' susceptible(result, age = c(1, 2))
##'
##' ## Extract the number of susceptible individuals in the first age
##' ## category in the first and third node after each time step in
##' ## the simulation
##' susceptible(result, i = c(1, 3), age = 1)
setGeneric("susceptible",
           function(model, ...) standardGeneric("susceptible"))

##' Infected
##'
##' Extracts the number of infected
##' @rdname infected-methods
##' @docType methods
##' @param model The \code{model} to extract the infected from
##' @param ... Additional arguments affecting the measure
##' @param age For models with age categories, the age category to
##' extract.
##' @param i Indices specifying the nodes to include when extracting
##' the number of infected. Default is NULL, which includes all nodes.
##' @param by The number to increment the sequence of time points
##' starting from 1. Default is 1, which gives the number of
##' infected at every time point.
##' @keywords methods
##' @export
##' @examples
##' ## Create a 'SISe' demo model with 5 nodes and initialize
##' ## it to run over 10 days.
##' model <- demo_model(nodes = 5, days = 10, model = "SISe")
##'
##' ## Run the model and save the result
##' result <- run(model)
##'
##' ## Extract the number of infected individuals in each
##' ## node after each time step in the simulation
##' infected(result)
##'
##' ## Extract the number of infected individuals in the
##' ## first node after each time step in the simulation
##' infected(result, i = 1)
##'
##' ## Extract the number of infected individuals in the
##' ## first and third node after each time step in the simulation
##' infected(result, i = c(1, 3))
##'
##' ## Extract the number of infected individuals in the first
##' ## and third node after every other time step in the simulation
##' infected(result, i = c(1, 3), by = 2)
##'
##' ## Create a 'SISe3' demo model with 5 nodes and initialize
##' ## it to run over 10 days.
##' model <- demo_model(nodes = 5, days = 10, model = "SISe3")
##'
##' ## Run the model and save the result
##' result <- run(model)
##'
##' ## Extract the sum all of infected individuals in all age
##' ## categories in each node after each time step in the simulation
##' infected(result)
##'
##' ## Extract the number of infected individuals in the first age
##' ## category in each node after each time step in the simulation
##' infected(result, age = 1)
##'
##' ## Extract the sum of infected individuals in the first and
##' ## second age category in each node after each time step in
##' ## the simulation
##' infected(result, age = c(1, 2))
##'
##' ## Extract the number of infected individuals in the first age
##' ## category in the first and third node after each time step in
##' ## the simulation
##' infected(result, i = c(1, 3), age = 1)
setGeneric("infected",
           function(model, ...) standardGeneric("infected"))

##' Prevalence
##'
##' Calculate the proportion infected individuals
##' @rdname prevalence-methods
##' @docType methods
##' @param model The \code{model} to calculated the prevalence from
##' @param ... Additional arguments affecting the measure
##' @param i Indices specifying the nodes to include in the
##' calculation of the prevalence. If \code{wnp = TRUE}, then
##' specifying which nodes to extract prevalence for. Default is NULL,
##' which includes all nodes.
##' @param age For models with age categories, the age category to
##' include in the calculation. Default is that all age categories are
##' included.
##' @param wnp Determine within-node prevalence. Default is FALSE.
##' @param by The number to increment the sequence of time points
##' starting from 1. Default is 1, which gives the prevalence at every
##' time point.
##' @keywords methods
##' @export
##' @examples
##' ## Create a 'SISe' demo model with 5 nodes and initialize
##' ## it to run over 10 days.
##' model <- demo_model(nodes = 5, days = 10, model = "SISe")
##'
##' ## Run the model and save the result
##' result <- run(model)
##'
##' ## Extract the prevalence of infected nodes after each time
##' ## step in the simulation
##' prevalence(result)
##'
##' ## Extract the prevalence of infected nodes after each time
##' ## step in the simulation when including only the first,
##' ## second and third node in the population at risk.
##' prevalence(result, i = 1:3)
##'
##' ## Extract the prevalence of infected nodes after every other
##' ## time step in the simulation when including only the first,
##' ## second and third node in the population at risk.
##' prevalence(result, i = 1:3, by = 2)
##'
##' ## Extract the within-node prevalence of infected individuals
##' ## in each node after each time step in the simulation
##' prevalence(result, wnp = TRUE)
##'
##' ## Extract the within-node prevalence of infected individuals
##' ## in the first and third node after each time step in the
##' ## simulation
##' prevalence(result, wnp = TRUE, i = c(1, 3))
##'
##' ## Extract the within-node prevalence of infected individuals
##' ## in the first and third node after every other time step in
##' ## the simulation
##' prevalence(result, wnp = TRUE, i = c(1, 3), by = 2)
##'
##' ## Create a 'SISe3' demo model with 5 nodes and initialize
##' ## it to run over 10 days.
##' model <- demo_model(nodes = 5, days = 10, model = "SISe3")
##'
##' ## Run the model and save the result
##' result <- run(model)
##'
##' ## Extract the prevalence of infected nodes after each time
##' ## step in the simulation
##' prevalence(result)
##'
##' ## Extract the within-node prevalence of infected
##' ## individuals in the third age category after each
##' ## time step in the simulation
##' prevalence(result, wnp = TRUE, age = 3)
setGeneric("prevalence",
           function(model, ...) standardGeneric("prevalence"))
