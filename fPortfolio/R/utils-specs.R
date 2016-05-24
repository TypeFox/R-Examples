
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA
    

####################################################################################################
# FUNCTION:
####################################################################################################
                                                  
                     
# Arguments:          
#                       start objective lower upper linCons funCons amplCons             control ...
#                                          =0    =1                                      =list() 
#                               class()              list()  list()   list()              list()                   
#                                                                                                
# LP Solver                                                                                      
# glpkLP      default       .   numeric     x     x    yes       .        .       glpkLPControl    .
# symphonyLP        .       .   numeric     x     x    yes       .        .   symphonyLPControl    .
# amplLP            .       .   numeric     x     x    yes       .        .       amplLPControl    .
#                                                                                                
# QP Solver                                                                                      
# quadprogQP  default       .      list     x     x    yes       .        .   quadprogQPControl    .
# ipoptQP           .       .      list     x     x    yes       .        .      ipoptQPControl    .
# amplQP            .       .      list     x     x    yes       .        .       amplQPControl    .
#                                                                                                                                                        
# NLP Solver                                                                                     
# solnpNLP    default     yes  function     x     x    yes     yes        .     solnpNLPControl    .
# donlp2NLP         .     yes  function     x     x    yes     yes        .    donlp2NLPControl    .
# nlminb2NLP        .     yes  function     x     x    yes     yes        .   nlminb2NLPControl    .        
# amplNLP           .  length character     x     x      .       .      yes      amplNLPControl  yes      
#


# AMPL Solver tested: ipopt, donlp2, mosek, ...


####################################################################################################


# Value:
#   
#   value <- list(
#       call, control, dots, 
#       opt,  
#       solution, objective, status, message, solver,
#   class(value) <- c("solver", "list")


####################################################################################################


# Methods:

# print.solver <- function(x, ...) x
# summary.solver <- function(x, ...) x
  
# getSolution <- function(x) x$solution
# getObjective <- function(x) x$objective
# getStatus <- function(x) x$status
# getMessage <- function(x) x$message

  
####################################################################################################

