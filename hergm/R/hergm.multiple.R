###########################################################################
# Copyright 2009 Nobody                                                   #
#                                                                         #
# This file is part of hergm.                                             #
#                                                                         # 
#    hergm is free software: you can redistribute it and/or modify        #
#    it under the terms of the GNU General Public License as published by #
#    the Free Software Foundation, either version 3 of the License, or    #
#    (at your option) any later version.                                  #
#                                                                         # 
#    hergm is distributed in the hope that it will be useful,             #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
#    GNU General Public License for more details.                         #
#                                                                         #
#    You should have received a copy of the GNU General Public License    #
#    along with hergm.  If not, see <http://www.gnu.org/licenses/>.       #
#                                                                         # 
###########################################################################

hergm.relabel_1 <- function(max_number, indicator, number_runs)
# Relabeling algorithm, which aims to minimize posterior expected loss
# input: number of categories, indicators, number of runs
# output: minimum and minimizer of posterior expected loss 
{
  loss <- vector(length = number_runs)
  minimum_loss <- Inf
  for (i in 1:number_runs)
    {
    if (number_runs > 1) cat("\n------\nRun ", i, "\n------\n", sep="")
    output <- hergm.min_loss_1(max_number, indicator, 100) # Loss function of Schweinberger and Handcock (2015)
    loss[i] <- output$loss
    if (output$loss < minimum_loss)
      {
      minimum_loss <- output$loss
      min_output <- output
      }
    }
  if (number_runs > 1) cat("\n", "Minimum loss: ", min_output$loss, "\n", sep="")
  min_output
}

hergm.relabel_2 <- function(max_number, indicator)
# Relabeling algorithm, which aims to minimize posterior expected loss
# input: number of categories, indicators
# output: minimum and minimizer of posterior expected loss 
{
  min_output <- hergm.min_loss_2(max_number, indicator) # Loss function of Peng and Carvalho (2015); note: the algorithm converges to the same minimum in each run, therefore multiple runs are not necessary
  min_output
}
 
