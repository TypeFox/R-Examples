##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
#  Copyright © 2008,2010 The University of Edinburgh                     #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  any later version.                                                    #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program. If not, see <http://www.gnu.or/licenses/>.   #
#                                                                        #
##########################################################################

#================= Sample data =================

matrixx = matrix(sin(1:2000), ncol=50)

big_matrixx = matrix(sin(1:10000), ncol=50)

listt = list(matrix(sin(1:10000), ncol=200), matrix(sin(1:200), ncol=25), matrix(sin(1:50), ncol=5))

ffobjectt = ff(sin(1:10000), vmode="double", dim=c(200,50))

simple_list = list(1:10,3:13)  # Does not work.

logic_list = list(c(TRUE,FALSE,FALSE,TRUE)) # Does not work.

integers = 1:10  # Does not work.

filename_ = "test.out"

#===============================================

#Compare results for mean function applied over rows
test.compare_results_matrix_rows <- function()
  {

    expected_result = apply(matrixx, 1, mean)
    papply_result = papply(matrixx, mean, 1)   

    checkEquals(as.vector(expected_result), unlist(papply_result), " papply mean function on matrix over rows")
    
  }

#Compare results for mean ff function applied over rows
test.compare_results_ff_matrix_rows <- function()
{
	DEACTIVATED()
    expected_result = apply(big_matrixx, 1, mean)
    papply_result = papply(ffobjectt, mean, 1, out_filename=filename_)
	
    checkEquals(as.vector(expected_result), papply_result[], " papply mean function on ff matrix over rows")
    
}

#Compare results for mean function applied over columns
test.compare_results_matrix_columns <- function()
  {
    
    expected_result = apply(matrixx, 2, mean)
    papply_result = papply(matrixx, mean, 2)   

    checkEquals(as.vector(expected_result), unlist(papply_result), " papply mean function on matrix over columns")
    
  }

#Compare results for mean function applied over columns to ff matrix
test.compare_results_ff_matrix_columns <- function()
{
    DEACTIVATED()
    expected_result = apply(big_matrixx, 2, mean)
	papply_result <- papply(ffobjectt, mean, 2, out_filename=filename_)
	
    checkEquals(as.vector(expected_result), papply_result[], " papply mean function on ff matrix over columns")
    
}

#Compare results for sin function applied over rows
test.compare_results_matrix_rows_sin <- function()
  {

    expected_result = apply(matrixx, 1, sin)
    papply_result = papply(matrixx, sin, 1)

    checkEquals(as.vector(expected_result), unlist(papply_result), " papply sin function on matrix over rows")
    
  }

#Compare results for sin ff function applied over rows
test.compare_results_ff_matrix_rows_sin <- function()
{
	DEACTIVATED()
    expected_result = apply(big_matrixx, 1, sin)
    papply_result = papply(ffobjectt, sin, 1, out_filename=filename_)
	
    checkEquals(as.vector(expected_result), papply_result[], " papply sin function on ff matrix over rows")
    
}

#Compare results for sin function applied over columns
test.compare_results_matrix_columns_sin <- function()
  {
    
    expected_result = apply(matrixx, 2, sin)
    papply_result = papply(matrixx, sin, 2)   

    checkEquals(as.vector(expected_result), unlist(papply_result), " papply sin function on matrix over columns")
    
  }

#Check result when a function definition is passed as an argument
test.compare_results_function_definition_line <- function()
  {

    expected_result = apply(matrixx, 1, function(x) return(x*2))
    papply_result = papply(matrixx, function(x) return(x*2))

    checkEquals(as.vector(expected_result), unlist(papply_result), " papply function definition (line) as an argument")
    
  }

#Check result when a function definition is passed as an argument (multiple lines)
test.compare_results_function_definition_multiline <- function()
  {

    expected_result = apply(matrixx, 2, function(x)
      {
      return(x*2)
      }
      )
    papply_result =  papply(matrixx, function(x)
      {
      return(x*2)
      }, 2
      )

    checkEquals(as.vector(expected_result), unlist(papply_result), " papply function definition (multiple lines) as an argument") 
    
  }

#==== Test parallel lapply ====#

test.compare_results_list <- function()
  {

    expected_result = lapply(listt, mean)
    papply_result = papply(listt,mean)

    checkEquals(expected_result, papply_result, " papply, list of matrices as a data argument, mean function")

  }
test.error_message_for_ffobject <- function()
{
	expected_message = "papply only supports a matrix or a list of matrices as input."
	checkException(papply(ffobjectt, mean), "papply, ff object, mean function")
	checkTrue(as.logical(grep(expected_message, geterrmessage())), "Expected error message when ff object passed into papply")
}


test.compare_results_simple_list <- function()
{
	expected_message = "papply only supports a matrix or a list of matrices as input."
	checkException(papply(simple_list, mean), "papply, integers as a data argument, mean function")
	checkTrue(as.logical(grep(expected_message, geterrmessage())), "Expected error message when integers passed into papply")
}

test.compare_results_logic_list <- function()
{
	expected_message = "papply only supports a matrix or a list of matrices as input."
	checkException(papply(logic_list, mean), "papply, integers as a data argument, mean function")
	checkTrue(as.logical(grep(expected_message, geterrmessage())), "Expected error message when integers passed into papply")
}

test.compare_results_integers <- function()
{
	expected_message = "papply only supports a matrix or a list of matrices as input."
	checkException(papply(integers,mean), "papply, integers as a data argument, mean function")
	checkTrue(as.logical(grep(expected_message, geterrmessage())), "Expected error message when integers passed into papply")
}
