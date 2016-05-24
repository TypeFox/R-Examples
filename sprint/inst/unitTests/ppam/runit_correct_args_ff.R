##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
#  Copyright © 2008,2009 The University of Edinburgh                     #
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


# = =============================================================== =
# =  Massive unit test to check all possible combinations of input  =
# =  parameters and make sure that the output matches the output    =
# =  from the serial version.                                       =
# = =============================================================== =

# === Set up the test suite ====

size_of_rows <- 1000
size_of_columns <- 50
n_clusters <- 12
my_medoids <- NULL

test.ff.matrix <- function() {

    # Generate random data
    input_dataset <- matrix(rnorm(size_of_rows * size_of_columns), ncol=size_of_columns)

    # Create different object types accepted by parallel pam
    data_symmetric_matrix <- 1-cor(t(input_dataset))
    data_distance_matrix <- as.dist(data_symmetric_matrix)
    data_binary_file <- ff(data_symmetric_matrix, vmode="double", dim=c(size_of_rows, size_of_rows))
    
    # Execute original version
    original_pam_result <- pam(data_distance_matrix, n_clusters)

     # Execute parallel version passing different object types on input
    ppam_result_binary_file <- ppam(data_binary_file, n_clusters)

    # Compare clustering results
    
    invisible(checkEqualsNumeric(original_pam_result$clustering,
                                 ppam_result_binary_file$clustering))

    # Compare medoids

    invisible(checkEqualsNumeric(original_pam_result$medoids,
                                 ppam_result_binary_file$medoids))

    # Check clustering information and statistics

    invisible(checkEqualsNumeric(original_pam_result$clusinfo,
                                 ppam_result_binary_file$clusinfo))

    # Check isolation (Isolation, Isolation, Isolation...)

    invisible(checkEqualsNumeric(original_pam_result$isolation,
                                 ppam_result_binary_file$isolation))

    # Check clustering objective

    invisible(checkEqualsNumeric(original_pam_result$objective,
                                 ppam_result_binary_file$objective))
    
}


