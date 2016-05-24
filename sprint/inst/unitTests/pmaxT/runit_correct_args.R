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

test.correct_args <- function() {

    size_of_rows <- 10

    # Load test data and class label
    data(golub)
    classlabel_1 <- golub.cl
    classlabel_2 <- rep(c(0,1),19)
    classlabel_3 <- rep(0:18,2)

    maxT_test <- c("t", "t.equalvar", "wilcoxon", "f", "pairt", "blockf")
    maxT_side <- c("abs", "upper", "lower")
    maxT_fixed <- c("y", "n")
    maxT_nonpara <- c("y", "n")

    all_options <- .comp_combinations(maxT_test, maxT_side, maxT_fixed, maxT_nonpara)

    temp_output_sink <- tempfile(pattern =  "_sink_" , tmpdir = getwd())

    # Suspend quiting on stop
    options(error = expression(NULL))

    # Sink all output
    sink(file=temp_output_sink, append=FALSE)

    # -------------------------------------------------------------------------------------------------------
    #  For loop to check the random permutations generator with all option combinations (both fixed and not)
    # -------------------------------------------------------------------------------------------------------
    for(i in 1:dim(all_options)[1] ) {

        # Get a random set of rows
        smallgd <- golub[sample(1:dim(golub)[1], size_of_rows),]

        if( any(all_options[i,1] == c("t", "f", "wilcoxon", "t.equalvar")) )
            classlabel <- classlabel_1
        if ( all_options[i,1] == "pairt" )
            classlabel <- classlabel_2
        if ( all_options[i,1] == "blockf" )
            classlabel <- classlabel_3

        # Check if we are indeed using the random generator for these tests
        # Default B is 10000
        newB <- sprint:::getmaxB(classlabel, test=all_options[i,1], 10000)
        invisible(checkEqualsNumeric(newB[2], 0))

        # Execute parallel version
        res_from_pmaxT <- pmaxT(smallgd, classlabel, test=all_options[i,1], side=all_options[i,2], fixed.seed.sampling=all_options[i,3], nonpara=all_options[i,4])

        # Execute serial version
        res_from_maxT <- mt.maxT(smallgd, classlabel, test=all_options[i,1], side=all_options[i,2], fixed.seed.sampling=all_options[i,3], nonpara=all_options[i,4])

        invisible(checkEqualsNumeric(res_from_maxT[,2], res_from_pmaxT[,2]))
       
        invisible(checkEqualsNumeric(res_from_maxT[,3], res_from_pmaxT[,3]))

        invisible(checkEqualsNumeric(res_from_maxT[,4], res_from_pmaxT[,4]))
    }

    # -----------------------------------------------------------------------------------------------------------
    #  For loop to check the "complete permutations" generator with all option combinations (both fixed and not)
    # -----------------------------------------------------------------------------------------------------------

    classlabel_1 <- golub.cl[15:33]
    classlabel_2 <- rep(c(0,1), 14)
    classlabel_3 <- rep(c(0:4), 2)

    for(i in 1:dim(all_options)[1] ) {

        if( any(all_options[i,1] == c("t", "f", "wilcoxon", "t.equalvar")) ) {
            classlabel <- classlabel_1
            # Get a random set of rows
            smallgd <- golub[sample(1:dim(golub)[1], size_of_rows), sample(1:dim(golub)[2], 19)]
        }
        if ( all_options[i,1] == "pairt" ) {
            classlabel <- classlabel_2
            # Get a random set of rows
            smallgd <- golub[sample(1:dim(golub)[1], size_of_rows), sample(1:dim(golub)[2], 28)]
        }
        if ( all_options[i,1] == "blockf" ){
            classlabel <- classlabel_3
            # Get a random set of rows
            smallgd <- golub[sample(1:dim(golub)[1], size_of_rows), sample(1:dim(golub)[2], 10)]
        }

        # Check if we are indeed using the random generator for these tests
        newB <- sprint:::getmaxB(classlabel, test=all_options[i,1], B=0)
        invisible(checkEqualsNumeric(newB[2], 1))

        # Execute parallel version
        res_from_pmaxT <- pmaxT(smallgd, classlabel, test=all_options[i,1], side=all_options[i,2], fixed.seed.sampling=all_options[i,3], nonpara=all_options[i,4])

        # Execute serial version
        res_from_maxT <- mt.maxT(smallgd, classlabel, test=all_options[i,1], side=all_options[i,2], fixed.seed.sampling=all_options[i,3], nonpara=all_options[i,4])

        invisible(checkEqualsNumeric(res_from_maxT[,2], res_from_pmaxT[,2]))
       
        invisible(checkEqualsNumeric(res_from_maxT[,3], res_from_pmaxT[,3]))

        invisible(checkEqualsNumeric(res_from_maxT[,4], res_from_pmaxT[,4]))
    }


    # Remove sink
    sink(file=NULL)

    # Delete sink file
    unlink(temp_output_sink)

    # Enable stop functionality
    options(error = NULL)

}

# = =============================================================== =
# =  Function to return all possible combinations of the arguments  =
# =  passed as input parameters.                                    =
# = =============================================================== =
.comp_combinations <- function(..., depth=NA)
{
    in_args <- list(...)

    if ( is.na(depth) ) {
        depth <- 1
    }

    if ( depth == length(in_args) ) {
        res <- in_args[[depth]]
        total <- 1
        for( i in 1:(depth-1) )
            total <- total * length( in_args[[i]] )
        return(rep(res, total))
    }
    else {
        res <- Recall(..., depth=depth+1)
        temp_vec <- in_args[[depth]]

        total_after <- 1
        for( i in (depth+1):length(in_args) )
            total_after <- total_after * length( in_args[[i]] )

        temp_1 <- NULL
        for( i in 1:length(temp_vec) ) {
            temp_1 <- c(temp_1, as.vector( rep(temp_vec[i], total_after) ))
        }

        if ( depth != 1 ) {
            total <- 1
            for( i in 1:(depth-1) )
                total <- total * length( in_args[[i]] )
            temp_1 <- rep(temp_1, total)
        } else {
            temp = cbind(temp_1, res);
            colnames(temp) <- NULL
            return(temp)
        }

        return(cbind(temp_1, res))
    }
}

