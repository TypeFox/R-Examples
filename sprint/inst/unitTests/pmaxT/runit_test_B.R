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


# ================================================================
#  Run the parallel and serial version of maxT and compare the   =
#  results. The test sets are creted by randomly selecting rows  =
#  from the "golub" data set provided by the multtest package    =
# ================================================================

test.check_B_comp <- function() {

    # Sink error messages. I know it's gonna yield but don't need to see it
#    temp_output_sink <- tempfile(pattern =  "_sink_" , tmpdir = getwd())

    # Suspend quiting on stop
    options(error = expression(NULL))

    # Load test data and class label
    data(golub)
    classlabel <- golub.cl

    # * *************************************************************** *
    # *  Check these tests first, they have the same computation for B  *
    # * *************************************************************** *
    all_test <- c("t", "f", "wilcoxon", "t.equalvar")

    for( next_test in all_test) {

        # All these checks should cause getmaxB to fail because of restrictions
        # Test for: (1) zero (which means complete permutations)
        #           (2) above the maximum limit
        for( i in c(0, sprint:::.BLIM+1) ) {
            invisible(suppressWarnings(newB <- sprint:::getmaxB(classlabel, next_test, i)))
            # Check that the value of B is *not* equal to 0
            invisible(checkTrue(is.na(newB[1])))
            invisible(checkTrue(is.na(newB[2])))
        }

        # All these checks should cause getmaxB to return a valid value
        # Test for: (1) just one permutation
        #           (2) a random acceptable number
        #           (3) a random acceptable number
        #           (4) a random acceptable number
        #           (5) exactly as many as the maximum limit
        for( i in c(1, 12485, 34871, 70835, sprint:::.BLIM) ) {
            newB <- sprint:::getmaxB(classlabel, next_test, i)
            # Check that the value of B is *not* equal to 0
            invisible(checkTrue(!is.na(newB[1])))
            invisible(checkTrue(!is.na(newB[2])))
        }
        # Additional check with a classlabel small enough to support complete permutations
        newB <- sprint:::getmaxB(classlabel[23:32], next_test, i)
        invisible(checkTrue(!is.na(newB[1])))
        invisible(checkTrue(!is.na(newB[2])))

    }

    # * *************** *
    # *  Check "pairt"  *
    # * *************** *
    all_test <- c("pairt")
    classlabel<-rep(c(0,1),19)
    classlabel_ext <- c(classlabel, classlabel)
    for( next_test in all_test) {

        # All these checks should cause getmaxB to fail because of restrictions
        # Test for: (1) zero (which means complete permutationsn not doable with the classlabel tested)
        #           (2) above the maximum limit
        for( i in c(0, sprint:::.BLIM+1) ) {
            invisible(suppressWarnings(newB <- sprint:::getmaxB(classlabel_ext, next_test, i)))
            # Check that the value of B is NA
            invisible(checkTrue(is.na(newB[1])))
            invisible(checkTrue(is.na(newB[2])))
        }

        # All these checks should cause getmaxB to return a valid value
        # Test for: (1) complete permutations (doable for the classlabel tested)
        #           (2) just one permutations
        #           (3) a random acceptable number
        #           (4) a random acceptable number
        #           (5) exactly as many as the maximum limit
        for( i in c(0, 1, 10865, 23375, sprint:::.BLIM) ) {
            newB <- sprint:::getmaxB(classlabel, next_test, i)
            # Check that the value of B is *not* equal to NA
            invisible(checkTrue(!is.na(newB[1])))
            invisible(checkTrue(!is.na(newB[2])))
        }
    }

    # * **************** *
    # *  Check "blockf"  *
    # * **************** *
    all_test <- c("blockf")
    classlabel <- rep(0:18,2)
    classlabel_small <- rep(0:4,2)
    # Run for different values of B between 0 and the maximum allowed limit
    for( next_test in all_test) {

        # All these checks should cause getmaxB to fail because of restrictions
        # Test for: (1) zero (which means complete permutations)
        #           (2) above the maximum limit
        for( i in c(0, sprint:::.BLIM+1) ) {
            invisible(suppressWarnings(newB <- sprint:::getmaxB(classlabel, next_test, i)))
            # Check that the value of B is *not* equal to 0
            invisible(checkTrue(is.na(newB[1])))
            invisible(checkTrue(is.na(newB[2])))
        }

        # All these checks should cause getmaxB to return a valid value
        # Test for: (1) complete permutations (doable for the classlabel tested)
        #           (2) just one permutations
        #           (3) a random acceptable number
        #           (4) a random acceptable number
        #           (5) exactly as many as the maximum limit
        for( i in c(0, 1, 45639, 24183, sprint:::.BLIM) ) {
            newB <- sprint:::getmaxB(classlabel_small, next_test, i)
            # Check that the value of B is *not* equal to 0
            invisible(checkTrue(!is.na(newB[1])))
            invisible(checkTrue(!is.na(newB[2])))
        }
    }


    # Delete sink file
#    unlink(temp_output_sink)

    # Enable stop functionality
    options(error = NULL, warning = NULL)

}

