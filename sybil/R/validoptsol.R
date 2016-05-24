#  validoptsol.R
#  FBA and friends with R.
#
#  Copyright (C) 2010-2014 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of sybil.
#
#  Sybil is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Sybil is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybil.  If not, see <http://www.gnu.org/licenses/>.


################################################
# Function: .validoptsol
#
# Validity checking of an object of class optsol
#
# Returns TRUE if the model is valid, otherwise
# a character String containing a description of
# the error.


.validoptsol <- function(object) {

    if (!is(object, "optsol")) {
        return("needs an object of class optsol!")
    }
    
#    if (length(solver(object)) != 1) {
#        return("solver must have a length of 1!")
#    }
#    # if ((length(solver(object)) != 1) || (length(method(object)) != 1)) {
#    #     return("solver and method must have a length of 1!")
#    # }
#    
#    if (length(lp_dir(object)) == 0) {
#    
#        if (length(solver(object)) != 1) {
#            return("solver must have a length of 1!")
#        }
#        if (length(method(object)) != 1) {
#            return("method must have a length of 1!")
#        }
#    
#    }
#    else {
#    
#        if (length(lp_dir(object)) != 1) {
#            return("lp_dir must have a length of 1!")
#        }
        if (length(lp_num_cols(object)) != 1) {
            return("lp_num_cols must have a length of 1!")
        }
        if (length(lp_num_rows(object)) != 1) {
            return("lp_num_rows must have a length of 1!")
        }
        
        num_of_prob <- length(lp_obj(object))
        
        if (length(lp_obj(object)) != num_of_prob) {
            return("wrong length of lp_obj!")
        }
        if (length(lp_ok(object)) != num_of_prob) {
            return("wrong length of lp_ok!")
        }
        if (length(lp_stat(object)) != num_of_prob) {
            return("wrong length of lp_stat!")
        }
    
#    }
    return(TRUE)
}
