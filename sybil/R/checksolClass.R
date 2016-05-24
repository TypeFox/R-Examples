#  checksolClass.R
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


# checksolClass


#------------------------------------------------------------------------------#
#                   definition of the class checksol                           #
#------------------------------------------------------------------------------#

setClass("checksol",
         representation(
              num_of_prob    = "integer",
              exit_code      = "integer",
              exit_num       = "integer",
              exit_meaning   = "character",
              status_code    = "integer",
              status_num     = "integer",
              status_meaning = "character"
        )
)


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

checksol <- function() {
    new("checksol")
}


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# num_of_prob
setMethod("num_of_prob", signature(object = "checksol"),
          function(object) {
              return(object@num_of_prob)
          }
)

setReplaceMethod("num_of_prob", signature = (object = "checksol"),
                 function(object, value) {
                     object@num_of_prob <- value
                     return(object)
                 }
)


# exit_code
setMethod("exit_code", signature(object = "checksol"),
          function(object) {
              return(object@exit_code)
          }
)

setReplaceMethod("exit_code", signature = (object = "checksol"),
                 function(object, value) {
                     object@exit_code <- value
                     return(object)
                 }
)


# exit_num
setMethod("exit_num", signature(object = "checksol"),
          function(object) {
              return(object@exit_num)
          }
)

setReplaceMethod("exit_num", signature = (object = "checksol"),
                 function(object, value) {
                     object@exit_num <- value
                     return(object)
                 }
)


# exit_meaning
setMethod("exit_meaning", signature(object = "checksol"),
          function(object) {
              return(object@exit_meaning)
          }
)

setReplaceMethod("exit_meaning", signature = (object = "checksol"),
                 function(object, value) {
                     object@exit_meaning <- value
                     return(object)
                 }
)


# status_code
setMethod("status_code", signature(object = "checksol"),
          function(object) {
              return(object@status_code)
          }
)

setReplaceMethod("status_code", signature = (object = "checksol"),
                 function(object, value) {
                     object@status_code <- value
                     return(object)
                 }
)


# status_num
setMethod("status_num", signature(object = "checksol"),
          function(object) {
              return(object@status_num)
          }
)

setReplaceMethod("status_num", signature = (object = "checksol"),
                 function(object, value) {
                     object@status_num <- value
                     return(object)
                 }
)



# status_meaning
setMethod("status_meaning", signature(object = "checksol"),
          function(object) {
              return(object@status_meaning)
          }
)

setReplaceMethod("status_meaning", signature = (object = "checksol"),
                 function(object, value) {
                     object@status_meaning <- value
                     return(object)
                 }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

setMethod("show", signature(object = "checksol"),
    function(object) {
        cat("Return code:\n")
        cat(" Code    #       meaning\n")
        tmp <- sprintf(" %-8i%-8i%s\n",
                       exit_code(object),
                       exit_num(object),
                       exit_meaning(object))
        cat(tmp, sep = "")

        cat("\n")
        cat("Solution status:\n")
        cat(" Code    #       meaning\n")

        tmp <- sprintf(" %-8i%-8i%s\n",
                       status_code(object),
                       status_num(object),
                       status_meaning(object))
        cat(tmp, sep = "")
        num <- num_of_prob(object)
        if (num > 1) {
            cat("\n", num, " optimizations were performed.\n", sep = "")
        }
    }
)
