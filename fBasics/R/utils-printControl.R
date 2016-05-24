
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received A copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                 CONTROL ATTRIBUTES:
#  print.control             Prints unlisted control attributes  
# FUNCTION:                 DESCRIPTION:
#  .print                    Used in regression package
################################################################################


print.control <-
function(x, ...)
{
    # Return Value:
    print(unlist(x))
}


################################################################################

.print <- 
function(x, ...)
{
    UseMethod(".print")
}


# ------------------------------------------------------------------------------


.plot <- 
function(x, ...)
{
    UseMethod(".plot")
}


# ------------------------------------------------------------------------------


.summary <- 
function(object, ...)
{
    UseMethod(".summary")
}


# ------------------------------------------------------------------------------


.predict <- 
function(object, ...)
{
    UseMethod(".predict")
}


################################################################################

