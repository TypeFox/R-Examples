###############################################################################
##
## MCResultAnalytical.R
##
## Definition of class MCResultAnalytical
## Class of mcreg result objects that contain analytical results.
##
## Copyright (C) 2011 Roche Diagnostics GmbH
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

###############################################################################
## Class definition
###############################################################################

setClass(
    Class = "MCResultAnalytical",
    representation = representation(
        ## global mean (weighted mean)
        xmean = "numeric"),
        contains = "MCResult"
)


###############################################################################
## Method registration
###############################################################################

setMethod(f="initialize",signature="MCResultAnalytical",definition=MCResultAnalytical.initialize)

setMethod("calcResponse",signature=c(.Object="MCResultAnalytical"),definition=MCResultAnalytical.calcResponse)

setMethod(f="printSummary",signature=c(.Object="MCResultAnalytical"),definition=MCResultAnalytical.printSummary)
