# Part of the mi package for multiple imputation of missing data
# Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015 Trustees of Columbia University
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


setGeneric("change", def = function(data, y, to, what, ...) standardGeneric("change"))
setGeneric("change_family", def = function(data, y, to, ...) standardGeneric("change_family"))
setGeneric("change_imputation_method", def = function(data, y, to, ...) standardGeneric("change_imputation_method"))
setGeneric("change_link", def = function(data, y, to, ...) standardGeneric("change_link"))
setGeneric("change_model", def = function(data, y, to, ...) standardGeneric("change_model"))
setGeneric("change_size", def = function(data, y, to, ...) standardGeneric("change_size"))
setGeneric("change_transformation", def = function(data, y, to, ...) standardGeneric("change_transformation"))
setGeneric("change_type", def = function(data, y, to, ...) standardGeneric("change_type"))
setGeneric("complete", def = function(y, m, ...) standardGeneric("complete"))
setGeneric("fit_model", def = function(y, data, ...) standardGeneric("fit_model"))
setGeneric("get_parameters", def = function(object, ...) standardGeneric("get_parameters"))
setGeneric("hist", def = function(x, ...) standardGeneric("hist"))
setGeneric("mi", def = function(y, model, ...) standardGeneric("mi"))
setGeneric("missing_variable", def = function(y, type, ...) standardGeneric("missing_variable"))
setGeneric("missing_data.frame", def = function(y, ...) standardGeneric("missing_data.frame"))
## FIXME: acount for the other stuff in the original AllGeneric.R
