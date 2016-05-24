#
#    _____      _ ____  ____
#   / ___/_____(_) __ \/ __ )
#   \__ \/ ___/ / / / / __  |
#  ___/ / /__/ / /_/ / /_/ / 
# /____/\___/_/_____/_____/  
#
#
#
# BEGIN_COPYRIGHT
#
# This file is part of SciDB.
# Copyright (C) 2008-2014 SciDB, Inc.
#
# SciDB is free software: you can redistribute it and/or modify
# it under the terms of the AFFERO GNU General Public License as published by
# the Free Software Foundation.
#
# SciDB is distributed "AS-IS" AND WITHOUT ANY WARRANTY OF ANY KIND,
# INCLUDING ANY IMPLIED WARRANTY OF MERCHANTABILITY,
# NON-INFRINGEMENT, OR FITNESS FOR A PARTICULAR PURPOSE. See
# the AFFERO GNU General Public License for the complete license terms.
#
# You should have received a copy of the AFFERO GNU General Public License
# along with SciDB.  If not, see <http://www.gnu.org/licenses/agpl-3.0.html>
#
# END_COPYRIGHT
#

# A general SciDB array class for R. It's a hybrid S4 class with some S3
# methods. The class can represent SciDB arrays and array promises.
#
# A scidb object is fully defined by:
# name = any SciDB expression that can produce an array 
# schema = the corresponding SciDB array schema for 'name' above
# gc = environment
#      If gc$remove = TRUE, remove SciDB array when R gc is run on object.
#      The gc environment also stores dependencies required by array promises.
#
# Objects includes a few additional convenience slots:
# dimensions = a character vector of SciDB dimension names
# attributes = character vector of array attribute names
#              (derived from the schema)

setClassUnion("numericOrNULL", c("numeric", "NULL")) 
setClass("scidb",
         representation(name="character",
                        schema="character",
                        attributes="character",
                        dimensions="character",
                        logical_plan="character",
                        gc="environment"),
         S3methods=TRUE)

setClass("scidbdf",
         representation(name="character",
                        schema="character",
                        attributes="character",
                        dimensions="character",
                        logical_plan="character",
                        gc="environment"),
         S3methods=TRUE)

setClassUnion("scidb_or_scidbdf", c("scidb", "scidbdf")) 
setClassUnion("MNSN", c("missing", "NULL", "scidb", "scidbdf", "numeric"))
