# 
# Copyright (c) 2010, 2014, IBM Corp. All rights reserved. 
# 		
# This program is free software: you can redistribute it and/or modify 
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or 
# (at your option) any later version. 
#
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU General Public License for more details. 
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>. 
#

#remote data frame
setClass("ida.data.frame", representation(table="character",where="character",cols="character",colDefs="list"))
setClass("ida.data.frame.rows", representation(where="character"))

#object storage
setClass("ida.list", representation(tableName="character"));

#Column expressions
setClass("ida.col.def", representation(term="character",table="ida.data.frame",type="character",aggType="character"))
