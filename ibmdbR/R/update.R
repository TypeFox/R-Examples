# 
# Copyright (c) 2013, 2014, IBM Corp. All rights reserved. 
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

idaUpdate <- function (db2Conn, updf, dfrm, idaIndex = "", conType = "odbc") {
  if (conType == "odbc") {
    if (idaIndex != ""){
      sqlUpdate(db2Conn, dfrm, updf, index = idaIndex,fast=ifelse(idaIsOracleMode(), F, T))
    } else {
      sqlUpdate(db2Conn, dfrm, updf, fast = ifelse(idaIsOracleMode(), F, T))
    }
  }
}