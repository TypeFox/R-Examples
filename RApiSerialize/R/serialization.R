
##  RApiSerialize -- Packge to provide Serialization as in the R API 
##
##  Copyright (C) 2014         Dirk Eddelbuettel 
##
##  This file is part of RApiSerialize.
##
##  RApiSerialize is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 2 of the License, or
##  (at your option) any later version.
##
##  RApiSerialize is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with RApiSerialize.  If not, see <http://www.gnu.org/licenses/>.

serializeToRaw <- function(obj) {
    .Call("serializeToRaw", obj, PACKAGE="RApiSerialize")
}

unserializeFromRaw <- function(obj) {
    .Call("unserializeFromRaw", obj, PACKAGE="RApiSerialize")
}
