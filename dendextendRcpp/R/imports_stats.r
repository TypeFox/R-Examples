# Copyright (C) Tal Galili
#
# This file is part of dendextend.
#
# dendextend is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# dendextend is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

 
# This file includes hidden functions from "stats" that have been imported into this package
# in case their usage would be changed in future versions of R.

## Unexported object imported by a ':::' call: 'stats:::labels.dendrogram'

# stats:::labels.dendrogram
stats_labels.dendrogram <- function (object, ...) {
	unlist(dendrapply(object, function(n) attr(n, "label")))
}


