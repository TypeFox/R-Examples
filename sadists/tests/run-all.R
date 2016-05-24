# /usr/bin/r
#
# Copyright 2015-2016 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav
#
# This file is part of sadists.
#
# sadists is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# sadists is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with sadists.  If not, see <http://www.gnu.org/licenses/>.
#
# Created: 2016.01.04
# Copyright: Steven E. Pav, 2016
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

# because Hadley says it should be like this.
# see https://github.com/hadley/devtools/wiki/Testing

library(testthat)
library(sadists)

#test_package("sadists")
test_check("sadists")
