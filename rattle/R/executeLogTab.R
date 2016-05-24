# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2014-07-24 21:30:01 gjw>
#
# Execute Log Tab
#
# Copyright (c) 2014 Togaware Pty Ltd
#
# This file is part of Rattle.
#
# Rattle is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# Rattle is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rattle. If not, see <http://www.gnu.org/licenses/>.

executeLogTab <- function()
{
  log.text <- getTextviewContent("log_textview")
  eval(parse(text=log.text))
}

