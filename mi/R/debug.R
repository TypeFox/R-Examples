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

## FIXME: Make damn sure .MI_DEBUG is FALSE before pushing to CRAN
.MI_DEBUG <- FALSE 

## FIXME: Also, make all if(.MI_DEBUG) statements one-liners in other files
## and do sed s/if(.MI_DEBUG)/#if(.MI_DEBUG)/g *.R

if(TRUE && .MI_DEBUG) { # define multi-line debugging functions in here
  options(error = recover)
}
