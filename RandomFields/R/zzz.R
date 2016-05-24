
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  



## sudo apt-get install tcl8.5-dev
## sudo apt-get install tk8.5-dev
## see R Installation and Administration
##./configure --with-tcl-config=/usr/lib/tcl8.5/tclConfig.sh --with-tk-config=/usr/lib/tk8.5/tkConfig.sh
# sudo tar xvf ~/TMP/bwidget-1.9.5.tar 


# readline <- function(...) return("")

#.onLoad <- function(lib, pkg) {
#}

#.onAttach <- function (lib, pkg) {
#  #packageStartupMessage("\n")
#}

#.onUnload <- function(lib, pkg){
#  }
#Implementierung von Cox & Isham's non-separable model

#individuelle Praeferenzliste:
#  1) erste praeferenz
#  2) if # pkte < 10 / 50 -> Gauss??
#  4) if nugget/variance > 0.01 -> CE (and not spectral, e.g.)
#  3) if C(maxh)/C(0) > 0.05  -> TBM else CE


.onUnload <- function(lib, pkg){
  RFoptions(storing=FALSE) ## delete everything
}
 
