# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created: pre 01-01-2013
# last modification: 10-10-2014
# Copyright (C) 2014
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
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

print.mht <- print.mht.order <-
function(x,...) {
    # 0. part
    print(x$call)

    # 2a. part
    alpha=as.numeric(colnames(x$coefficients))  	     
    cat("","\n")
    for(i in 1:length(alpha))
{    cat("Results for alpha=",alpha[i],"\n")
    cat("Number of relevant variables:",x$kchap[i],"\n")
  	cat("which are the following columns of object$data$X:",x$relevant_var[i,],"\n\n")

}
}

