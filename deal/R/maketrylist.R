## maketrylist.R
## Author          : Claus Dethlefsen
## Created On      : Fri Jan 11 10:54:00 2002
## Last Modified By: Claus Dethlefsen
## Last Modified On: Wed Jul 23 13:35:14 2003
## Update Count    : 196
## Status          : Unknown, Use with caution!
###############################################################################
##
##    Copyright (C) 2002  Susanne Gammelgaard Bøttcher, Claus Dethlefsen
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program; if not, write to the Free Software
##    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
######################################################################


maketrylist <- function(initnw,data,prior=jointprior(network(data)),
                        timetrace=FALSE) {
    
    if (timetrace) {t1 <- proc.time();cat("[Maketrylist ")}
    
    tryl <- networkfamily(data,initnw,prior,timetrace=timetrace)$trylist

    if (timetrace) {
        t2 <- proc.time()
        cat((t2-t1)[1],"]\n")
    }
    tryl
}


