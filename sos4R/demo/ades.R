# Copyright (C) 2011 by 52 North Initiative for Geospatial Open Source Software GmbH, Contact: info@52north.org
# This program is free software; you can redistribute and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation. This program is distributed WITHOUT ANY WARRANTY; even without the implied WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program (see gpl-2.0.txt). If not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA or visit the Free Software Foundation web page, http://www.fsf.org.
# Author: Daniel Nuest (daniel.nuest@uni-muenster.de)
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r
library("sos4R")

##############################################################################
# ADES
ades <- SOS(url = "http://sosades.brgm.fr/REST/sos",
		method = SosSupportedConnectionMethods()[["GET"]])
print(object.size(ades), units=c("Mb"))
# 7.8 Mb, large capabilities file!

ades.off <- sosOfferingIds(ades)

# TODO

###################################
# Demo finished, try another one! #
###################################
