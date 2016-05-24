################################################################################
# Copyright (C) 2010 by 52 North                                               #
# Initiative for Geospatial Open Source Software GmbH                          #
#                                                                              #
# Contact: Andreas Wytzisk                                                     #
# 52 North Initiative for Geospatial Open Source Software GmbH                 #
# Martin-Luther-King-Weg 24                                                    #
# 48155 Muenster, Germany                                                      #
# info@52north.org                                                             #
#                                                                              #
# This program is free software; you can redistribute and/or modify it under   #
# the terms of the GNU General Public License version 2 as published by the    #
# Free Software Foundation.                                                    #
#                                                                              #
# This program is distributed WITHOUT ANY WARRANTY; even without the implied   #
# WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU #
# General Public License for more details.                                     #
#                                                                              #
# You should have received a copy of the GNU General Public License along with #
# this program (see gpl-2.0.txt). If not, write to the Free Software           #
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA or #
# visit the Free Software Foundation web page, http://www.fsf.org.             #
#                                                                              #
# Author: Daniel Nuest (daniel.nuest@uni-muenster.de)                          #
# Created: 2010-06-18                                                          #
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r #
#                                                                              #
################################################################################

# load required libraries
library("XML")
library("RCurl")
library("sp")

# load required source files for testing
if(!exists(".sos4Rpath"))
	.sos4Rpath = "/home/daniel/Dokumente/2010_SOS4R/workspace/sos4R"
	
source(paste(.sos4Rpath, "R", "Constants.R",  sep = "/"))

source(paste(.sos4Rpath, "R", "Class-OWS.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "Class-GML.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "Class-SWE.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "Class-OM.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "Class-SA.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "Class-OGC.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "Class-SOS.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "Class-SML.R",  sep = "/"))

source(paste(.sos4Rpath, "R", "Generic-methods.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "OWS-methods.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "OWS-methods-parsing.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "SOS-methods-parsing.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "SOS-methods-plotting.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "OM-methods.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "OM-methods-parsing.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "OM-methods-coercion.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "SA-methods.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "GML-methods.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "SWE-methods.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "GML-methods-parsing.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "SA-methods-parsing.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "SWE-methods-parsing.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "OGC-methods.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "SML-methods.R",  sep = "/"))

source(paste(.sos4Rpath, "R", "PrintShowStructureSummmary-methods.R",  sep = "/"))

source(paste(.sos4Rpath, "R", "SOS-methods-coercion.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "SML-methods-coercion.R",  sep = "/"))

source(paste(.sos4Rpath, "R", "SOS-methods-accessor.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "SOS-methods-util.R",  sep = "/"))
source(paste(.sos4Rpath, "R", "SML-methods-util.R",  sep = "/"))

source(paste(.sos4Rpath, "R", "SOS-methods.R",  sep = "/"))

source(paste(.sos4Rpath, "R", "Defaults.R",  sep = "/"))
