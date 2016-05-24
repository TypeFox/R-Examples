################################################################################
# Copyright (C) 2013 by 52°North                                               #
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
# Created: 2013-03-06                                                          #
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r #
#                                                                              #
################################################################################

context("utils")

dcps <- c("Post" = "http://url/with/endpoint/one", "Post" = "url.to/endpoint/two",
		"Get" = "some.thing.com/different/")

test_that(".sosFilterDCPs works", {
	expect_that(length(.sosFilterDCPs(dcp = dcps, pattern = "*")), equals(3))
	expect_that(.sosFilterDCPs(dcp = dcps, pattern = list("POX" = "/endpoint"))[[2]],
			is_equivalent_to("url.to/endpoint/two"))
	expect_that(.sosFilterDCPs(dcp = dcps, pattern = list("POX" = "/endpoint")),
			equals(c("Post" = "http://url/with/endpoint/one", "Post" = "url.to/endpoint/two")))
	expect_equivalent(.sosFilterDCPs(dcps, list("POX" = "/one")),
			"http://url/with/endpoint/one")
})
