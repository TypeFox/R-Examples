####**********************************************************************
####**********************************************************************
####
####  SPIKE AND SLAB 1.1.5
####
####  Copyright 2013, Cleveland Clinic Foundation
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 2
####  of the License, or (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  You should have received a copy of the GNU General Public
####  License along with this program; if not, write to the Free
####  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
####  Boston, MA  02110-1301, USA.
####
####  ----------------------------------------------------------------
####  Written and Developed by:
####  ----------------------------------------------------------------
####    Hemant Ishwaran, Ph.D.
####    Director of Statistical Methodology
####    Professor, Division of Biostatistics
####    Clinical Research Building, Room 1058
####    1120 NW 14th Street
####    University of Miami, Miami FL 33136
####
####    email:  hemant.ishwaran@gmail.com
####    URL:    http://web.ccs.miami.edu/~hishwaran
####    --------------------------------------------------------------
####    J. Sunil Rao, Ph. D.
####    Professor and Director of the Division of Biostatistics, 
####    Department of Epidemiology & Public Health
####    Clinical Research Bldg, R-669
####    1120 NW 14th Street, Room 1056
####    Miami, FL 33136
####    email:  rao.jsunil@gmail.com
####    URL:    http://biostat.med.miami.edu/people/primary-faculty/sunil-rao
####  ----------------------------------------------------------------
####  Maintained by:
####    Udaya B. Kogalur, Ph.D.
####    Adjunct Staff
####    Dept of Quantitative Health Sciences
####    Cleveland Clinic Foundation
####    
####    Kogalur & Company, Inc.
####    5425 Nestleway Drive, Suite L1
####    Clemmons, NC 27012
####
####    email:  ubk@kogalur.com
####    URL:    http://www.kogalur.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************

spikeslab.news <- function(...) {
    newsfile <- file.path(system.file(package="spikeslab"), "NEWS")
    file.show(newsfile)
}
