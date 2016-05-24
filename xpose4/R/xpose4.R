# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

xpose4 <- function() {

  
  ## check that classes are present
  if (length(findClass("xpose.data")) < 1) {
    createXposeClasses()
  }


  ##
  ## THIS messes up menu system!  Leave it out!
  ##
  ## Set error handling options such that it always return to the main
  ## menu. This may not necessarily be the best and a more elaborate
  ## error handling may be needed. See help for 'stop', 'try' and
  ## 'invokeRestart'.
  ##oldopts <- options(error=main.menu)
  ##on.exit(options(oldopts))
          
  cat("
              Welcome to Xpose!

       Xpose is a population analysis model
       building aid for NONMEM developed by:

       Andrew C. Hooker, Justin J. Wilkins,
       Mats O. Karlsson and E. Niclas Jonsson

       Pharmacometrics research group, 
       Department of Pharmaceutical Biosciences,
       Uppsala University, Sweden.
		                        
Version: Xpose 4.5.0, 2014-05-16
				                        
http://xpose.sourceforge.net

Please report bugs to Andrew Hooker (andrew.hooker@farmbio.uu.se)!

Xpose, Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
Mats O. Karlsson and E. Niclas Jonsson.  Copyright (C) 2009-2010
Andrew C. Hooker, Mats O. Karlsson and E. Niclas Jonsson.

Xpose is free software and comes with ABSOLUTELY NO WARRANTY.
Xpose is made available under the terms of the GNU Lesser General
Public License (LGPL), version 3 or later. You are welcome to redistribute
it under the conditions described therein.  http://www.gnu.org/licenses/

")
  
  ## Get the data
  change.xp.obj()

  ## Start the menus
  main.menu()
  
  ##return()
}
