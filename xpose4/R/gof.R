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
gof <- function(runno=NULL,save=FALSE,onefile=FALSE,
                saveType="pdf",pageWidth=7.6,pageHeight=4.9,
                structural = TRUE,residual=TRUE,covariate=FALSE,
                iiv=FALSE,iov=FALSE,all=FALSE,myTrace=xpPage) {

  ## This is a template function for creating goodness of fit plots using
  ## Xpose specific. Type ?gof at the R-prompt for more help.

#########################################################
## Start the graphics device and define save-filenames ##
## (probably no need to change this).                  ##
#########################################################
  gofSetup(runno,save,onefile,saveType,pageWidth,pageHeight)
 
#############################################
### Set up the data and create the graphs ###
#############################################

## Note! With lattice it is necessary to issue print(plotobject) to
## make the plot appear on the graphics device., e.g.
## xplot <- xyplot(DV~TIME,data)
## print(xplot)

  ## Read the data and do any modifications
  xpdb <- xpose.data(runno)

  ## Create structural model diagnostics
  if(structural || all) {

    xplot <- dv.vs.pred.ipred(xpdb,page=myTrace)
    print(xplot)
  }

  ## Create residual model diagnostics
  if(residual || all) {

    xplot <- absval.wres.vs.pred(xpdb,page=NULL)
    print(xplot)

  }

  ## Create covariate model diagnostics
  if(covariate || all) {

  }

  ## Create iiv model diagnostics
  if(iiv || all) {

  }

  ## Create iov model diagnostics
  if(iov || all) {

  }


#######################
### Finishing tasks ###
#######################

  ## Turn off device if a file device, i.e. save=TRUE
  if(save) dev.off()
  invisible()
}

