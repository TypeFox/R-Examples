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

vpc.npc.menu <- function() {

  choices <- c("Return to previous menu ->",
               "Numerical predictive check plot",
               "Visual predictive check (VPC) plot",
               "Categorical VPC plot",
               "Categorical and continuous VPC plot",
               "* Settings"
               )

  title="\nVISUAL AND NUMERICAL PREDICTIVE CHECK PLOTS MENU\n  \\main\\Visual and numerical predictive check plots"

  pick <- menu(choices,title=title)

  
  run.npc.coverage <- function(){
    cat("\nPlease type the name of the npc results file from PsN\n",
        "Relative or full paths to the file may be used:\n")

    ans <- readline()
    npc.info <- as.character(ans)
    cat("\nRunning command:\n",
        "  npc.coverage(npc.info=\"",npc.info,"\")\n",sep="")
    print(npc.coverage(npc.info=npc.info))
  }

  run.xpose.VPC <- function(){
    cat("\nPlease type the name of the vpc results file from PsN\n",
        "Relative or full paths to the file may be used:\n")

    ans <- readline()
    vpc.info <- as.character(ans)

    cat("\nPlease type the name of the vpctab file from PsN\n",
        "Relative or full paths to the file may be used:\n")

    ans <- readline()
    vpctab <- as.character(ans)
    cat("\nRunning command:\n",
        "  xpose.VPC(vpc.info=\"",vpc.info,"\", vpctab=\"",vpctab,"\")\n",sep="")
    print(xpose.VPC(vpc.info=vpc.info,vpctab=vpctab))
  }

  run.xpose.VPC.categorical <- function(){
    cat("\nPlease type the name of the vpc results file from PsN\n",
        "Relative or full paths to the file may be used:\n")

    ans <- readline()
    vpc.info <- as.character(ans)

    cat("\nPlease type the name of the vpctab file from PsN\n",
        "Relative or full paths to the file may be used:\n")

    ans <- readline()
    vpctab <- as.character(ans)
    cat("\nRunning command:\n",
        "  xpose.VPC.categorical(vpc.info=\"",vpc.info,"\", vpctab=\"",vpctab,"\")\n",sep="")
    print(xpose.VPC.categorical(vpc.info=vpc.info,vpctab=vpctab))
  }

  run.xpose.VPC.both <- function(){
    cat("\nPlease type the name of the vpc results file from PsN\n",
        "Relative or full paths to the file may be used:\n")

    ans <- readline()
    vpc.info <- as.character(ans)

    cat("\nPlease type the name of the vpctab file from PsN\n",
        "Relative or full paths to the file may be used:\n")

    ans <- readline()
    vpctab <- as.character(ans)
    cat("\nRunning command:\n",
        "  xpose.VPC.both(vpc.info=\"",vpc.info,"\", vpctab=\"",vpctab,"\")\n",sep="")
    print(xpose.VPC.both(vpc.info=vpc.info,vpctab=vpctab))
  }


  qx <- 0
  switch(pick+1,
         qx <- 2,
         qx <- 1,
         run.npc.coverage(),
         run.xpose.VPC(),
         run.xpose.VPC.categorical(),
         run.xpose.VPC.both(),
         cat("Not yet implemented, please use command line for this feature!\n",
             "See '?xpose.VPC' or '?npc.coverage'\n")
         )

  if(qx == 2) {
    return(invisible(2))
  } else {
    if(qx == 1) {
      return(invisible(0))
    } else {
      Recall()
    }
  }


}
