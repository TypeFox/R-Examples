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

"xpose.calculate.cwres" <-
  function(object,
           cwres.table.prefix="cwtab",
           tab.suffix="",
           sim.suffix="sim",
           est.tab.suffix=".est",
           deriv.tab.suffix=".deriv",
           old.file.convention=FALSE,
           id="ALL",
           printToOutfile=TRUE,
           onlyNonZero=FALSE,
           classic=FALSE,
           ...) {
    
    
    cwres <- compute.cwres(run.number=object@Runno,
                           tab.prefix=cwres.table.prefix,
                           est.tab.suffix=est.tab.suffix,
                           deriv.tab.suffix=deriv.tab.suffix,
                           old.file.convention=old.file.convention,
                           id=id,
                           printToOutfile=printToOutfile,
                           onlyNonZero=onlyNonZero,
                           sim.suffix="",
                           ...)

    if (!is.null(cwres)){
      object@Data$CWRES=as.vector(cwres)
      object@Prefs@Xvardef$cwres="CWRES"    
    }
    
    
    ## simulation CWRES
    if(!is.null(object@Nsim)){
      if(object@Nsim==1){
        cwres.sim <- compute.cwres(run.number=object@Runno,
                                   tab.prefix=cwres.table.prefix,
                                   est.tab.suffix=est.tab.suffix,
                                   deriv.tab.suffix=deriv.tab.suffix,
                                   old.file.convention=old.file.convention,
                                   id=id,
                                   printToOutfile=printToOutfile,
                                   onlyNonZero=onlyNonZero,
                                   sim.suffix=sim.suffix,
                                   ...)

        if (!is.null(cwres.sim)){
          object@SData$CWRES=as.vector(cwres.sim)
          object@Prefs@Xvardef$cwres="CWRES"    
        } else {
          cat("Table file needed to compute CWRES not present\n")
          cat("For simulated data\n\n")
          cat("CWRES not calculated for simulated data\n\n")
          #cat("Simulated data will be dropped from xpose object.\n\n")
        }
      } else {
        cat("CWRES cannot be calculated for table files\n")
        cat("with multiple simulations in a single file\n\n")
        cat("CWRES not calculated for simulated data\n\n")
        #cat("Simulated data will be dropped from xpose object.\n\n")
      }
    }
    

        
    if (classic==TRUE) {
      #.cur.db <- object
      c1<-call("assign",pos = 1, ".cur.db",object)
      eval(c1)
    }
    
    invisible(object)

  }

