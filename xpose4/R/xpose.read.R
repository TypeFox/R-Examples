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

"xpose.read" <-
  function(object,
           file = "xpose.ini"
           ) {
           
  # x <- object
  options(warn = -1)
  # read ini file
  prefs <- read.table(file, col.names=c("Option","Value"))
  prefs.m <- as.matrix(prefs)
  
  # iterate and assign
  for (i in 1:nrow(prefs.m)) {
  
    # General
    
    if (prefs.m[i,1] == "Miss") {
      object@Prefs@Miss = as.numeric(prefs.m[i,2])
    }
    if (prefs.m[i,1] == "Cat.level") {
      object@Prefs@Cat.level = as.numeric(prefs.m[i,2])
    }
    if (prefs.m[i,1] == "DV.Cat.level") {
      object@Prefs@DV.Cat.level = as.numeric(prefs.m[i,2])
    }
    
    # Plotting 
    if (prefs.m[i,1] == "type") {
      object@Prefs@Graph.prefs$type = as.character(prefs.m[i,2])
    }
    if (prefs.m[i,1] == "cex") {
      object@Prefs@Graph.prefs$cex = as.numeric(prefs.m[i,2])
    }
    if (prefs.m[i,1] == "lty") {
      object@Prefs@Graph.prefs$lty = as.numeric(prefs.m[i,2])
    }   
    if (prefs.m[i,1] == "lwd") {
      object@Prefs@Graph.prefs$lwd = as.numeric(prefs.m[i,2])
    }
    if (prefs.m[i,1] == "col") {
      if ((prefs.m[i,2]!="") && (!is.na(as.numeric(prefs.m[i,2])))) {
        object@Prefs@Graph.prefs$col = as.numeric(prefs.m[i,2])
      } else {
        object@Prefs@Graph.prefs$col = as.character(prefs.m[i,2])
      }
    }
    if (prefs.m[i,1] == "pch") {
      object@Prefs@Graph.prefs$pch = as.numeric(prefs.m[i,2])
    }   
    if (prefs.m[i,1] == "grid") {
      object@Prefs@Graph.prefs$grid = as.logical(prefs.m[i,2])
    }   
    if (prefs.m[i,1] == "aspect") {
      object@Prefs@Graph.prefs$aspect = as.numeric(prefs.m[i,2])
    } 
    
    # Conditioning

#    if (prefs.m[i,1] == "ordby") {
#      if (grep("NULL", prefs.m[i,2]) != 0) {
#        object@Prefs@Graph.prefs$ordby = NULL
#      } else {
#        object@Prefs@Graph.prefs$ordby = as.character(prefs.m[i,2])
#      }
#    }
    if (prefs.m[i,1] == "byordfun") {
      object@Prefs@Graph.prefs$byordfun = as.character(prefs.m[i,2])
    }
    if (prefs.m[i,1] == "shingnum") {
      object@Prefs@Graph.prefs$shingnum = as.numeric(prefs.m[i,2])
    }   
    if (prefs.m[i,1] == "shingol") {
      object@Prefs@Graph.prefs$shingol = as.numeric(prefs.m[i,2])
    }    
        
    # abline
    
#    if (prefs.m[i,1] == "abline") {
#      if (grep("NULL", prefs.m[i,2]) != 0) {
#        object@Prefs@Graph.prefs$abline = NULL
#      } else {
#        object@Prefs@Graph.prefs$abline = as.character(prefs.m[i,2])
#      }
#    }
    if (prefs.m[i,1] == "ablcol") {
      if ((prefs.m[i,2]!="") && (!is.na(as.numeric(prefs.m[i,2])))) {
        object@Prefs@Graph.prefs$ablcol = as.numeric(prefs.m[i,2])
      } else {
        object@Prefs@Graph.prefs$ablcol = as.character(prefs.m[i,2])
      }
    }
    if (prefs.m[i,1] == "abllty") {
      object@Prefs@Graph.prefs$abllty = as.numeric(prefs.m[i,2])
    }   
    if (prefs.m[i,1] == "abllwd") {
      object@Prefs@Graph.prefs$abllwd = as.numeric(prefs.m[i,2])
    }    
    
    # lmline
    
#    if (prefs.m[i,1] == "lmline") {
#      if (grep("NULL", prefs.m[i,2]) != 0) {
#        object@Prefs@Graph.prefs$lmline = NULL
#      } else {
#        object@Prefs@Graph.prefs$lmline = as.character(prefs.m[i,2])
#      }
#    }
    if (prefs.m[i,1] == "lmcol") {
      if ((prefs.m[i,2]!="") && (!is.na(as.numeric(prefs.m[i,2])))) {
        object@Prefs@Graph.prefs$lmcol = as.numeric(prefs.m[i,2])
      } else {
        object@Prefs@Graph.prefs$lmcol = as.character(prefs.m[i,2])
      }
    }
    if (prefs.m[i,1] == "lmlty") {
      object@Prefs@Graph.prefs$lmlty = as.numeric(prefs.m[i,2])
    }   
    if (prefs.m[i,1] == "lmlwd") {
      object@Prefs@Graph.prefs$lmlwd = as.numeric(prefs.m[i,2])
    }    
    
    # smooth
    
    if (prefs.m[i,1] == "smooth") {
      if (grep("NULL", prefs.m[i,2]) != 0) {
        object@Prefs@Graph.prefs$smooth = NULL
      } else {
        object@Prefs@Graph.prefs$smooth = as.character(prefs.m[i,2])
      }
    }
    if (prefs.m[i,1] == "smcol") {
      if ((prefs.m[i,2]!="") && (!is.na(as.numeric(prefs.m[i,2])))) {
        object@Prefs@Graph.prefs$smcol = as.numeric(prefs.m[i,2])
      } else {
        object@Prefs@Graph.prefs$smcol = as.character(prefs.m[i,2])
      }
    }
    if (prefs.m[i,1] == "smlty") {
      object@Prefs@Graph.prefs$smlty = as.numeric(prefs.m[i,2])
    }   
    if (prefs.m[i,1] == "smlwd") {
      object@Prefs@Graph.prefs$smlwd = as.numeric(prefs.m[i,2])
    }  
    if (prefs.m[i,1] == "smspan") {
      object@Prefs@Graph.prefs$smspan = as.numeric(prefs.m[i,2])
    }   
    if (prefs.m[i,1] == "smdegr") {
      object@Prefs@Graph.prefs$smdegr = as.numeric(prefs.m[i,2])
    } 
    
    # suline                          

#    if (prefs.m[i,1] == "suline") {
#      if (grep("NULL", prefs.m[i,2]) != 0) {
#        object@Prefs@Graph.prefs$suline = NULL
#      } else {
#        object@Prefs@Graph.prefs$suline = as.character(prefs.m[i,2])
#      }
#    }
    if (prefs.m[i,1] == "sucol") {
      if ((prefs.m[i,2]!="") && (!is.na(as.numeric(prefs.m[i,2])))) {
        object@Prefs@Graph.prefs$sucol = as.numeric(prefs.m[i,2])
      } else {
        object@Prefs@Graph.prefs$sucol = as.character(prefs.m[i,2])
      }
    }
    if (prefs.m[i,1] == "sulty") {
      object@Prefs@Graph.prefs$sulty = as.numeric(prefs.m[i,2])
    }   
    if (prefs.m[i,1] == "sulwd") {
      object@Prefs@Graph.prefs$sulwd = as.numeric(prefs.m[i,2])
    }  
    if (prefs.m[i,1] == "suspan") {
      object@Prefs@Graph.prefs$suspan = as.numeric(prefs.m[i,2])
    }   
    if (prefs.m[i,1] == "sudegr") {
      object@Prefs@Graph.prefs$sudegr = as.numeric(prefs.m[i,2])
    } 
    
    # Labelling
    
    if (prefs.m[i,1] == "ids") {
      object@Prefs@Graph.prefs$ids = as.logical(prefs.m[i,2])
    }
#    if (prefs.m[i,1] == "idsmode") {
#      if (grep("NULL", prefs.m[i,2]) != 0) {
#        object@Prefs@Graph.prefs$idsmode = NULL
#      } else {
#        object@Prefs@Graph.prefs$idsmode = as.character(prefs.m[i,2])
#      }
#    }
    if (prefs.m[i,1] == "idsext") {
      object@Prefs@Graph.prefs$idsext = as.numeric(prefs.m[i,2])
    }   
    if (prefs.m[i,1] == "idscex") {
      object@Prefs@Graph.prefs$idscex = as.numeric(prefs.m[i,2])
    }  
    if (prefs.m[i,1] == "idsdir") {
      object@Prefs@Graph.prefs$idsdir = as.character(prefs.m[i,2])
    }   
    
    # Dilution
    
    if (prefs.m[i,1] == "dilfrac") {
      object@Prefs@Graph.prefs$dilfrac = as.numeric(prefs.m[i,2])
    }   
#    if (prefs.m[i,1] == "diltype") {
#      if (grep("NULL", prefs.m[i,2]) != 0) {
#        object@Prefs@Graph.prefs$diltype = NULL
#      } else {
#        object@Prefs@Graph.prefs$diltype = as.character(prefs.m[i,2])
#      }
#    }
    if (prefs.m[i,1] == "dilci") {
      object@Prefs@Graph.prefs$dilci = as.numeric(prefs.m[i,2])
    }      
    
    # Prediction intervals

    if (prefs.m[i,1] == "PIuplty") {
      object@Prefs@Graph.prefs$PIuplty = as.numeric(prefs.m[i,2])
    }
    if (prefs.m[i,1] == "PIdolty") {
      object@Prefs@Graph.prefs$PIdolty = as.numeric(prefs.m[i,2])
    }
    if (prefs.m[i,1] == "PImelty") {
      object@Prefs@Graph.prefs$PImelty = as.numeric(prefs.m[i,2])
    }   
    if (prefs.m[i,1] == "PIuptyp") {
      object@Prefs@Graph.prefs$PIuptyp = as.character(prefs.m[i,2])
    }  
    if (prefs.m[i,1] == "PIdotyp") {
      object@Prefs@Graph.prefs$PIdotyp = as.character(prefs.m[i,2])
    }   
    if (prefs.m[i,1] == "PImetyp") {
      object@Prefs@Graph.prefs$PImetyp = as.character(prefs.m[i,2])
    }  
    if (prefs.m[i,1] == "PIupcol") {
      if ((prefs.m[i,2]!="") && (!is.na(as.numeric(prefs.m[i,2])))) {
        object@Prefs@Graph.prefs$PIupcol = as.numeric(prefs.m[i,2])
      } else {
        object@Prefs@Graph.prefs$PIupcol = as.character(prefs.m[i,2])
      }
    }
    if (prefs.m[i,1] == "PIdocol") {
      if ((prefs.m[i,2]!="") && (!is.na(as.numeric(prefs.m[i,2])))) {
        object@Prefs@Graph.prefs$PIdocol = as.numeric(prefs.m[i,2])
      } else {
        object@Prefs@Graph.prefs$PIdocol = as.character(prefs.m[i,2])
      }
    }
    if (prefs.m[i,1] == "PImecol") {
      if ((prefs.m[i,2]!="") && (!is.na(as.numeric(prefs.m[i,2])))) {
        object@Prefs@Graph.prefs$PImecol = as.numeric(prefs.m[i,2])
      } else {
        object@Prefs@Graph.prefs$PImecol = as.character(prefs.m[i,2])
      }
    }   
    if (prefs.m[i,1] == "PIuplwd") {
      object@Prefs@Graph.prefs$PIuplwd = as.numeric(prefs.m[i,2])
    }  
    if (prefs.m[i,1] == "PIdolwd") {
      object@Prefs@Graph.prefs$PIdolwd = as.numeric(prefs.m[i,2])
    }   
    if (prefs.m[i,1] == "PImelwd") {
      object@Prefs@Graph.prefs$PImelwd = as.numeric(prefs.m[i,2])
    }    
    if (prefs.m[i,1] == "PIuplimit") {
      # object@Prefs@Graph.prefs$PIdolwd = prefs.m[i,2]
      PIuplimit = as.numeric(prefs.m[i,2])
    }   
    if (prefs.m[i,1] == "PIdolimit") {
      # object@Prefs@Graph.prefs$PImelwd = prefs.m[i,2]
      PIdolimit = as.numeric(prefs.m[i,2])
    }    
    if (prefs.m[i,1] == "PIarcol") {
      if ((prefs.m[i,2]!="") && (!is.na(as.numeric(prefs.m[i,2])))) {
        object@Prefs@Graph.prefs$PIarcol = as.numeric(prefs.m[i,2])
      } else {
        object@Prefs@Graph.prefs$PIarcol = as.character(prefs.m[i,2])
      }
    }
    
    #B&W plots
    
    if (prefs.m[i,1] == "bwhoriz") {
      object@Prefs@Graph.prefs$bwhoriz = as.logical(prefs.m[i,2])
    }
    if (prefs.m[i,1] == "bwratio") {
      object@Prefs@Graph.prefs$bwratio = as.numeric(prefs.m[i,2])
    }
    if (prefs.m[i,1] == "bwvarwid") {
      object@Prefs@Graph.prefs$bwvarwid = as.logical(prefs.m[i,2])
    }   
    if (prefs.m[i,1] == "bwdotpch") {
      object@Prefs@Graph.prefs$bwdotpch = as.numeric(prefs.m[i,2])
    }  
    if (prefs.m[i,1] == "bwdotcol") {
      if ((prefs.m[i,2]!="") && (!is.na(as.numeric(prefs.m[i,2])))) {
        object@Prefs@Graph.prefs$bwdotcol = as.numeric(prefs.m[i,2])
      } else {
        object@Prefs@Graph.prefs$bwdotcol = as.character(prefs.m[i,2])
      }
    }
    if (prefs.m[i,1] == "bwdotcex") {
      object@Prefs@Graph.prefs$bwdotcex = as.numeric(prefs.m[i,2])
    }
    if (prefs.m[i,1] == "bwrecfill") {
      if ((prefs.m[i,2]!="") && (!is.na(as.numeric(prefs.m[i,2])))) {
        object@Prefs@Graph.prefs$bwrecfill = as.numeric(prefs.m[i,2])
      } else {
        object@Prefs@Graph.prefs$bwrecfill = as.character(prefs.m[i,2])
      }
    }
    if (prefs.m[i,1] == "bwreccol") {
      if ((prefs.m[i,2]!="") && (!is.na(as.numeric(prefs.m[i,2])))) {
        object@Prefs@Graph.prefs$bwreccol = as.numeric(prefs.m[i,2])
      } else {
        object@Prefs@Graph.prefs$bwreccol = as.character(prefs.m[i,2])
      }
    }
    if (prefs.m[i,1] == "bwreclty") {
      object@Prefs@Graph.prefs$bwreclty = as.numeric(prefs.m[i,2])
    }  
    if (prefs.m[i,1] == "bwreclwd") {
      object@Prefs@Graph.prefs$bwreclwd = as.numeric(prefs.m[i,2])
    }  
    if (prefs.m[i,1] == "bwumbcol") {
      if ((prefs.m[i,2]!="") && (!is.na(as.numeric(prefs.m[i,2])))) {
        object@Prefs@Graph.prefs$bwumbcol = as.numeric(prefs.m[i,2])
      } else {
        object@Prefs@Graph.prefs$bwumbcol = as.character(prefs.m[i,2])
      }
    }  
    if (prefs.m[i,1] == "bwumblty") {
      object@Prefs@Graph.prefs$bwumblty = as.numeric(prefs.m[i,2])
    }
    if (prefs.m[i,1] == "bwumblwd") {
      object@Prefs@Graph.prefs$bwumblwd = as.numeric(prefs.m[i,2])
    }
    if (prefs.m[i,1] == "bwoutpch") {
      object@Prefs@Graph.prefs$bwoutpch = as.numeric(prefs.m[i,2])
    }   
    if (prefs.m[i,1] == "bwoutcol") {
      if ((prefs.m[i,2]!="") && (!is.na(as.numeric(prefs.m[i,2])))) {
        object@Prefs@Graph.prefs$bwoutcol = as.numeric(prefs.m[i,2])
      } else {
        object@Prefs@Graph.prefs$bwoutcol = as.character(prefs.m[i,2])
      }
    }
    if (prefs.m[i,1] == "bwoutcex") {
      object@Prefs@Graph.prefs$bwoutcex = as.numeric(prefs.m[i,2])
    }  
    
    # Histograms
    
    if (prefs.m[i,1] == "hiborder") {
      if ((prefs.m[i,2]!="") && (!is.na(as.numeric(prefs.m[i,2])))) {
        object@Prefs@Graph.prefs$hiborder = as.numeric(prefs.m[i,2])
      } else {
        object@Prefs@Graph.prefs$hiborder = as.character(prefs.m[i,2])
      }
    }
    if (prefs.m[i,1] == "hicol") {
      if ((prefs.m[i,2]!="") && (!is.na(as.numeric(prefs.m[i,2])))) {
        object@Prefs@Graph.prefs$hicol = as.numeric(prefs.m[i,2])
      } else {
        object@Prefs@Graph.prefs$hicol = as.character(prefs.m[i,2])
      }
    } 
    if (prefs.m[i,1] == "hidcol") {
      if ((prefs.m[i,2]!="") && (!is.na(as.numeric(prefs.m[i,2])))) {
        object@Prefs@Graph.prefs$hidcol = as.numeric(prefs.m[i,2])
      } else {
        object@Prefs@Graph.prefs$hidcol = as.character(prefs.m[i,2])
      }
    }
    if (prefs.m[i,1] == "hilty") {
      object@Prefs@Graph.prefs$hilty = as.numeric(prefs.m[i,2])
    }
    if (prefs.m[i,1] == "hilwd") {
      object@Prefs@Graph.prefs$hilwd = as.numeric(prefs.m[i,2])
    }  
    if (prefs.m[i,1] == "hidlty") {
      object@Prefs@Graph.prefs$hidlty = as.numeric(prefs.m[i,2])
    }
    if (prefs.m[i,1] == "hidlwd") {
      object@Prefs@Graph.prefs$hidlwd = as.numeric(prefs.m[i,2])
    }   

  } 
  
  object@Prefs@Graph.prefs$PIlimits = c(PIdolimit, PIuplimit)  
  
  options(warn = 1)
  return(object)
  
}