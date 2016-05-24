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

"xpose.write" <-
  function(object,
           file = "xpose.ini"
           ) {
           
  # build matrix (66 options)
  prefs <- matrix(1:150, ncol=2)
  
  # General
  prefs[1,]  <- c("Miss",          object@Prefs@Miss)   
  prefs[2,]  <- c("Cat.levels",    object@Prefs@Cat.levels) 
  prefs[3,]  <- c("DV.Cat.levels", object@Prefs@DV.Cat.levels) 
  
  # Plotting    
  prefs[4,]  <- c("type",          object@Prefs@Graph.prefs$type)
  prefs[5,]  <- c("cex",           object@Prefs@Graph.prefs$cex)
  prefs[6,]  <- c("lty",           object@Prefs@Graph.prefs$lty)
  prefs[7,]  <- c("lwd",           object@Prefs@Graph.prefs$lwd)
  prefs[8,]  <- c("col",           object@Prefs@Graph.prefs$col)   
  prefs[9,]  <- c("pch",           object@Prefs@Graph.prefs$pch)  
  prefs[10,] <- c("grid",          object@Prefs@Graph.prefs$grid)  
  prefs[11,] <- c("aspect",        object@Prefs@Graph.prefs$aspect) 
  
  # Conditioning
  prefs[12,] <- c("byordfun",      object@Prefs@Graph.prefs$byordfun)
  prefs[13,] <- c("shingnum",      object@Prefs@Graph.prefs$shingnum)
  prefs[14,] <- c("shingol",       object@Prefs@Graph.prefs$shingol)
  
  # abline
  prefs[15,]  <- c("ablcol",       object@Prefs@Graph.prefs$ablcol)
  prefs[16,]  <- c("abllty",       object@Prefs@Graph.prefs$abllty)
  prefs[17,]  <- c("abllwd",       object@Prefs@Graph.prefs$abllwd)   
  
  # lmline
  prefs[18,]  <- c("lmcol",        object@Prefs@Graph.prefs$lmcol)
  prefs[19,]  <- c("lmlty",        object@Prefs@Graph.prefs$lmlty)
  prefs[20,]  <- c("lmlwd",        object@Prefs@Graph.prefs$lmlwd)  
  
  # smooth
  prefs[21,]  <- c("smooth",       object@Prefs@Graph.prefs$smooth)
  prefs[22,]  <- c("smcol",        object@Prefs@Graph.prefs$sucol)
  prefs[23,]  <- c("smlty",        object@Prefs@Graph.prefs$sulty)
  prefs[24,]  <- c("smlwd",        object@Prefs@Graph.prefs$sulwd)   
  prefs[25,]  <- c("smspan",       object@Prefs@Graph.prefs$suspan)
  prefs[26,]  <- c("smdegr",       object@Prefs@Graph.prefs$sudegr)  
  
  # suline
  prefs[27,]  <- c("sucol",        object@Prefs@Graph.prefs$sucol)
  prefs[28,]  <- c("sulty",        object@Prefs@Graph.prefs$sulty)
  prefs[29,]  <- c("sulwd",        object@Prefs@Graph.prefs$sulwd)   
  prefs[30,]  <- c("suspan",       object@Prefs@Graph.prefs$suspan)
  prefs[31,]  <- c("sudegr",       object@Prefs@Graph.prefs$sudegr)  
  
  # Labelling
  prefs[32,]  <- c("ids",          object@Prefs@Graph.prefs$ids)
  prefs[33,]  <- c("idsext",       object@Prefs@Graph.prefs$idsext)
  prefs[34,]  <- c("idscex",       object@Prefs@Graph.prefs$idscex)   
  prefs[35,]  <- c("idsdir",       object@Prefs@Graph.prefs$idsdir)
  
  # Dilution
  prefs[36,]  <- c("dilfrac",      object@Prefs@Graph.prefs$dilfrac)   
  prefs[37,]  <- c("dilci",        object@Prefs@Graph.prefs$dilci)  
  
  # Prediction intervals
  prefs[38,]  <- c("PIuplty",      object@Prefs@Graph.prefs$PIuplty)
  prefs[39,]  <- c("PIdolty",      object@Prefs@Graph.prefs$PIdolty)
  prefs[40,]  <- c("PImelty",      object@Prefs@Graph.prefs$PImelty)
  prefs[41,]  <- c("PIuptyp",      object@Prefs@Graph.prefs$PIuptyp)   
  prefs[42,]  <- c("PIdotyp",      object@Prefs@Graph.prefs$PIdotyp)
  prefs[43,]  <- c("PImetyp",      object@Prefs@Graph.prefs$PImetyp)  
  prefs[44,]  <- c("PIupcol",      object@Prefs@Graph.prefs$PIupcol)
  prefs[45,]  <- c("PIdocol",      object@Prefs@Graph.prefs$PIdocol)
  prefs[46,]  <- c("PImecol",      object@Prefs@Graph.prefs$PImecol)
  prefs[47,]  <- c("PIuplwd",      object@Prefs@Graph.prefs$PIuplwd)
  prefs[48,]  <- c("PIdolwd",      object@Prefs@Graph.prefs$PIdolwd)
  prefs[49,]  <- c("PImelwd",      object@Prefs@Graph.prefs$PImelwd)  
  prefs[50,]  <- c("PIuplimit",    object@Prefs@Graph.prefs$PIlimits[2])
  prefs[51,]  <- c("PIdolimit",    object@Prefs@Graph.prefs$PIlimits[1]) 
  prefs[52,]  <- c("PIarcol",      object@Prefs@Graph.prefs$PImelwd)   
  
  # B&W plots
  prefs[53,]  <- c("bwhoriz",      object@Prefs@Graph.prefs$bwhoriz)
  prefs[54,]  <- c("bwratio",      object@Prefs@Graph.prefs$bwratio)
  prefs[55,]  <- c("bwvarwid",     object@Prefs@Graph.prefs$bwvarwid)
  prefs[56,]  <- c("bwdotpch",     object@Prefs@Graph.prefs$bwdotpch)   
  prefs[57,]  <- c("bwdotcol",     object@Prefs@Graph.prefs$bwdotcol)
  prefs[58,]  <- c("bwdotcex",     object@Prefs@Graph.prefs$bwdotcex)  
  prefs[59,]  <- c("bwrecfill",    object@Prefs@Graph.prefs$bwrecfill)   
  prefs[60,]  <- c("bwreccol",     object@Prefs@Graph.prefs$bwreccol)
  prefs[61,]  <- c("bwreclty",     object@Prefs@Graph.prefs$bwreclty) 
  prefs[62,]  <- c("bwreclwd",     object@Prefs@Graph.prefs$bwreclwd)
  prefs[63,]  <- c("bwumbcol",     object@Prefs@Graph.prefs$bwumbcol)
  prefs[64,]  <- c("bwumblty",     object@Prefs@Graph.prefs$bwumblty) 
  prefs[65,]  <- c("bwumblwd",     object@Prefs@Graph.prefs$bwumblwd)
  prefs[66,]  <- c("bwoutpch",     object@Prefs@Graph.prefs$bwoutpch)   
  prefs[67,]  <- c("bwoutcol",     object@Prefs@Graph.prefs$bwoutcol)
  prefs[68,]  <- c("bwoutcex",     object@Prefs@Graph.prefs$bwoutcex)   
  
  # Histograms
  prefs[69,]  <- c("hiborder",     object@Prefs@Graph.prefs$hiborder)
  prefs[70,]  <- c("hicol",        object@Prefs@Graph.prefs$hicol)
  prefs[71,]  <- c("hilty",        object@Prefs@Graph.prefs$hilty)
  prefs[72,]  <- c("hilwd",        object@Prefs@Graph.prefs$hilwd)   
  prefs[73,]  <- c("hidcol",       object@Prefs@Graph.prefs$hidcol)
  prefs[74,]  <- c("hidlty",       object@Prefs@Graph.prefs$hidlty)  
  prefs[75,]  <- c("hidlwd",       object@Prefs@Graph.prefs$hidlwd)   
  
  for (i in 1:nrow(prefs)) {
    if (prefs[i,1] == prefs[i,2]) {
      prefs[i,2] = "NULL"
    }
  }

  # save matrix
  write.table(prefs, file = file, sep = "\t", col.names = FALSE, row.names=FALSE,
              quote=FALSE)
}