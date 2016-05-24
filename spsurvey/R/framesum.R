framesum <- function(att.frame, design, type.frame="finite", stratum=NULL,
   mdcaty=NULL, auxvar=NULL, units.in="Number", scale=1, units.out="Number") {

################################################################################
# Function: framesum
# Purpose: Summarize frame size for a survey design
# Programmers: Tony Olsen, Tom Kincaid
# Date: April 27, 2005
# Last Revised: February 27, 2006
# Description:
#   This function summarizes the frame for a survey design.  When type.frame
#   equals "finite", summary is a count of number of units in att.frame for
#   cross-tabulation of stratum, mdcaty, and auxvar.  When type.frame equals
#   "linear" or "area", summary is the sum of length or area for units for
#   cross- tabulation of stratum, mdcaty, and auxvar.  Note that length and area
#   are taken from length_mdm and area_mdm, which are calculated by the function
#   read.dbf and added to att.frame.  If argument mdcaty or argument stratum
#   equals NULL or if both arguments equal NULL, then the cross-tabulation is
#   performed without use of the design variable(s).  
# Arguments
#   att.frame = a data frame composed of attributes associated with elements in
#     the frame, which must contain the columns used for stratum and mdcaty (if
#     required by the survey design).
#   design = named list of stratum design specifications which are also lists.
#     Stratum names must be subset of values in stratum argument.  Each stratum
#     list has four list components:
#       panel = named vector of sample sizes for each panel in stratum
#       seltype = the type of random selection, which must be one of following:
#         "Equal" - equal probability selection, "Unequal" - unequal probability
#         selection by the categories specified in caty.n and mdcaty, or
#         "Continuous" - unequal probability selection proportional to auxiliary
#         variable mdcaty
#       caty.n = if seltype equals "Unequal", a named vector of sample sizes for
#         each category specified by mdcaty, where sum of the sample sizes must
#         equal sum of the panel sample sizes, and names must be a subset of
#         values in mdcaty
#       over = number of replacement sites ("oversample" sites) for the entire
#         design, which is set equal to 0 if none are required.
#     Example design: 
#       design=list(Stratum1=list(panel=c(PanelOne=50), seltype="Equal",
#         over=10), Stratum2=list(panel=c(PanelOne=50, PanelTwo=50),
#         seltype="Unequal", caty.n=c(CatyOne=25, CatyTwo=25, CatyThree=25,
#         CatyFour=25), over=75))
#   type.frame = the type of frame, which must be one of following: "finite",
#     "linear", or "area".  The default is "finite".
#   stratum = name of the column from att.frame that identifies stratum
#     membership for each element in the frame.  If stratum equals NULL, the
#     design is unstratified.  The default is NULL.
#   mdcaty = name of the column from att.frame that identifies the unequal
#     probability category for each element in the frame.  The default is NULL.
#   auxvar = a vector containing the names of columns from att.frame that
#     identify auxiliary variables to be used to summarize frame size.  The
#     default is NULL. 
#   units.in  = a character string giving the name of units used to measure size
#     in the frame.  The default is "Number".
#   scale = the scale factor used to change units.in to units.out.  For example,
#     use 1000 to change "Meters" to "Kilometers".  The default is 1.
#   units.out = a character string giving the name of units used to measure size
#     in the results.  The default is "Number".
# Results:
#   A list containing two components named DesignSize and AuxVarSize.
#   DesignSize summarizes the frame for survey design variables, and AuxVarSize
#   summarizes the frame for auxiliary variables.  DesignSize is either a table
#   (for type.frame equals "finite") or an array (for type.frame equals "linear"
#   or "area") that contains the cross-tabulation of frame extent for design
#   variables multidensity category (mdcaty) and stratum, where extent of the
#   frame is the number of sites for type.frame equals "finite", the sum of
#   resource length for type.frame equals "linear", or the sum of resource area
#   for type.frame equals "area".  AuxVarSize is a list containing a component
#   for each auxiliary variable, where each component of the list is one of the
#   following: (1) if the type of random selection does not equal "Continuous"
#   for any stratum, each component is either a table (for type.frame equals
#   "finite") or an array (for type.frame equals "linear" or "area") that
#   contains the cross-tabulation of frame extent for mdcaty, stratum, and the
#   auxiliary variable or (2) if type of random selection equals "Continuous"
#   for all strata, each component is either a table (finite frame) or an array
#   (linear or area frame) containing the cross-tabulation of frame extent for
#   stratum and the auxiliary variable.  In addition the output list plus
#   labeling information is printed to the console.
# Other Functions Required:
#   vecprint - takes an input vector and outputs a character string with
#     line breaks inserted
# Examples:
#   attframe <- read.dbf("shapefile")
#   design <- list(Stratum1=list(panel=c(PanelOne=50), seltype="Equal",
#      over=10), Stratum2=list(panel=c(PanelOne=50, PanelTwo=50),
#      seltype="Unequal", caty.n=c(CatyOne=25, CatyTwo=25, CatyThree=25,
#      CatyFour=25), over=75))
#   framesum(att.frame=attframe, design=design, type.frame="area",
#      stratum="stratum", mdcaty="mdcaty", auxvar=c("ecoregion",
#      "state"), units.in="Meters", scale=1000, units.out="Kilometers")
################################################################################

# Determine whether stratum and mdcaty are present

   stratum.ind <- ifelse(is.null(stratum), FALSE, TRUE)
   mdcaty.ind <- ifelse(is.null(mdcaty), FALSE, TRUE)

# Determine whether the type of random selection is "Continuous" for any stratum

   seltype.ind <- FALSE
   for(i in 1:length(design)) {
      if(design[[i]]$seltype == "Continuous")
         seltype.ind <- TRUE
   }

# Summarize a finite frame

   if(type.frame == "finite") {

# Summarize design variables

      if(!seltype.ind && mdcaty.ind) {
         if(stratum.ind) {
            DesignSize <- addmargins(table(att.frame[,mdcaty],
               att.frame[,stratum], dnn=c(mdcaty, stratum)))
            cat(paste("Frame Summary: Number of Sites Classified by ", mdcaty, " (Multidensity Category) \nand ", stratum, " (Stratum)\n\n", sep=""))
            print(DesignSize)
         } else {
            DesignSize <- addmargins(table(att.frame[,mdcaty], dnn=mdcaty))
            cat(paste("Frame Summary: Number of Sites Classified by ", mdcaty, " (Multidensity Category)\n\n", sep=""))
            print(DesignSize)
         }
      } else {
         if(stratum.ind) {
            DesignSize <- addmargins(table(att.frame[,stratum], dnn=stratum))
            cat(paste("Frame Summary: Number of Sites Classified by ", stratum, " (Stratum)\n\n", sep=""))
            print(DesignSize)
         } else {
            DesignSize <- nrow(att.frame)
            cat(paste("Frame Summary: Number of Sites in the Frame is ", DesignSize, "\n", sep=""))
         }
      }

# Summarize auxiliary variables

      if(is.null(auxvar)) {
         AuxVarSize <- NULL

      } else {
         temp <- match(auxvar, names(att.frame), nomatch=0)
         if(any(temp == 0)) {
            temp.str <- vecprint(auxvar[temp == 0])
            stop(paste("\nThe following values in the vector of auxiliary variable names do not occur among the columns in the att.frame data frame:\n", temp.str, sep=""))
         }

         AuxVarSize <- list()
         for(i in auxvar) {
            if(!is.factor(att.frame[,i]))
               att.frame[,i] <- as.factor(att.frame[,i])
            if(!seltype.ind && mdcaty.ind) {
               if(stratum.ind) {
                  AuxVarSize[[i]] <- addmargins(table(att.frame[,mdcaty],
                     att.frame[,stratum], att.frame[,i], dnn=c(mdcaty,
                     stratum, i)))
                  cat(paste("\n\nFrame Summary: Number of Sites Classified by ", mdcaty, " (Multidensity Category), \n", stratum, " (Stratum), and ", i, " (Auxiliary Variable)\n\n", sep=""))
                  print(AuxVarSize[[i]])
               } else {
                  AuxVarSize[[i]] <- addmargins(table(att.frame[,mdcaty],
                     att.frame[,i], dnn=c(mdcaty, i)))
                  cat(paste("\n\nFrame Summary: Number of Sites Classified by ", mdcaty, " (Multidensity Category), \nand ", i, " (Auxiliary Variable)\n\n", sep=""))
                  print(AuxVarSize[[i]])
               }
            } else {
               if(stratum.ind) {
                  AuxVarSize[[i]] <- addmargins(table(att.frame[,stratum],
                     att.frame[,i], dnn=c(stratum, i)))
                  cat(paste("\n\nFrame Summary: Number of Sites Classified by ", stratum, " (Stratum), \nand ", i, " (Auxiliary Variable)\n\n", sep=""))
                  print(AuxVarSize[[i]])
               } else {
                  AuxVarSize[[i]] <- addmargins(table(att.frame[,i], dnn=i))
                  cat(paste("\n\nFrame Summary: Number of Sites Classified by ", i, " (Auxiliary Variable)\n\n", sep=""))
                  print(AuxVarSize[[i]])
               }
            }
         }
      }

# Summarize a linear frame

   } else if(type.frame == "linear") {

# Summarize design variables

      if(!seltype.ind && mdcaty.ind) {
         if(stratum.ind) {
            allsize <- sum(att.frame$length_mdm)
            mdsize <- tapply(att.frame$length_mdm, att.frame[,mdcaty], sum)
            stsize <- tapply(att.frame$length_mdm, att.frame[,stratum], sum)
            mdstsize <- tapply(att.frame$length_mdm, list(att.frame[,mdcaty],
               att.frame[,stratum]), sum)
            DesignSize <- rbind(cbind(mdstsize, Sum=mdsize), Sum=c(stsize,
               allsize))/scale
            names(dimnames(DesignSize)) <- list(mdcaty, stratum)
            cat(paste("Frame Summary: Resource Length Classified by ", mdcaty, " (Multidensity Category) \nand ", stratum, " (Stratum)\n\n", sep=""))
            print(DesignSize)
         } else {
            allsize <- sum(att.frame$length_mdm)
            mdsize <- tapply(att.frame$length_mdm, att.frame[,mdcaty], sum)
            DesignSize <- as.array(c(mdsize, Sum=allsize)/scale)
            names(dimnames(DesignSize)) <- mdcaty
            cat(paste("Frame Summary: Resource Length Classified by ", mdcaty, " (Multidensity Category)\n\n", sep=""))
            print(DesignSize)
         }
      } else {
         if(stratum.ind) {
            allsize <- sum(att.frame$length_mdm)
            stsize <- tapply(att.frame$length_mdm, att.frame[,stratum], sum)
            DesignSize <- as.array(c(stsize, Sum=allsize)/scale)
            names(dimnames(DesignSize)) <- stratum
            cat(paste("Frame Summary: Resource Length Classified by ", stratum, " (Stratum)\n\n", sep=""))
            print(DesignSize)
         } else {
            DesignSize <- sum(att.frame$length_mdm)/scale
            cat(paste("Frame Summary: Total Length of the Resource is " , DesignSize, "\n", sep=""))
         }
      }

# Summarize auxiliary variables

      if(is.null(auxvar)) {
         AuxVarSize <- NULL
      } else {
         temp <- match(auxvar, names(att.frame), nomatch=0)
         if(any(temp == 0)) {
            temp.str <- vecprint(auxvar[temp == 0])
            stop(paste("\nThe following values in the vector of auxiliary variable names do not occur among the columns in the att.frame data frame:\n", temp.str, sep=""))
         }

         AuxVarSize <- list()
         for(i in auxvar) {
            if(!is.factor(att.frame[,i]))
               att.frame[,i] <- as.factor(att.frame[,i])
            if(!seltype.ind && mdcaty.ind) {
               if(stratum.ind) {
                  lev <- levels(att.frame[,i])
                  temp.array <- array(dim=c(length(levels(att.frame[,mdcaty]))
                     + 1, length(levels(att.frame[,stratum])) + 1,
                     length(lev)))
                  dimnames(temp.array) <- list(c(levels(att.frame[,mdcaty]),
                     "sum"), c(levels(att.frame[,stratum]), "sum"), lev)
                  names(dimnames(temp.array)) <- list(mdcaty, stratum, i)
                  k <- 1
                  for(j in lev) {
                     temp.frame <- att.frame[att.frame[,i] == j,]
                     allsize <- sum(temp.frame$length_mdm)
                     mdsize <- tapply(temp.frame$length_mdm,
                        temp.frame[,mdcaty], sum)
                     stsize <- tapply(temp.frame$length_mdm,
                        temp.frame[,stratum], sum)
                     mdstsize <- tapply(temp.frame$length_mdm,
                        list(temp.frame[,mdcaty], temp.frame[,stratum]), sum)
                     temp.array[,,k] <- rbind(cbind(mdstsize, Sum=mdsize),
                        Sum=c(stsize, allsize))/scale
                     k <- k+1
                  }
                  AuxVarSize[[i]] <- temp.array
                  cat(paste("\n\nFrame Summary: Resource Length Classified by ", mdcaty, " (Multidensity Category), \n", stratum, " (Stratum), and ", i, " (Auxiliary Variable)\n\n", sep=""))
                  print(AuxVarSize[[i]])
               } else {
                  allsize <- sum(att.frame$length_mdm)
                  mdsize <- tapply(att.frame$length_mdm, att.frame[,mdcaty],
                     sum)
                  axsize <- tapply(att.frame$length_mdm, att.frame[,i], sum)
                  mdaxsize <- tapply(att.frame$length_mdm,
                     list(att.frame[,mdcaty], att.frame[,i]), sum)
                  temp <- rbind(cbind(mdaxsize, Sum=mdsize), Sum=c(axsize,
                     allsize))/scale
                  names(dimnames(temp)) <- list(mdcaty, i)
                  AuxVarSize[[i]] <- temp
                  cat(paste("\n\nFrame Summary: Resource Length Classified by ", mdcaty, " (Multidensity Category), \nand ", i, " (Auxiliary Variable)\n\n", sep=""))
                  print(AuxVarSize[[i]])
               }
            } else {
               if(stratum.ind) {
                  allsize <- sum(att.frame$length_mdm)
                  stsize <- tapply(att.frame$length_mdm, att.frame[,stratum],
                     sum)
                  axsize <- tapply(att.frame$length_mdm, att.frame[,i], sum)
                  staxsize <- tapply(att.frame$length_mdm,
                     list(att.frame[,stratum], att.frame[,i]), sum)
                  temp <- rbind(cbind(staxsize, Sum=stsize), Sum=c(axsize,
                     allsize))/scale
                  names(dimnames(temp)) <- list(stratum, i)
                  AuxVarSize[[i]] <- temp
                  cat(paste("\n\nFrame Summary: Resource Length Classified by ", stratum, " (Stratum), \nand ", i, " (Auxiliary Variable)\n\n", sep=""))
                  print(AuxVarSize[[i]])
               } else {
                  allsize <- sum(att.frame$length_mdm)
                  axsize <- tapply(att.frame$length_mdm, att.frame[,i], sum)
                  temp <- as.array(c(axsize, Sum=allsize)/scale)
                  names(dimnames(temp)) <- i
                  AuxVarSize[[i]] <- temp
                  cat(paste("\n\nFrame Summary: Resource Length Classified by ", i, " (Auxiliary Variable)\n\n", sep=""))
                  print(AuxVarSize[[i]])
               }
            }
         }
      }

# Summarize an area frame

   } else if(type.frame == "area") {
   
# Summarize design variables

      if(!seltype.ind && mdcaty.ind) {
         if(stratum.ind) {
            allsize <- sum(att.frame$area_mdm)
            mdsize <- tapply(att.frame$area_mdm, att.frame[,mdcaty], sum)
            stsize <- tapply(att.frame$area_mdm, att.frame[,stratum], sum)
            mdstsize <- tapply(att.frame$area_mdm, list(att.frame[,mdcaty],
               att.frame[,stratum]), sum)
            DesignSize <- rbind(cbind(mdstsize, Sum=mdsize), Sum=c(stsize,
               allsize))/scale
            names(dimnames(DesignSize)) <- list(mdcaty, stratum)
            cat(paste("Frame Summary: Resource Area Classified by ", mdcaty, " (Multidensity Category) \nand ", stratum, " (Stratum)\n\n", sep=""))
            print(DesignSize)
         } else {
            allsize <- sum(att.frame$area_mdm)
            mdsize <- tapply(att.frame$area_mdm, att.frame[,mdcaty], sum)
            DesignSize <- as.array(c(mdsize, Sum=allsize)/scale)
            names(dimnames(DesignSize)) <- mdcaty
            cat(paste("Frame Summary: Resource Area Classified by ", mdcaty, " (Multidensity Category)\n\n", sep=""))
            print(DesignSize)
         }
      } else {
         if(stratum.ind) {
            allsize <- sum(att.frame$area_mdm)
            stsize <- tapply(att.frame$area_mdm, att.frame[,stratum], sum)
            DesignSize <- as.array(c(stsize, Sum=allsize)/scale)
            names(dimnames(DesignSize)) <- stratum
            cat(paste("Frame Summary: Resource Area Classified by ", stratum, " (Stratum)\n\n", sep=""))
            print(DesignSize)
         } else {
            DesignSize <- sum(att.frame$area_mdm)/scale
            cat(paste("Frame Summary: Total Area of the Resource is " , DesignSize, "\n", sep=""))
         }
      }

# Summarize auxiliary variables

      if(is.null(auxvar)) {
         AuxVarSize <- NULL
      } else {
         temp <- match(auxvar, names(att.frame), nomatch=0)
         if(any(temp == 0)) {
            temp.str <- vecprint(auxvar[temp == 0])
            stop(paste("\nThe following values in the vector of auxiliary variable names do not occur among the columns in the att.frame data frame:\n", temp.str, sep=""))
         }

         AuxVarSize <- list()
         for(i in auxvar) {
            if(!is.factor(att.frame[,i]))
               att.frame[,i] <- as.factor(att.frame[,i])
            if(!seltype.ind && mdcaty.ind) {
               if(stratum.ind) {
                  lev <- levels(att.frame[,i])
                  temp.array <- array(dim=c(length(levels(att.frame[,mdcaty]))
                     + 1, length(levels(att.frame[,stratum])) + 1,
                     length(lev)))
                  dimnames(temp.array) <- list(c(levels(att.frame[,mdcaty]),
                     "sum"), c(levels(att.frame[,stratum]), "sum"), lev)
                  names(dimnames(temp.array)) <- list(mdcaty, stratum, i)
                  k <- 1
                  for(j in lev) {
                     temp.frame <- att.frame[att.frame[,i] == j,]
                     allsize <- sum(temp.frame$area_mdm)
                     mdsize <- tapply(temp.frame$area_mdm,
                        temp.frame[,mdcaty], sum)
                     stsize <- tapply(temp.frame$area_mdm,
                        temp.frame[,stratum], sum)
                     mdstsize <- tapply(temp.frame$area_mdm,
                        list(temp.frame[,mdcaty], temp.frame[,stratum]), sum)
                     temp.array[,,k] <- rbind(cbind(mdstsize, Sum=mdsize),
                        Sum=c(stsize, allsize))/scale
                     k <- k+1
                  }
                  AuxVarSize[[i]] <- temp.array
                  cat(paste("\n\nFrame Summary: Resource Area Classified by ", mdcaty, " (Multidensity Category), \n", stratum, " (Stratum), and ", i, " (Auxiliary Variable)\n\n", sep=""))
                  print(AuxVarSize[[i]])
               } else {
                  allsize <- sum(att.frame$area_mdm)
                  mdsize <- tapply(att.frame$area_mdm, att.frame[,mdcaty],
                     sum)
                  axsize <- tapply(att.frame$area_mdm, att.frame[,i], sum)
                  mdaxsize <- tapply(att.frame$area_mdm,
                     list(att.frame[,mdcaty], att.frame[,i]), sum)
                  temp <- rbind(cbind(mdaxsize, Sum=mdsize), Sum=c(axsize,
                     allsize))/scale
                  names(dimnames(temp)) <- list(mdcaty, i)
                  AuxVarSize[[i]] <- temp
                  cat(paste("\n\nFrame Summary: Resource Area Classified by ", mdcaty, " (Multidensity Category), \nand ", i, " (Auxiliary Variable)\n\n", sep=""))
                  print(AuxVarSize[[i]])
               }
            } else {
               if(stratum.ind) {
                  allsize <- sum(att.frame$area_mdm)
                  stsize <- tapply(att.frame$area_mdm, att.frame[,stratum],
                     sum)
                  axsize <- tapply(att.frame$area_mdm, att.frame[,i], sum)
                  staxsize <- tapply(att.frame$area_mdm,
                     list(att.frame[,stratum], att.frame[,i]), sum)
                  temp <- rbind(cbind(staxsize, Sum=stsize), Sum=c(axsize,
                     allsize))/scale
                  names(dimnames(temp)) <- list(stratum, i)
                  AuxVarSize[[i]] <- temp
                  cat(paste("\n\nFrame Summary: Resource Area Classified by ", stratum, " (Stratum), \nand ", i, " (Auxiliary Variable)\n\n", sep=""))
                  print(AuxVarSize[[i]])
               } else {
                  allsize <- sum(att.frame$area_mdm)
                  axsize <- tapply(att.frame$area_mdm, att.frame[,i], sum)
                  temp <- as.array(c(axsize, Sum=allsize)/scale)
                  names(dimnames(temp)) <- i
                  AuxVarSize[[i]] <- temp
                  cat(paste("\n\nFrame Summary: Resource Area Classified by ", i, " (Auxiliary Variable)\n\n", sep=""))
                  print(AuxVarSize[[i]])
               }
            }
         }
      }
   }

# Create the output list

   rslt <- list(DesignSize=DesignSize, AuxVarSize=AuxVarSize)
   attr(rslt, "units") <- units.out

# Return the list

   invisible(rslt)
}
