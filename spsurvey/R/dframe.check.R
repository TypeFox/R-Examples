dframe.check <- function(sites, design, subpop, data.cat, data.cont,
   data.risk, design.names) {

################################################################################
# Function: dframe.check
# Programmer: Tom Kincaid
# Date: September 26, 2003
# Last Revised: October 10, 2012
# Description:
#   This function checks site IDs, the sites data frame, the subpop data frame,
#      the data.cat data frame, the data.cont data frame, the data.ar data
#      frame, and the data.rr data frame to assure valid contents.  If they do
#      not exist, then the sites data frame and the subpop data frame are
#      created.
#   Input:
#      design = the design data frame.
#      sites = the sites data frame.
#      subpop = the subpop data frame.
#      data.cat = the data.cat data frame of categorical response variables.
#      data.cont = the data.cont data frame of continuous response variables.
#      data.risk = the data.ar or data.rr data frame of categorical response and
#                  stressor variables.
#      design.names = names for the design data frame.
#   Output:
#      A list consisting of the sites data frame, design data frame, subpop data
#      frame, data.cat data frame, and data.cont data frame.
#   Other Functions Required:
#      vecprint - takes an input vector and outputs a character string with
#         line breaks inserted
################################################################################

# Check the sites data frame for contents

   if(is.null(sites)) {
      sites <- data.frame(siteID=design$siteID,
                          use.sites=rep(TRUE, nrow(design)))
      temp <- is.na(sites[,1])
      if(any(temp)) {
         temp.str <- vecprint(seq(nrow(sites))[temp])
         stop(paste("\nThe following rows in the sites data frame contain missing site ID values:\n", temp.str, sep=""))
      }
      temp <- sapply(split(sites[,1], sites[,1]), length)
      if(any(temp > 1)) {
         temp.str <- vecprint(names(temp)[temp > 1])
         stop(paste("The following site ID values in the sites data frame occur more than \nonce:\n", temp.str, sep=""))
      }
      siteID <- sites$siteID
   } else {
      if(!is.data.frame(sites))
         stop("\nThe sites argument must be a data frame.")
      if(ncol(sites) != 2)
         stop("\nThe sites argument must contain exactly two variables.")
      temp <- is.na(sites[,1])
      if(any(temp)) {
         temp.str <- vecprint(seq(nrow(sites))[temp])
         stop(paste("\nThe following rows in the sites data frame contain missing site ID values:\n", temp.str, sep=""))
      }
      temp <- is.na(sites[,2])
      if(any(temp)) {
         temp.str <- vecprint(seq(nrow(sites))[temp])
         stop(paste("\nThe following rows in the sites data frame contain missing logical variable values:\n", temp.str, sep=""))
      }
      if(!is.logical(sites[,2]))
         stop("\nThe second variable in the sites data frame is not a logical variable.")
      siteID <- uniqueID(sites[,1])[sites[,2]]
      sites <- sites[sites[,2],]
      temp <- sapply(split(sites[,1], sites[,1]), length)
      if(any(temp > 1)) {
         temp.str <- vecprint(names(temp)[temp > 1])
         stop(paste("The following site ID values in the sites data frame occur more than \nonce:\n", temp.str, sep=""))
      }
   }
   names(sites)[1] <- design.names[1]

# Check the design data frame for contents

   temp <- is.na(design$siteID)
   if(any(temp)) {
      temp.str <- vecprint(seq(nrow(design))[temp])
      stop(paste("\nThe following rows in the design data frame contain missing site ID values:\n", temp.str, sep=""))
   }
   temp <- match(siteID, uniqueID(design$siteID), nomatch=0)
   if(any(temp == 0)) {
      temp.str <- vecprint(unique(siteID[temp == 0]))
      stop(paste("\nThe following site ID values in the sites data frame do not occur among the \nsite ID values in the design data frame:\n", temp.str, sep=""))
   }
   design <- design[temp,]

# Check the subpop data frame for contents

   if(is.null(subpop)) {
      subpop <- data.frame(siteID=siteID,
         all.sites=rep("All Sites", nrow(sites)))
   } else {
      if(!is.data.frame(subpop))
         stop("\nThe subpop argument must be a data frame.")
      if(ncol(subpop) < 2)
         stop("\nThe subpop argument must contain at least two variables.")
      temp <- is.na(subpop[,1])
      if(any(temp)) {
         temp.str <- vecprint(seq(nrow(subpop))[temp])
         stop(paste("\nThe following rows in the subpop data frame contain missing site ID values:\n", temp.str, sep=""))
      }
      temp <- match(siteID, uniqueID(subpop[,1]), nomatch=0)
      if(any(temp == 0)) {
         temp.str <- vecprint(unique(siteID[temp == 0]))
         stop(paste("\nThe following site ID values in the sites data frame do not occur among the \nsite ID values in the subpop data frame:\n", temp.str, sep=""))
      }
      subpop <- subpop[temp,]
   }
   names(subpop)[1] <- design.names[1]

# Check the data.cat data frame for contents

   if(!is.null(data.cat)) {
      if(!is.data.frame(data.cat))
         stop("\nThe data.cat argument must be a data frame.")
      temp <- is.na(data.cat[,1])
      if(any(temp)) {
         temp.str <- vecprint(seq(nrow(data.cat))[temp])
         stop(paste("\nThe following rows in the data.cat data frame contain missing site ID values:\n", temp.str, sep=""))
      }
      temp <- match(siteID, uniqueID(data.cat[,1]), nomatch=0)
      if(any(temp == 0)) {
         temp.str <- vecprint(unique(siteID[temp == 0]))
         stop(paste("\nThe following site ID values in the sites data frame do not occur among the \nsite ID values in the data.cat data frame:\n", temp.str, sep=""))
      }
      data.cat <- data.cat[temp,]
      names(data.cat)[1] <- design.names[1]
   }

# Check the data.cont data frame for contents

   if(!is.null(data.cont)) {
      if(!is.data.frame(data.cont))
         stop("\nThe data.cont argument must be a data frame.")
      temp <- is.na(data.cont[,1])
      if(any(temp)) {
         temp.str <- vecprint(seq(nrow(data.cont))[temp])
         stop(paste("\nThe following rows in the data.cont data frame contain missing site ID values:\n", temp.str, sep=""))
      }
      temp <- match(siteID, uniqueID(data.cont[,1]), nomatch=0)
      if(any(temp == 0)) {
         temp.str <- vecprint(unique(siteID[temp == 0]))
         stop(paste("\nThe following site ID values in the sites data frame do not occur among the \nsite ID values in the data.cont data frame:\n", temp.str, sep=""))
      }
      data.cont <- data.cont[temp,]
      names(data.cont)[1] <- design.names[1]
   }

# Check the data.risk data frame for contents

   if(!is.null(data.risk)) {
      if(!is.data.frame(data.risk))
         stop("\nThe data.risk argument must be a data frame.")
      temp <- is.na(data.risk[,1])
      if(any(temp)) {
         temp.str <- vecprint(seq(nrow(data.risk))[temp])
         stop(paste("\nThe following rows in the data.risk data frame contain missing site ID values:\n", temp.str, sep=""))
      }
      temp <- match(siteID, uniqueID(data.risk[,1]), nomatch=0)
      if(any(temp == 0)) {
         temp.str <- vecprint(unique(siteID[temp == 0]))
         stop(paste("\nThe following site ID values in the sites data frame do not occur among the \nsite ID values in the data.risk data frame:\n", temp.str, sep=""))
      }
      data.risk <- data.risk[temp,]
      names(data.risk)[1] <- design.names[1]
   }

# Return the list

   list(sites=sites, design=design, subpop=subpop, data.cat=data.cat,
        data.cont=data.cont, data.risk=data.risk)
}
