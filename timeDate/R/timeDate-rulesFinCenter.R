
# This R package is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This R package is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this R package; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                 DESCRIPTION:
#  rulesFinCenter            Returns DST rules for a financial center
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
rulesFinCenter <-
    function(FinCenter = "")
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Show the day light saving rules for a financial center

    # Arguments:
    #   FinCenter - a character string with the the location of the
    #       financial center named as "continent/city".

    # FUNCTION:
    if (FinCenter == "")
        FinCenter <- getRmetricsOptions("myFinCenter")

    # Check:
    if (any(FinCenter %in% c("GMT", "UTC", " ")))
        stop("There are no DST rules for GMT FinCenter!")

    # Internal Function for Conversion from Ical Tables:
    #    if (FALSE) {
    #    rulesFinCenter2 =
    #    function(FinCenter = myFinCenter) {
    #        # A function implemented by Diethelm Wuertz
    #        # Description:
    #        #   Show the day light saving rules for a financial center
    #        # Arguments:
    #        #   FinCenter - a character string with the the location of the
    #        #       financial center named as "continent/city".
    #        # Value:
    #        #   Returns a printed list of DST rules.
    #        # Example:
    #        #   > rulesFinCenter("Zurich")
    #        #               ruleChanges offSet
    #        #   1   1894-05-31 23:30:16   3600
    #        #   2   1940-11-01 23:00:00   7200
    #        #   3   1940-12-30 22:00:00   3600
    #        #   5   1941-10-04 22:00:00   3600
    #        #   6   1942-05-03 01:00:00   7200
    #        #   7   1942-10-03 22:00:00   3600
    #        #   8   1980-12-31 23:00:00   3600
    #        #   9   1981-03-29 01:00:00   7200
    #        #   ...
    #        # Note:
    #        #   Important, the "TZ" environment variable must set
    #        #   to "GMT" in your Windows Environment!
    #
    #        # Set Timezone to GMT:
    #        myTZ = Sys.getenv("TZ")
    #        Sys.setenv(TZ = "GMT")
    #        if (FinCenter == "") FinCenter = "GMT"
    #
    #        # Read the Rules:
    #        # Get IcalPath from .FirstLib
    #        file = paste0(IcalPath, FinCenter)
    #        zfile <- zip.file.extract(file, "Rdata.zip")
    #        ical = read.table(zfile, skip = 2)
    #
    #        # GMT Offsets:
    #        hm = as.integer(ical[,6])
    #        sg = sign(hm)
    #        hm = abs(hm)
    #        h = floor(hm/100)
    #        hms.off = sg * ( floor(hm/100)*3600 + (hm - 100*h)*60 + 0 )
    #        hms.off
    #
    #        # When have the rules changed?
    #        months.num = 1:12
    #        names(months.num) = c(
    #            "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    #            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    #        Y = as.integer(ical[,4])
    #        m = as.integer(months.num[as.character(ical[,3])])
    #        d = as.integer(ical[,2])
    #        CCYYMMDD = as.character(Y*10000+100*m+d)
    #        hms = unlist(strsplit(as.character(ical[,5]), ":"))
    #        hms = matrix(as.integer(hms), byrow=TRUE, ncol=3)
    #        hms = 1000000 + 10000*hms[,1] + 100*hms[,2] + hms[,3]
    #        hhmmss = substr(as.character(hms), 2, 7)
    #        ruleChangesGMT = strptime(paste(CCYYMMDD, hhmmss), "%Y%m%d %H%M%S")
    #        attr(ruleChangesGMT, "tzone") <- "GMT"
    #
    #        # Return Value:
    #        Sys.setenv(TZ = myTZ)
    #        data.frame(ruleChanges = as.character(ruleChangesGMT),
    #            offSet = hms.off) }
    #    }
    #    ## Instead:

    # Match City:
    fccity <- strsplit(FinCenter, "/")[[1]]
    City <- fccity[length(fccity)]
    fun <- match.fun(City)

    # Return Value:
    fun()
}


################################################################################

