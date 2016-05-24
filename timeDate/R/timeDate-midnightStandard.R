
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

# Copyrights (C)
# for this R-port:
#   1999 - Diethelm Wuertz, GPL
#   2007 - Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@phys.ethz.ch>
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:                 DESCRIPTION:
#  midnightStandard          Corrects midnight standard called by 'timeDate'
# DEPRECATED:
#  .midnightStandard
################################################################################

## # YC :midnigStandard2 returns object in POSIXct and avoid
## # wasting time in strptime

if(getRversion() < "2.15") {
# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
paste0 <- function(...) paste(..., sep = '')
}

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
midnightStandard2 <- function(charvec, format)
{
    # A function written by Diethelm Wuertz
    # and entirely rewritten by Martin Maechler
    # modifications for speed improvement by Yohan Chalabi

    # Description:
    #   Midnight Standard & conversion to isoFormat:

    # FUNCTION:

    if(all(is.na(charvec)))
        return(as.POSIXct(charvec))
    ## Motivation: strptime() {et al}  cannot deal with "24:00:00"
    ##         In that case, subtract 1 seconds convert and re-add it

    # Missing Format:
    if (missing(format)) format <- whichFormat(charvec)

    ## convert to strptime and inspect NA's returned by strptime
    ans <- as.POSIXct(strptime(charvec, format, tz = "GMT"))
    if (any(idx <- is.na(ans))) {

        # inspect problematic dates
        charvec <- charvec[idx]

        # Format:
        rng.nch <- range(nchar(charvec))
        if(rng.nch[1] != rng.nch[2])
            stop("'charvec' has non-NA entries of different number of characters")
        nch <- rng.nch[1]
        n <- length(charvec)
        s <- numeric(n)

        ## Do two common formats *fast* (for large n), and then use
        ## flexible approach:

        # ISO-8601 Midnight Standard:
        if (length(grep("%H:%M:%S", format, fixed = TRUE)) == 1) {
            ii <- grep("24:00:00", charvec, fixed=TRUE, useBytes=TRUE)
            if (length(ii) > 0) {
                s[ii] <- 1
                charvec[ii] <- gsub("24:00:00", "23:59:59", charvec[ii], fixed=TRUE)
            }
        } else if (length(grep("%H%M%S$", format)) == 1) {
            ## format *ends* in  %H%M%S, i.e. last 6 chars are time
            ch.time <- substr(charvec, nch-6+1, nch)
            if (length(ii <- grep("240000$", ch.time)) > 0) {
                s[ii] <- 1
                charvec[ii] <- paste(substr(charvec[ii], 1, nch-6),
                                     gsub("240000$", "235959", ch.time[ii]), sep = "")
            }
        } else {
            ## Very general approach, to work for any valid format:
            forms <- c("%Y", "%m", "%d",  "%H","%M","%S")
            nums  <- c("2003","01","31",  "23","59","58") # pairwise different
            fDate <- format
            for(i in seq_along(forms)) {
                ## make sure, we don't have nums[i] already :
                if(length(grep(nums[i], fDate, fixed=TRUE)))
                    fDate <- gsub(nums[i], paste(rep("x", nchar(nums[i])), collapse=""),
                                  fDate, fixed=TRUE)
                fDate <- sub(forms[i], nums[i], fDate, fixed=TRUE)
            }
            ## in the ISO case, now have  fDate == "2001-01-31 23:59:58"
            names(nums) <- forms
            ## at which character positions in charvec do I need to look for %H, ... :
            iHMS <- sapply(nums[c("%H","%M","%S")], regexpr, text=fDate, fixed=TRUE)
            if(iHMS["%H"] >= 1) {
                ## have "%H" -- otherwise, nothing to do!
                has.S <- iHMS["%S"] >= 1
                has.M <- iHMS["%M"] >= 1
                if(has.S && !has.M) stop("invalid format: has '%S' but no '%M'")
                ## 3 remaining cases:  (H,M,S), (H,M), (H)
                m. <- 1 + has.M + has.S # in {1,2,3}
                HMStab <- matrix(unlist(lapply(iHMS[seq_len(m.)],
					       function(ic)
					       substr(charvec, start=ic, stop=ic+1L))),
				 n, m.)
                twenty4 <- paste0("24", if(has.M)"00", if(has.S)"00")
                isMidN <- twenty4 == apply(HMStab, 1, paste, collapse='')
                if(any(isMidN)) {
                    ## need midnight correction
                    s[isMidN] <- 1
                    ## now *need* seconds, so we can subtract and add 1 sec :
                    if(!has.S) {
                        if(!has.M) {
                            iHMS["%M"] <- nchar(fDate) + 1
                            format <-  paste0(format,  "%M")
                            fDate  <-  paste0(fDate,   "00")
                            charvec <- paste0(charvec, "00")
                        }
                        iHMS["%S"] <- nchar(fDate) + 1
                        format <-  paste0(format,  "%S")
                        charvec <- paste0(charvec, "00")
                    }
                    substr(charvec[isMidN], iHMS["%H"], iHMS["%H"]+1) <- "23"
                    substr(charvec[isMidN], iHMS["%M"], iHMS["%M"]+1) <- "59"
                    substr(charvec[isMidN], iHMS["%S"], iHMS["%S"]+1) <- "59"
                }
            }
        }

        ## Convert "charvec" to standard ISO format:
        ## YC: added tz = "GMT" to avoid confusion when DST is active
        ans[idx] <- s + as.POSIXct(strptime(charvec, format, tz = "GMT"))
    }

    # Return Value:
    ans
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
midnightStandard <- function(charvec, format)
{
    # YC: uses now the faster midngithStandard2() function
    # but still return a character

    # Description:
    #   Midnight Standard & conversion to isoFormat:

    # FUNCTION:

    # Missing Format:
    if (missing(format)) format <- whichFormat(charvec)
    ans <- midnightStandard2(charvec, format = format)
    ans <- format(ans, format = "%Y-%m-%d %H:%M:%S")

    # Return Value:
    ans
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
.midnightStandard <- midnightStandard


################################################################################

