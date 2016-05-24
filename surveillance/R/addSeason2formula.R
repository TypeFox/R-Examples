################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Conveniently add sine-cosine terms to a model formula
###
### Copyright (C) 2010 Michaela Paul, 2013-2015 Sebastian Meyer
### $Revision: 1299 $
### $Date: 2015-03-27 14:19:29 +0100 (Fre, 27. Mär 2015) $
################################################################################

## for S = 1, 'sin(2*pi * t/period) + cos(2*pi * t/period)' is added to 'f'
addSeason2formula <- function (
    f = ~1,       # formula to enhance
    S = 1,        # number of sine/cosine pairs
    period = 52,  # periodicity of the sinusoidal wave
    timevar = "t" # name of the time variable
){
    ## check arguments
    stopifnot(inherits(f, "formula"),
              is.vector(S, mode = "numeric"), S >= 0,
              isScalar(period))
    
    ## return unchanged formula if S = 0
    if (max(S) == 0)
        return(f)

    ## character representation of old formula
    ftext <- paste0(deparse(f), collapse = "")
    
    ## add sine-cosine terms
    if (length(S) == 1L) {
        for (i in seq_len(S)) {
            ftext <- paste0(ftext,
                " + sin(", 2*i, "*pi*", timevar, "/", period, ")",
                " + cos(", 2*i, "*pi*", timevar, "/", period, ")")
        }
    } else {
        ## unit-specific seasonality for hhh4() via the special fe() function
        for (i in seq_len(max(S))) {
            which <- paste0(i <= S, collapse = ",")
            ftext <- paste0(ftext,
                " + fe(sin(",2*i,"*pi*",timevar,"/",period,"), which=c(",which,"))",
                " + fe(cos(",2*i,"*pi*",timevar,"/",period,"), which=c(",which,"))")
        }
    }

    ## convert back to a formula
    as.formula(ftext, env = .GlobalEnv)
}
