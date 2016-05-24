###############################################################################
## Fit the core NMMAPS model with PM10 and mortality
## Copyright (C) 2004, Roger D. Peng <rpeng@jhsph.edu> 
##     
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
###############################################################################

## `fitCity' should be used in conjunction with the `tsModelSpec'
## package and the `NMMAPSdata' package.  It fits the basic/core
## NMMAPS model with a smooth function of time, smooth functions of
## time for age categories 2, 3, and smooth functions of temp and dew
## point.  Data should be preprocessed with `basicNMMAPS' (or
## something similar).

## The `formula' argument is used to specify the response and the
## pollutant portion of the model.


fitCity <- function(data, formula, df.Time = 7 * 14, df.AgeTime = 1 * 14,
                    df.Temp = 6, df.Dew = 3, smooth = c("ns", "bs"),
                    control = glm.control(epsilon = 1e-10, maxit = 1000)) {
        smooth <- match.arg(smooth)
        dataname <- as.name(substitute(data))

        ## Parse formula to get cause of death and pollutant terms
        cause <- all.vars(formula)[1]
        pollterms <- attr(terms(formula), "term.labels")
        modelFormula <- setupFormula(cause, pollterms, df.Time, df.AgeTime,
                                     df.Temp, df.Dew, smooth)

        ## Fit the model!
        call <- substitute(glm(modelFormula, family = quasipoisson, data = data,
                               control = control, na.action = na.exclude),
                           list(modelFormula = modelFormula, data = dataname,
                                control = substitute(control)))
        fit <- eval.parent(call)
        structure(fit, class = c("tsModel", class(fit)), pollterms = pollterms)
}

setupFormula <- function(cause, pollutant, df.Time, df.AgeTime, df.Temp,
                         df.Dew, smooth) {
        f <- substitute(~ dow + agecat + smooth(tmpd, df.Temp)
                        + smooth(rmtmpd, df.Temp) + smooth(dptp, df.Dew)
                        + smooth(rmdptp, df.Dew) + smooth(time, df.Time)
                        + I(smooth(time, df.AgeTime) * Age2Ind)
                        + I(smooth(time, df.AgeTime) * Age3Ind),
                        list(df.Time = df.Time, df.AgeTime = df.AgeTime,
                             df.Temp = df.Temp, df.Dew = df.Dew,
                             smooth = as.name(smooth)))

        ## Tack on the pollutant(s) to the end of the formula
        rhs <- as.character(f)
        rhs[-1] <- paste(rhs[-1], paste(pollutant, collapse = "+"), sep = "+")

        ## Add LHS and coerce to `formula'
        as.formula(paste(cause, paste(rhs, collapse = "")))
}

