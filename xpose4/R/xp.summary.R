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

"xp.summary" <-
  function(gamobj=NULL) {

      if(is.null(gamobj)){
          gamobj <- check.gamobj()
          if(is.null(gamobj)){
              return()
          } else {
          }
      } else {
          c1 <- call("assign",pos=1, "current.gam", gamobj,immediate=T)
          eval(c1)
      }

    cat("\nSUMMARY")
    print(summary(eval(parse(text="current.gam"))))

    cat("\nPATH TO FINAL MODEL\n")
    print(eval(parse(text="current.gam$anova")))

    cat("\nCOEFFICIENTS\n")
    print(coefficients(eval(parse(text="current.gam"))))

    cat("\nPRERUN RESULTS\n")
    cat("Dispersion:",eval(parse(text="current.gam$dispersion")),"\n")

    cat("\nDATA\n")
    cat("Subset expression:",eval(parse(text="current.gam$subset")),"\n")
    cat("Only first value of covariate considered\n")
    cat("for each individual:",eval(parse(text="current.gam$onlyfirst")),"\n")
    cat("Covariates normalized to median:",eval(parse(text="current.gam$medianNorm")),"\n")

    return(invisible())

  }
