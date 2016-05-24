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

"cwres.menu" <- function() {

  choices <- c("Return to previous menu ->",
               "Compute CWRES",
               "Compare CWRES and weighted residuals",
               "Basic goodness of fit plots using CWRES",
               "CWRES vs independent variable",
               "CWRES vs independent variable (BW)",
               "CWRES vs PRED",
               "CWRES vs PRED (BW)",
               "Distribution of CWRES (hist)",
               "Distribution of CWRES (QQ)",
               "Individual distributions of CWRES (hist)",
               "Individual distributions of CWRES (QQ)",
               "Absolute value of CWRES vs predictions",
               "Covariates vs absolute value of CWRES (BW)",
               "Absolute value of CWRES vs pred|covariates",
               "Autocorrelation of CWRES"
               )

  title="\nCONDITIONAL WEIGHTED RESIDUALS (CWRES) MENU\n  \\main\\Conditional weighted residuals"

  pick <- menu(choices,title=title)

  qx <- 0
  switch(pick+1,
         qx <- 2,
         qx <- 1,
         xpose.calculate.cwres(eval(parse(text=".cur.db")),classic=TRUE),
         print(cwres.wres.vs.idv(eval(parse(text=".cur.db")))),
         print(basic.gof(eval(parse(text=".cur.db")))),
         print(cwres.vs.idv(eval(parse(text=".cur.db")))),
         print(cwres.vs.idv.bw(eval(parse(text=".cur.db")))),
         print(cwres.vs.pred(eval(parse(text=".cur.db")))),
         print(cwres.vs.pred.bw(eval(parse(text=".cur.db")))),
         print(cwres.dist.hist(eval(parse(text=".cur.db")))),
         print(cwres.dist.qq(eval(parse(text=".cur.db")))),
         print(ind.plots.cwres.hist(eval(parse(text=".cur.db")))),
         print(ind.plots.cwres.qq(eval(parse(text=".cur.db")))),
         print(absval.cwres.vs.pred(eval(parse(text=".cur.db")))),
         print(absval.cwres.vs.cov.bw(eval(parse(text=".cur.db")),bins=9)),
         print(absval.cwres.vs.pred.by.cov(eval(parse(text=".cur.db")))),
         print(autocorr.cwres(eval(parse(text=".cur.db"))))
         )
  
  if(qx == 2) {
    return(invisible(2))
  } else {
    if(qx == 1) {
      return(invisible(0))
    } else {
      Recall()
    }
  } 
  
  
}
