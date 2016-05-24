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

residual.diagnostics.menu <-
  function() {

    choices <- c("Return to previous menu ->",
                 "Distribution of weighted residuals (hist)",
                 "Distribution of weighted residuals (QQ)",
                 "Individual distributions of weighted residuals (hist)",
                 "Individual distributions of weighted residuals (QQ)",
                 "Absolute value of weighted residuals/IWRES vs predictions/IPRED",
                 "Absolute value of weighted residuals vs predictions",
                 "Absolute value of IWRES vs IPRED",
                 "Covariates vs absolute value of weighted residuals (BW)",
                 "Absolute value of weighted residuals vs pred|covariates",
                 "Absolute value of IWRES vs ipred|covariates",
                 "Autocorrelation of weighted residuals"
                 )
      title="\nRESIDUAL DIAGNOSTICS MENU\n  \\main\\goodness of fit plots\\Residual error model diagnostics"

    pick <- menu(choices,title=title)

    if(is.null(check.vars(c("cwres"),eval(parse(text=".cur.db")),silent=TRUE))) {
      wres <- "wres"
    }else{
      wres <- "cwres"
    }

    qx <- 0
    switch(pick+1,
           qx <- 2,
           qx <- 1,
           ##print(wres.dist.hist(eval(parse(text=".cur.db")))),
           print(eval(parse(text=paste(wres,".dist.hist(.cur.db)",sep="")))),

           ##print(wres.dist.qq(eval(parse(text=".cur.db")))),
           print(eval(parse(text=paste(wres,".dist.qq(.cur.db)",sep="")))),
           
           ##print(ind.plots.wres.hist(eval(parse(text=".cur.db")))),
           print(eval(parse(text=paste("ind.plots.",wres,".hist(.cur.db)",sep="")))),
           
           ##print(ind.plots.wres.qq(eval(parse(text=".cur.db")))),
           print(eval(parse(text=paste("ind.plots.",wres,".qq(.cur.db)",sep="")))),
           
           ##print(absval.iwres.wres.vs.ipred.pred(eval(parse(text=".cur.db")))),
           print(eval(parse(text=paste("absval.iwres.",wres,".vs.ipred.pred(.cur.db)",sep="")))),
           
           ##print(absval.wres.vs.pred(eval(parse(text=".cur.db")))),
           print(eval(parse(text=paste("absval.",wres,".vs.pred(.cur.db)",sep="")))),
           
           print(absval.iwres.vs.ipred(eval(parse(text=".cur.db")))),

           ##print(absval.wres.vs.cov.bw(eval(parse(text=".cur.db")),bins=9)),
           print(eval(parse(text=paste("absval.",wres,".vs.cov.bw(.cur.db),bins=9)",sep="")))),

           ##print(absval.wres.vs.pred.by.cov(eval(parse(text=".cur.db")))),
           print(eval(parse(text=paste("absval.",wres,".vs.pred.by.cov(.cur.db)",sep="")))),

           print(absval.iwres.vs.ipred.by.cov(eval(parse(text=".cur.db")))),

           ##print(autocorr.wres(eval(parse(text=".cur.db")))),
           print(eval(parse(text=paste("autocorr.",wres,"(.cur.db)",sep=""))))
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
