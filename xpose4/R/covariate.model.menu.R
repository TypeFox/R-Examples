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

"covariate.model.menu" <-
  function() {
  
    choices <- c("Return to previous menu ->",
                 "Numerically summarize the covariates",
                 "Scatterplot matrix of covariates",
                 "Parameters vs covariates",
                 "* Parameters vs covariates + model prediction",
                 "Weighted residuals vs covariates",
                 "GAM",
                 "Bootstrap of the GAM", "Bootstrap of the SCM",
                 "* Tree"
                 )

    title="\nCOVARIATE MODEL MENU\n  \\main\\covariate model"
    
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
           cov.summary(eval(parse(text=".cur.db")),out.file=".ask"),
           print(cov.splom(eval(parse(text=".cur.db")))),
           print(parm.vs.cov(eval(parse(text=".cur.db")))),
           cat("\nNot implemented yet!\n"),

           ##print(wres.vs.cov(eval(parse(text=".cur.db")))),
           print(eval(parse(text=paste(wres,".vs.cov(.cur.db)",sep="")))),

           qx <- gam.menu(),
           qx <- bootgam.menu(),
           qx <- bootscm.menu(),
           cat("\nNot implemented yet!\n")
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
