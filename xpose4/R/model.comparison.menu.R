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

model.comparison.menu <-
  function() {

    choices <- c("Return to previous menu ->",
                 "Read in comparison model to Xpose",
                 "Read in second comparison model to Xpose",
                 "Basic model comparisons",
                 "Additional model comparisons",
                 "Delta OFV vs ID",
                 "Delta OFV vs Covariates",
                 "Delta OFV1 vs Delta OFV2",
                 "Delta PRED/IPRED/Weighted residuals vs covariates ->"
                 )

    title="\nMODEL COMPARISON MENU\n  \\main\\model comparison"
    
    pick <- menu(choices,title=title)

    if(is.null(check.vars(c("cwres"),eval(parse(text=".cur.db")),silent=TRUE))) {
      wres <- ""
    }else{
      wres <- ".cwres"
    }

    ref.db <- NULL
    if(exists(".ref.db")) ref.db <- eval(parse(text=".ref.db"))
    ref.db2 <- NULL
    if(exists(".ref.db2")) ref.db2 <- eval(parse(text=".ref.db2"))


    run.basic.model.comp  <- function(){
      cat("\nRunning command:\n",
          "basic.model.comp(.cur.db,object.ref=ref.db)\n",
          sep="")
      print(basic.model.comp(eval(parse(text=".cur.db")),object.ref=ref.db))
    }
    
    
    qx <- 0
    switch(pick+1,
           qx <- 2,
           qx <- 1,
           ref.db <- get.refrunno(),
           ref.db2 <- get.refrunno(database=".ref.db2"),
           run.basic.model.comp(),
           ##print(basic.model.comp(eval(parse(text=".cur.db")),object.ref=ref.db)),
           ##print(add.model.comp(eval(parse(text=".cur.db")))),
           print(eval(parse(text=paste("add.model.comp",
                              wres,"(.cur.db,object.ref=ref.db)",sep="")))),
           print(dOFV.vs.id(eval(parse(text=".cur.db")),ref.db)),
           print(dOFV.vs.cov(eval(parse(text=".cur.db")),ref.db)),
           print(dOFV1.vs.dOFV2(eval(parse(text=".cur.db")),ref.db,ref.db2)),
           qx <- model.comparison.covariates.menu()
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
