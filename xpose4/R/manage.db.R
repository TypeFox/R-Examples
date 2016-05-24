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

"manage.db"<-
  function(){
    
    ## Manage databases
    
    choices <- c("Return to the previous menu ->",
                 "Change run number/database", 
                 "List items in current database", 
                 "View the documentation for the current database",
                 "Change the documentation for the current database",
                 "Change xpose defined variables (id, idv, dv, etc.)",
                 "Change independent variable (idv)", 
                 "Change parameter scope", 
                 "Change covariate scope",
                 "Export variable definitions to a file",
                 "Import variable definitions from a file",
                                        #"Change ID variable", 
                                        #"Change dependent variable (DV)",
                                        #"Change PRED variable",
                                        #"Change IPRED variable",
                                        #"Change RES variable",
                                        #"Change WRES variable",
                                        #"Change IWRES variable",
                                        #"Change occasion variable",
                                        #"Change variable to label data points with", 
                 "Change the subset expression",
                 "Change missing data variable (MDV)",
                 "Change category threshold for variables",
                 "Change category threshold for DV",
                 "Change name of a variable",
                                        #"Change typical parameter scope",
                 "* Map random effects to typical parameters", 
                 "* Map typical parameters to parameters",
                 "* Copy variable definitions from another database"
                 )
    
    title="\nDATABASE MANAGEMENT MENU\n  \\main\\preferences\\Manage variables in the current database"
    
    pick <- menu(choices,title=title)
    
    qx <- 0
    switch(pick+1,
           qx <- 2,
           qx <- 1,
           change.xp.obj(),
           db.names(eval(parse(text=".cur.db"))),
           get.doc(eval(parse(text=".cur.db")), classic = T),
           set.doc(eval(parse(text=".cur.db")), classic = T),
           change.xvardef(eval(parse(text=".cur.db")), classic=T),
           change.xvardef(eval(parse(text=".cur.db")), classic=T, var="idv"),
           change.xvardef(eval(parse(text=".cur.db")), classic=T, var="parm"),
           change.xvardef(eval(parse(text=".cur.db")), classic=T, var="covariates"),
           export.variable.definitions(eval(parse(text=".cur.db"))),
           import.variable.definitions(eval(parse(text=".cur.db")), classic = T),
                                        #change.id(eval(parse(text=".cur.db")), classic = T),
                                        #change.idv(eval(parse(text=".cur.db")), classic = T),
                                        #change.dv(eval(parse(text=".cur.db")), classic = T),
                                        #change.pred(eval(parse(text=".cur.db")), classic = T),
                                        #change.ipred(eval(parse(text=".cur.db")), classic = T),
                                        #change.res(eval(parse(text=".cur.db")), classic = T),
                                        #change.wres(eval(parse(text=".cur.db")), classic = T),
                                        #change.iwres(eval(parse(text=".cur.db")), classic = T),
                                        #change.occ(eval(parse(text=".cur.db")), classic = T),
                                        #change.label(eval(parse(text=".cur.db")), classic = T),
           change.subset(eval(parse(text=".cur.db")), classic = T),
           change.miss(eval(parse(text=".cur.db")), classic = T),
           change.cat.levels(eval(parse(text=".cur.db")), classic = T),
           change.dv.cat.levels(eval(parse(text=".cur.db")), classic = T),
           change.var.name(eval(parse(text=".cur.db")), classic = T),
                                        #change.covs(eval(parse(text=".cur.db")), classic = T),
                                        #change.tvparms(eval(parse(text=".cur.db")), classic = T),
           cat("Not yet implemented!\n"),  # change.ranpar(eval(parse(text=".cur.db")), classic = T),
           cat("Not yet implemented!\n"),  # change.tvpar(eval(parse(text=".cur.db")), classic = T),
           cat("Not yet implemented!\n")  # copy.attr(eval(parse(text=".cur.db")), classic = T)
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

