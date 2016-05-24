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

"create.xpose.plot.classes" <- function()
{
  
  ## setClassUnion("character or NULL",c("character","NULL"))
  ## setClassUnion("character or numeric",c("character","numeric"))
  ## setClassUnion("numeric or NULL",c("numeric","NULL"))
  ## setClassUnion("data.frame or NULL",c("data.frame","NULL"))
  ## setClassUnion("list or NULL",c("list","NULL"))
  ## setClassUnion("lang or numeric",c("vector","numeric","list"))
  ## setClassUnion("logical or numeric",c("logical","numeric"))

  setClass("xpose.multiple.plot",where=.GlobalEnv,
           representation(plotList           = "list or NULL",
                          plotTitle          = "character or NULL",
                          prompt             = "logical",
                          new.first.window   = "logical",
                          max.plots.per.page = "numeric",
                          title              = "list",
                          ##title.x            = "xptmp",
                          ##title.y            = "lang or numeric",
                          ##title.just         = "vector",
                          ##title.gp           = "lang or numeric",
                          mirror             = "logical or numeric",
                          bql.layout         = "logical"
                          ),
           prototype(plotList = NULL,
                     plotTitle= NULL,
                     prompt   = FALSE,
                     new.first.window   = FALSE,
                     max.plots.per.page = 4,
                     title    = list(
                       title.x = unit(0.5, "npc"),
                       title.y = unit(0.5, "npc"),
                       title.gp= gpar(cex=1.2,fontface="bold"),#,font=2),
                       title.just = c("center","center")
                       ),
                     ##title.x            = unit(0.5, "npc"),
                     ##title.y            = c(unit(0.5, "npc")),
                     ##title.just         = c("center","center"),
                     ##title.gp           = list(cex=1.2,fontface="bold",font=2),
                     mirror             = FALSE,
                     bql.layout         = FALSE
                     )
           )
  

    ## Define the methods
  setMethod("print",signature(x="xpose.multiple.plot"),print.xpose.multiple.plot)
  setMethod("show","xpose.multiple.plot",function(object) print(x=object))
  invisible()
  
}
