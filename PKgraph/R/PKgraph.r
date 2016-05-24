#############################################################################################
## Project: PKgraph
## File: PKgraph.R
## Author: Xiaoyong Sun
## Date: 08/19/2009
## Goal: PKgraph
##        - interface
## Notes:
#############################################################################################


PKgraph <- function()
{

  mainHeight = getSubHeight()
  rightWidth = getSubWidth()

  assignInNamespace("PKW", gwindow(title="PKgraph", parent=c(100,50), height=mainHeight, width=rightWidth), "PKgraph")

  mb = gmenu(mbl)
  tb = gtoolbar(tbl)
  mainGroup = ggroup(horizontal=FALSE, spacing=0, cont=PKW, expand=TRUE)

###############################
## main layout
####################################

  add(mainGroup, mb)
  add(mainGroup, tb)
  bottomGroup = ggroup(horizontal=TRUE)
  add(mainGroup, bottomGroup, expand=TRUE)

  assignInNamespace("pmg.dialog.notebook",
                      gnotebook(closebuttons = TRUE, tearable = FALSE),
                      "PKgraph")
  size(pmg.dialog.notebook) <- c(rightWidth*0.6, mainHeight*0.8)

  assignInNamespace("pmg.dialog.notebook2",
                      gnotebook(closebuttons = TRUE, tearable = FALSE),
                      "PKgraph")
  size(pmg.dialog.notebook2) <- c(rightWidth*0.6, mainHeight*0.8)

  rightpane = gpanedgroup(pmg.dialog.notebook, horizontal=FALSE)

  assignInNamespace("pk.dirname", glabel(text=paste("Current directory: ", getwd())), "PKgraph")
  size(pk.dirname) <- c(rightWidth*0.3, mainHeight*0.05)
  assignInNamespace("pk.dir", gtable(as.character(dir()), sort.columns = 1:2, expand=TRUE, handler=openDataHandler), "PKgraph")

  leftpane=gpanedgroup(pk.dirname, pk.dir, horizontal=FALSE)
  pg = gpanedgroup(leftpane, rightpane)

  add(bottomGroup, pg, expand=TRUE)
  assignInNamespace("pmg.statusBar", gstatusbar("Ready", container=NULL), "PKgraph")
  add(mainGroup, pmg.statusBar)

  ## setup up value for main panel - right panel: notebook (right panel); dirname/dirtable (left panel)
  
  ## handler
  addhandlerclicked(pk.dir, handler=function(h,...)
                  {

                  })

}
