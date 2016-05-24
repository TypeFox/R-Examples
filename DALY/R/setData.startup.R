## Create 'data' window

setData.startup <-
function(n){
  ## Define some variables
  ar <- c(1, 5, 2, 6, 3, 7, 4, 8)

  countDist <- c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
  propDist <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)
  timeDist <- c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)

  selectDist <- cbind(countDist, propDist, timeDist, timeDist,
                      propDist, propDist, countDist, timeDist)

  ## Define active window state
  win <- paste("data", n, sep = "")
  DALYassign("active.windows", TRUE, item = win)

  ## Create 'data[n]' window
  DALYassign(win, tktoplevel(padx = 2, pady = 2))
  drawWindow(DALYget(win), DALYget(paste("outcome", n, "Name", sep = "")), n)

  ## Catch closing of main window
  tkbind(DALYget(win), "<Destroy>", function(){ DALYdestroy(win) })
  
  for (i in seq(8)){
    ## First update Tcl objects
    DALYupdate(paste(".dist", DALYget("txtLbl")[i], n, sep = ""))
	DALYupdate(paste(".strat", DALYget("txtLbl")[i], n, sep = ""))
	DALYupdate(paste(".", DALYget("txtlbl")[i], n, sep = ""))
	
	## Define 'dist' & 'strat' comboboxes
    DALYassign(paste("cbDistr", DALYget("txtLbl")[i], n, sep = ""), 
           defineCombo(DALYget(paste("Frame", ar[i], sep = "")),
                  DALYget(paste(".dist", DALYget("txtLbl")[i], n, sep = "")),
                  DALYget("distributions")[selectDist[, i]]))

    DALYassign(paste("cbStrat", DALYget("txtLbl")[i], n, sep = ""), 
           defineCombo(DALYget(paste("Frame", ar[i], sep = "")),
                  DALYget(paste(".strat", DALYget("txtLbl")[i], n, sep = "")),
                  DALYget("stratifications")))

	## Define table
    DALYassign(paste("table", DALYget("txtLbl")[i], n, sep = ""),
           defineTable(DALYget(paste("Frame", ar[i], sep = "")),
                  DALYget(paste(".", DALYget("txtlbl")[i], n, sep = ""))))
  }

  bindCombo(DALYget(paste("cbDistrInc", n, sep = "")),
            DALYget(paste("cbDistrOns", n, sep = "")),
            DALYget(paste("cbDistrDWt", n, sep = "")),
            DALYget(paste("cbDistrMrt", n, sep = "")), 
            DALYget(paste("cbDistrTrt", n, sep = "")),
            DALYget(paste("cbDistrDur", n, sep = "")),
            DALYget(paste("cbDistrDWn", n, sep = "")),
            DALYget(paste("cbDistrLxp", n, sep = "")), 
            DALYget(paste(".distInc", n, sep = "")),
            DALYget(paste(".distOns", n, sep = "")),
            DALYget(paste(".distDWt", n, sep = "")),
            DALYget(paste(".distMrt", n, sep = "")), 
            DALYget(paste(".distTrt", n, sep = "")),
            DALYget(paste(".distDur", n, sep = "")),
            DALYget(paste(".distDWn", n, sep = "")),
            DALYget(paste(".distLxp", n, sep = "")), 
            DALYget(paste("cbStratInc", n, sep = "")),
            DALYget(paste("cbStratOns", n, sep = "")),
            DALYget(paste("cbStratDWt", n, sep = "")),
            DALYget(paste("cbStratMrt", n, sep = "")), 
            DALYget(paste("cbStratTrt", n, sep = "")),
            DALYget(paste("cbStratDur", n, sep = "")),
            DALYget(paste("cbStratDWn", n, sep = "")),
            DALYget(paste("cbStratLxp", n, sep = "")), 
            DALYget(paste(".stratInc", n, sep = "")),
            DALYget(paste(".stratOns", n, sep = "")),
            DALYget(paste(".stratDWt", n, sep = "")),
            DALYget(paste(".stratMrt", n, sep = "")), 
            DALYget(paste(".stratTrt", n, sep = "")),
            DALYget(paste(".stratDur", n, sep = "")),
            DALYget(paste(".stratDWn", n, sep = "")),
            DALYget(paste(".stratLxp", n, sep = "")), 
            DALYget(paste("tableInc", n, sep = "")),
            DALYget(paste("tableOns", n, sep = "")),
            DALYget(paste("tableDWt", n, sep = "")),
            DALYget(paste("tableMrt", n, sep = "")), 
            DALYget(paste("tableTrt", n, sep = "")),
            DALYget(paste("tableDur", n, sep = "")),
            DALYget(paste("tableDWn", n, sep = "")),
            DALYget(paste("tableLxp", n, sep = "")), 
            DALYget(paste(".inc", n, sep = "")),
            DALYget(paste(".ons", n, sep = "")),
            DALYget(paste(".DWt", n, sep = "")),
            DALYget(paste(".mrt", n, sep = "")), 
            DALYget(paste(".trt", n, sep = "")),
            DALYget(paste(".dur", n, sep = "")),
            DALYget(paste(".DWn", n, sep = "")),
            DALYget(paste(".lxp", n, sep = "")))

  drawTables(DALYget(paste("tableInc", n, sep = "")),
             DALYget(paste("tableOns", n, sep = "")),
             DALYget(paste("tableDWt", n, sep = "")), 
             DALYget(paste("tableMrt", n, sep = "")),
             DALYget(paste("tableTrt", n, sep = "")),
             DALYget(paste("tableDur", n, sep = "")), 
             DALYget(paste("tableDWn", n, sep = "")),
             DALYget(paste("tableLxp", n, sep = "")),
             DALYget(paste("cbDistrInc", n, sep = "")), 
             DALYget(paste("cbDistrOns", n, sep = "")),
             DALYget(paste("cbDistrDWt", n, sep = "")),
             DALYget(paste("cbDistrMrt", n, sep = "")), 
             DALYget(paste("cbDistrTrt", n, sep = "")),
             DALYget(paste("cbDistrDur", n, sep = "")),
             DALYget(paste("cbDistrDWn", n, sep = "")), 
             DALYget(paste("cbDistrLxp", n, sep = "")),
             DALYget(paste("cbStratInc", n, sep = "")),
             DALYget(paste("cbStratOns", n, sep = "")), 
             DALYget(paste("cbStratDWt", n, sep = "")),
             DALYget(paste("cbStratMrt", n, sep = "")),
             DALYget(paste("cbStratTrt", n, sep = "")), 
             DALYget(paste("cbStratDur", n, sep = "")),
             DALYget(paste("cbStratDWn", n, sep = "")),
             DALYget(paste("cbStratLxp", n, sep = "")))
}