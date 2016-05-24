## Helper function for 'bindCombo'
## Update data table layout, based on selected dist/strat

setDist <-
function(dist, strat, tbl, var){
  d <- which(DALYget("distributions") == tclvalue(dist))
  s <- which(DALYget("stratifications") == tclvalue(strat))
  for (i in seq(6))
    var[[0, i]] <- NULL
  for (i in seq(5))
    var[[i, 0]] <- DALYget("ageGroups")[i]

  tkfocus(DALYget("Top"))
  for (x in seq(5)){
    for (y in seq(6)){
      tcl(.Tk.ID(tbl), "tag", "celltag", "alive",
          paste(x, ", ", y, sep = ""))
      tcl(.Tk.ID(tbl), "tag", "configure", "alive",
          state = "normal", bg = "white")
    }
  }

  if (d == 1){  # Beta-Pert
    if (s == 1){  # Age and Sex
      var[[0, 1]] <- "M.Mode"
      var[[0, 2]] <- "M.Min"
      var[[0, 3]] <- "M.Max"
      var[[0, 4]] <- "F.Mode"
      var[[0, 5]] <- "F.Min"
      var[[0, 6]] <- "F.Max"
    }
    if (s == 2){  # Age
      var[[0, 1]] <- "Mode"
      var[[0, 2]] <- "Min"
      var[[0, 3]] <- "Max"
      kill(tbl, var, c(1:5), c(4:6))
    }
    if (s == 3){  # Sex
      var[[0, 1]] <- "M.Mode"
      var[[0, 2]] <- "M.Min"
      var[[0, 3]] <- "M.Max"
      var[[0, 4]] <- "F.Mode"
      var[[0, 5]] <- "F.Min"
      var[[0, 6]] <- "F.Max"
      noGroups(var)
      kill(tbl, var, c(2:5), c(1:6))
    }
    if (s == 4){  # None
      var[[0, 1]] <- "Mode"
      var[[0, 2]] <- "Min"
      var[[0, 3]] <- "Max"
      noGroups(var)
      kill(tbl, var, c(1:5), c(4:6))
      kill(tbl, var, c(2:5), c(1:6))
    }
  }

  if (d == 2){  # Beta
    kill(tbl, var, c(1:5), c(5:6))
    if (s == 1){  # Age and Sex
      var[[0, 1]] <- "M.Alpha"
      var[[0, 2]] <- "M.Beta"
      var[[0, 3]] <- "F.Alpha"
      var[[0, 4]] <- "F.Beta"
    }
    if (s == 2){  # Age
      var[[0, 1]] <- "Alpha"
      var[[0, 2]] <- "Beta"
      kill(tbl, var, c(1:5), c(3:4))
    }
    if (s == 3){  # Sex
      var[[0, 1]] <- "M.Alpha"
      var[[0, 2]] <- "M.Beta"
      var[[0, 3]] <- "F.Alpha"
      var[[0, 4]] <- "F.Beta"
      noGroups(var)
      kill(tbl, var, c(2:5), c(1:6))
    }
    if (s == 4){  # None
      var[[0, 1]] <- "Alpha"
      var[[0, 2]] <- "Beta"
      noGroups(var)
      kill(tbl, var, c(1:5), c(3:4))
      kill(tbl, var, c(2:5), c(1:6))
    }
  }

  if (d == 3){  # Gamma
    kill(tbl, var, c(1:5), c(5:6))
    if (s == 1){  # Age and Sex
      var[[0, 1]] <- "M.Shape"
      var[[0, 2]] <- "M.Rate"
      var[[0, 3]] <- "F.Shape"
      var[[0, 4]] <- "F.Rate"
    }
    if (s == 2){  # Age
      var[[0, 1]] <- "Shape"
      var[[0, 2]] <- "Rate"
      kill(tbl, var, c(1:5), c(3:4))
    }
    if (s == 3){  # Sex
      var[[0, 1]] <- "M.Shape"
      var[[0, 2]] <- "M.Rate"
      var[[0, 3]] <- "F.Shape"
      var[[0, 4]] <- "F.Rate"
      noGroups(var)
      kill(tbl, var, c(2:5), c(1:6))
    }
    if (s == 4){  # None
      var[[0, 1]] <- "Shape"
      var[[0, 2]] <- "Rate"
      noGroups(var)
      kill(tbl, var, c(1:5), c(3:4))
      kill(tbl, var, c(2:5), c(1:6))
    }
  }

  if (d == 4){  # Normal
    kill(tbl, var, c(1:5), c(5:6))
    if (s == 1){  # Age and Sex
      var[[0, 1]] <- "M.Mean"
      var[[0, 2]] <- "M.Sigma"
      var[[0, 3]] <- "F.Mean"
      var[[0, 4]] <- "F.Sigma"
    }
    if (s == 2){  # Age
      var[[0, 1]] <- "Mean"
      var[[0, 2]] <- "Sigma"
      kill(tbl, var, c(1:5), c(3:4))
    }
    if (s == 3){  # Sex
      var[[0, 1]] <- "M.Mean"
      var[[0, 2]] <- "M.Sigma"
      var[[0, 3]] <- "F.Mean"
      var[[0, 4]] <- "F.Sigma"
      noGroups(var)
      kill(tbl, var, c(2:5), c(1:6))
    }
    if (s == 4){  # None
      var[[0, 1]] <- "Mean"
      var[[0, 2]] <- "Sigma"
      noGroups(var)
      kill(tbl, var, c(1:5), c(3:4))
      kill(tbl, var, c(2:5), c(1:6))
    }
  }

  if (d == 5){  # LogNormal.geom
    kill(tbl, var, c(1:5), c(5:6))
    if (s == 1){  # Age and Sex
      var[[0, 1]] <- "M.logMu"
      var[[0, 2]] <- "M.logSigma"
      var[[0, 3]] <- "F.logMu"
      var[[0, 4]] <- "F.logSigma"
    }
    if (s == 2){  # Age
      var[[0, 1]] <- "logMean"
      var[[0, 2]] <- "logSigma"
      kill(tbl, var, c(1:5), c(3:4))
    }
    if (s == 3){  # Sex
      var[[0, 1]] <- "M.logMu"
      var[[0, 2]] <- "M.logSigma"
      var[[0, 3]] <- "F.logMu"
      var[[0, 4]] <- "F.logSigma"
      noGroups(var)
      kill(tbl, var, c(2:5), c(1:6))
    }
    if (s == 4){  # None
      var[[0, 1]] <- "logMean"
      var[[0, 2]] <- "logSigma"
      noGroups(var)
      kill(tbl, var, c(1:5), c(3:4))
      kill(tbl, var, c(2:5), c(1:6))
    }
  }

  if (d == 6){  # LogNormal.arithm
    kill(tbl, var, c(1:5), c(5:6))
    if (s == 1){  # Age and Sex
      var[[0, 1]] <- "M.Mean"
      var[[0, 2]] <- "M.Sigma"
      var[[0, 3]] <- "F.Mean"
      var[[0, 4]] <- "F.Sigma"
    }
    if (s == 2){  # Age
      var[[0, 1]] <- "Mean"
      var[[0, 2]] <- "Sigma"
      kill(tbl, var, c(1:5), c(3:4))
    }
    if (s == 3){  # Sex
      var[[0, 1]] <- "M.Mean"
      var[[0, 2]] <- "M.Sigma"
      var[[0, 3]] <- "F.Mean"
      var[[0, 4]] <- "F.Sigma"
      noGroups(var)
      kill(tbl, var, c(2:5), c(1:6))
    }
    if (s == 4){  # None
      var[[0, 1]] <- "Mean"
      var[[0, 2]] <- "Sigma"
      noGroups(var)
      kill(tbl, var, c(1:5), c(3:4))
      kill(tbl, var, c(2:5), c(1:6))
    }
  }

  if (d == 7){  # Uniform
    kill(tbl, var, c(1:5), c(5:6))
    if (s == 1){  # Age and Sex
      var[[0, 1]] <- "M.Min"
      var[[0, 2]] <- "M.Max"
      var[[0, 3]] <- "F.Min" 
      var[[0, 4]] <- "F.Max"
    }
    if (s == 2){  # Age
      var[[0, 1]] <- "Min"
      var[[0, 2]] <- "Max"
      kill(tbl, var, c(1:5), c(3:4))
    }
    if (s == 3){  # Sex
      var[[0, 1]] <- "M.Min"
      var[[0, 2]] <- "M.Max"
      var[[0, 3]] <- "F.Min"
      var[[0, 4]] <- "F.Max"
      noGroups(var)
      kill(tbl, var, c(2:5), c(1:6))
    }
    if (s == 4){  # None
      var[[0, 1]] <- "Min"
      var[[0, 2]] <- "Max"
      noGroups(var)
      kill(tbl, var, c(1:5), c(3:4))
      kill(tbl, var, c(2:5), c(1:6))
    }
  }

  if (d == 8){  # Fixed
    kill(tbl, var, c(1:5), c(3:6))
    if (s == 1){  # Age and Sex
      var[[0, 1]] <- "Male"
      var[[0, 2]] <- "Female"
    }
    if (s == 2){  # Age
      var[[0, 1]] <- "Value"
      kill(tbl, var, c(1:5), 2)
    }
    if (s == 3){  # Sex
      var[[0, 1]] <- "Male"
      var[[0, 2]] <- "Female"
      noGroups(var)
      kill(tbl, var, c(2:5), c(1:6))
    }
    if (s == 4){  # None
      var[[0, 1]] <- "Value"
      noGroups(var)
      kill(tbl, var, c(1:5), 2)
      kill(tbl, var, c(2:5), c(1:6))
    }
  }
}