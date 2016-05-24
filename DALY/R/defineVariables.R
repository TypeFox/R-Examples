## Set standard values

defineVariables <-
function(){
  DALYassign("ageGroups", c("0-4", "5-14", "15-44", "45-59", "60+"))
  DALYassign("fixed", c("Age Group", "Male", "Female"))
  DALYassign("txtLBL", c("INCIDENCE", "TREATMENT", "ONSET",
                         "DURATION", "DWtreated", "DWuntreated",
                         "MORTALITY", "AvgAgeDeath"))
  DALYassign("txtLbl", c("Inc", "Trt", "Ons", "Dur",
                         "DWt", "DWn", "Mrt", "Lxp"))
  DALYassign("txtlbl", c("inc", "trt", "ons", "dur",
                         "DWt", "DWn", "mrt", "lxp"))
  DALYassign("distributions", c("Beta-Pert", "Beta", "Gamma", "Normal",
                                "LogNormal.geom", "LogNormal.arithm",
                                "Uniform", "Fixed"))
  DALYassign("stratifications", c("Age and Sex", "Age", "Sex", "None"))
  DALYassign("stdM", c(80.00, 79.36, 75.38, 70.40, 65.41, 60.44, 55.47,
                       50.51, 45.57, 40.64, 35.77, 30.99, 26.32, 21.81,
                       17.50, 13.58, 10.17, 7.45, 5.24, 3.54, 2.31))
  DALYassign("stdF", c(82.50, 81.84, 77.95, 72.99, 68.02, 63.08, 58.17,
                       53.27, 48.38, 43.53, 38.72, 33.99, 29.37, 24.83,
                       20.44, 16.20, 12.28, 8.90, 6.22, 4.25, 2.89))
  DALYassign("ages", c(0, 1, 5, 10, 15, 20, 25,
                       30, 35, 40, 45, 50, 55, 60,
                       65, 70, 75, 80, 85, 90, 95))

  ## 'pop' = Population matrix
  .pop <- tclArray()
  for(x in seq(0, 2))
    .pop[[0, x]] <- DALYget("fixed")[x+1]
  for(y in seq(5))
    .pop[[y, 0]] <- DALYget("ageGroups")[y]
  DALYassign(".pop", .pop)
  DALYassign("pop", matrix(nrow = 5, ncol = 2))

  ## 'LE' = Life Expectancy table
  .LE <- tclArray()
  .LE[[0, 0]] <- "Age"
  for(x in seq(2)) .LE[[0, x]] <- DALYget("fixed")[x+1]
  for(y in seq(21)) .LE[[y, 0]] <- DALYget("ages")[y]
  DALYassign(".LE", .LE)

  DALYassign("LE", matrix(nrow = 21, ncol = 2))
  setStdLE()

  ## Set 'data', 'dist' & 'strat'
  distList <- c(3, 2, 8, 8, 2, 2, 3, 8)
  for (i in seq(8)){
    for (j in seq(8)){
      ## 'txtlbl' + 'i' = parameters per outcome
      DALYassign(paste(".", DALYget("txtlbl")[j], i, sep = ""),
                 tclArray())
      DALYassign(paste(DALYget("txtlbl")[j], i, sep = ""),
                 matrix(nrow = 6, ncol = 5))

      ## 'strat' + 'txtLbl' + 'i' = stratification per parameter per outcome
      DALYassign(paste(".strat", DALYget("txtLbl")[j], i, sep = ""),
                 tclVar(DALYget("stratifications")[1]))
      DALYassign(paste("strat", DALYget("txtLbl")[j], i, sep = ""),
                 DALYget("stratifications")[1])

      ## 'distributions' + 'txtLbl' + 'i' = dist per parameter per outcome
      d <- distList[j]
      DALYassign(paste(".dist", DALYget("txtLbl")[j], i, sep = ""),
                 tclVar(DALYget("distributions")[d]))
      DALYassign(paste("dist", DALYget("txtLbl")[j], i, sep = ""),
                 DALYget("distributions")[d])
    }

    ## assign 'ageGroups' labels to '.txtlbl' + 'i'
    for (j in seq(5)){
      for(k in seq(8)){
        DALYeval(parse(text = paste(".", DALYget("txtlbl")[k], i,
                                    "[[", j, ",0]] <- '",
                                    DALYget("ageGroups")[j], "'", sep = "")))
      }
    }
  }

  ## assign 'disease' and 'outcome' names
  DALYassign("diseaseName", tclVar())
  for (i in seq(8))
    DALYassign(paste("outcome", i, "Name", sep = ""), tclVar())

  ## assign 'aw' and 'dr' variables
  DALYassign(".aw", tclVar("No"))
  DALYassign(".dr", tclVar("0"))

  ## assign 'option' variables
  DALYassign(".it", tclVar("20000"))
  DALYassign(".optOP", tclVar("Summed over age/sex classes"))
  DALYassign(".optOC", tclVar("Summed over outcomes"))
  DALYassign(".optRA", tclVar("Absolute"))
  DALYassign(".optHist", tclVar("1"))
  DALYassign("it", 20000)
  DALYassign("optOP", "Summed over age/sex classes")
  DALYassign("optOC", "Summed over outcomes")
  DALYassign("optRA", "Absolute")
  DALYassign("optHist", 1)
}