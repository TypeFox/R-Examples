## Set standard values

resetVariables <-
function(){
  ## 'pop' = Population matrix
  DALYassign("pop", matrix(nrow = 5, ncol = 2))
  DALYupdate(".pop")

  ## Set 'data', 'dist' & 'strat'
  distList <- c(3, 2, 8, 8, 2, 2, 3, 8)
  for (i in seq(8)){
    for (j in seq(8)){
      ## 'txtlbl' + 'i' = parameters per outcome
      DALYassign(paste(DALYget("txtlbl")[j], i, sep = ""),
                 matrix(nrow = 6, ncol = 5))
      DALYupdate(paste(".", DALYget("txtlbl")[j], i, sep = ""))

      ## 'strat' + 'txtLbl' + 'i' = stratification per parameter per outcome
      DALYassign(paste("strat", DALYget("txtLbl")[j], i, sep = ""),
                 DALYget("stratifications")[1])
      DALYupdate(paste(".strat", DALYget("txtLbl")[j], i, sep = ""))

      ## 'distributions' + 'txtLbl' + 'i' = dist per parameter per outcome
      d <- distList[j]
      DALYassign(paste("dist", DALYget("txtLbl")[j], i, sep = ""),
                 DALYget("distributions")[d])
      DALYupdate(paste(".dist", DALYget("txtLbl")[j], i, sep = ""))
    }
  }

  ## assign 'disease' and 'outcome' names
  DALYupdate("diseaseName", "")
  for (i in seq(8))
    DALYupdate(paste("outcome", i, "Name", sep = ""), "")

  ## assign 'aw' and 'dr' variables
  DALYupdate(".aw", "No")
  DALYupdate(".dr", "0")
}