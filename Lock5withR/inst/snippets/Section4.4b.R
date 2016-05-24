Cocaine <- rbind( 
  do(18) * data.frame( treatment = "Lithium", response = "Relapse"),
  do(6)  * data.frame( treatment = "Lithium", response = "No Relapse"),
  do(20) * data.frame( treatment = "Placebo", response = "Relapse"),
  do(4)  * data.frame( treatment = "Placebo", response = "No Relapse")
  )

