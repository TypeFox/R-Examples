print.tinpoAllBoots <-
function (x,...) {

  res.tinpoAllBoots <- x
  if (!inherits(res.tinpoAllBoots, "tinpoAllBoots")) stop ("no convenient data")
  cat("**TInPosition All Bootstrap output data**\n")
  cat("*Contains the following objects:\n\n")
  res <- array("", c(2, 2), list(1:2, c("name", "description")))
  
  res[1,] <- c("$fj.boot.data","Bootstrap data associated to ($fj; measures).")
  res[2,] <- c("$fi.boot.data","Bootstrap data associated to ($fi; groups).")
  
  print(res)

}
