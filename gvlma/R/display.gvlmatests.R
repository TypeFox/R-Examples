"display.gvlmatests" <-
function(gvlmaobj)
{
  Global <- gvlmaobj$GlobalTest
  alphalevel <- Global$LevelOfSignificance
  TestVals <- Global[grep("Stat", names(Global))]
  GlobalTests <- do.call("rbind",
                         lapply(TestVals,
                                function(x) as.data.frame(x))
                         )
  names(GlobalTests) <- list("Value", "p-value", "Decision")
  row.names(GlobalTests) <- c("Global Stat",
                              "Skewness",
                              "Kurtosis",
                              "Link Function",
                              "Heteroscedasticity")
  dec <- c("Assumptions acceptable.", 
           "Assumptions NOT satisfied!")
  GlobalTests$Decision <- dec[GlobalTests$Decision + 1]
  cat("\nASSESSMENT OF THE LINEAR MODEL ASSUMPTIONS")
  cat("\nUSING THE GLOBAL TEST ON 4 DEGREES-OF-FREEDOM:")
  cat("\nLevel of Significance = ", alphalevel,"\n")
  cat("\nCall:\n", deparse(gvlmaobj$GlobalTest$call), "\n\n")
  print(GlobalTests, digits = max(3, getOption("digits") - 3))
}

