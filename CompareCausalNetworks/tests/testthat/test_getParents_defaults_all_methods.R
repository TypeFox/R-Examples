library(CompareCausalNetworks)
context("All supported methods")

data("simData_unknownShiftInterventions")

X <- simData_unknownShiftInterventions$X
environment <- simData_unknownShiftInterventions$environment

methods <- c("ICP", "hiddenICP", "backShift", "pc", "LINGAM",
             "ges", "CAM", "rfci", "regression",
             "bivariateANM", "bivariateCAM")


# TODO: change all method names to spelling in original package?

for(method in methods){
  test_that(paste("Checks output type for", method), {
    
    expect_is(
      Ahat <- getParents(X, environment, method=method, alpha=0.1)
      , "Matrix")
    
    
    if(method %in% c("ICP", "hiddenICP")){
      expect_warning(
        Ahat <- getParents(X, environment, method=method, alpha=0.1)
      )
    }
    
    
  }
  )
}
# gies