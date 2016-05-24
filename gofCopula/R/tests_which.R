gofWhich = function(copula, d){
  if (d <= 1){stop("The dimension d must be 2 or even higher.")}
  if(is.element(copula, c("independence", "gaussian", "t", "clayton", "frank", "gumbel", "joe", 
                          "amh", "survival clayton", "survival gumbel", "survival joe", "90 clayton", 
                          "90 gumbel", "90 joe", "270 clayton", "270 gumbel", "270 joe", 
                          "bb1", "bb6", "bb7", "bb8", "survival bb1", "survival bb6", 
                          "survival bb7", "survival bb8", "90 bb1", "90 bb6", "90 bb7", 
                          "90 bb8", "270 bb1", "270 bb6", "270 bb7", "270 bb8")) == F){stop("This copula is for no test implemented.")}
  
  which1 = which2 = which3 = which4 = which5 = c()
  if (d >= 2 & is.element(copula, c("gaussian", "t", "clayton", "frank", "gumbel"))){#, "amh", "joe"))){
    which1 = c("gofRosenblattSnB", "gofRosenblattSnC", "gofADChisq", "gofADGamma", "gofSn")
  }
  if (d == 2 & is.element(copula, c("gaussian", "t", "clayton", "frank", "gumbel"))){#, "amh", "joe"))){
    which2 = c("gofRn")
  }
  if (d == 2 & is.element(copula, c("gaussian", "t", "clayton", "frank", "gumbel"))){
    which3= c("gofPIOSRn", "gofPIOSTn", "gofKernel")
  }
  if (d == 2 & is.element(copula, c(#"independence", 
                                    "gaussian", "t", "clayton", "frank", "gumbel"))){#, 
#                                    "joe", "survival clayton", "survival gumbel", "survival joe", "90 clayton", 
 #                                   "90 gumbel", "90 joe", "270 clayton", "270 gumbel", "270 joe"))){
    which4 = c("gofWhite")
  }
  if (d == 2 & is.element(copula, c("independence", 
                                    "gaussian", "clayton", "frank", "gumbel"))){#, 
#                                    "joe", "survival clayton", "survival gumbel", "survival joe", "90 clayton", 
#                                    "90 gumbel", "90 joe", "270 clayton", "270 gumbel", "270 joe", 
#                                    "bb1", "bb6", "bb7", "bb8", "survival bb1", "survival bb6", 
#                                    "survival bb7", "survival bb8", "90 bb1", "90 bb6", "90 bb7", 
#                                    "90 bb8", "270 bb1", "270 bb6", "270 bb7", "270 bb8"))){
    which5 = c("gofKendallCvM", "gofKendallKS")
  }
  res = c("gofHybrid", which1, which2, which3, which4, which5)
  res
}

gofWhichCopula = function(test){
  if (is.element(test, c("gofRosenblattSnB", "gofRosenblattSnC", "gofADChisq", "gofADGamma", "gofSn", "gofRn"))){
    res = c("gaussian", "t", "clayton", "frank", "gumbel")#, "amh", "joe")
  }
  if (is.element(test, c("gofPIOSRn", "gofPIOSTn", "gofKernel"))){
    res = c("gaussian", "t", "clayton", "frank", "gumbel")
  }
  if (is.element(test, c("gofWhite"))){
    res = c(#"independence", 
            "gaussian", "t", "clayton", "frank", "gumbel")#, "joe", 
#            "survival clayton", "survival gumbel", "survival joe", "90 clayton", 
#            "90 gumbel", "90 joe", "270 clayton", "270 gumbel", "270 joe")
  }
  if (is.element(test, c("gofKendallCvM", "gofKendallKS"))){
    res = c(#"independence", 
            "gaussian", "clayton", "frank", "gumbel")#, "joe", 
#            "survival clayton", "survival gumbel", "survival joe", "90 clayton", 
#            "90 gumbel", "90 joe", "270 clayton", "270 gumbel", "270 joe", 
#            "bb1", "bb6", "bb7", "bb8", "survival bb1", "survival bb6", 
#            "survival bb7", "survival bb8", "90 bb1", "90 bb6", "90 bb7", 
#            "90 bb8", "270 bb1", "270 bb6", "270 bb7", "270 bb8")
  }
  if (is.element(test, c("gofHybrid"))){
    print("The available copula depend on the used tests.")
  }
  res
}