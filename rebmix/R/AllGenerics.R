setGeneric("RNGMIX",
  function(model = "RNGMIX",
    Dataset.name = character(),
    rseed = -1,
    n = numeric(),
    Theta = list(), ...)
  standardGeneric("RNGMIX"))

setGeneric("REBMIX",
  function(model = "REBMIX",
    Dataset = list(),
    Preprocessing = character(),
    cmax = 15,
    Criterion = "AIC",
    pdf = character(),
    theta1 = numeric(),
    theta2 = numeric(),
    K = numeric(),
    y0 = numeric(),
    ymin = numeric(),
    ymax = numeric(),
    ar = 0.1,
    Restraints = "loose", ...)
  standardGeneric("REBMIX"))
  
setGeneric("coef",
  function(object = NULL,
    pos = 1, ...)
  standardGeneric("coef"))
  
setGeneric(".IC", function(x = NULL, Criterion = "AIC", pos = 1, ...) standardGeneric(".IC"))
  
setGeneric("logL", function(x = NULL, pos = 1, ...) standardGeneric("logL")) 
setGeneric("AIC", function(x = NULL, pos = 1, ...) standardGeneric("AIC"))     
setGeneric("AIC3", function(x = NULL, pos = 1, ...) standardGeneric("AIC3"))    
setGeneric("AIC4", function(x = NULL, pos = 1, ...) standardGeneric("AIC4"))    
setGeneric("AICc", function(x = NULL, pos = 1, ...) standardGeneric("AICc"))    
setGeneric("BIC", function(x = NULL, pos = 1, ...) standardGeneric("BIC"))    
setGeneric("CAIC", function(x = NULL, pos = 1, ...) standardGeneric("CAIC"))    
setGeneric("HQC", function(x = NULL, pos = 1, ...) standardGeneric("HQC"))    
setGeneric("MDL2", function(x = NULL, pos = 1, ...) standardGeneric("MDL2"))    
setGeneric("MDL5", function(x = NULL, pos = 1, ...) standardGeneric("MDL5"))    
setGeneric("AWE", function(x = NULL, pos = 1, ...) standardGeneric("AWE"))    
setGeneric("CLC", function(x = NULL, pos = 1, ...) standardGeneric("CLC"))   
setGeneric("ICL", function(x = NULL, pos = 1, ...) standardGeneric("ICL"))   
setGeneric("ICLBIC", function(x = NULL, pos = 1, ...) standardGeneric("ICLBIC"))   
setGeneric("PRD", function(x = NULL, pos = 1, ...) standardGeneric("PRD"))   
setGeneric("SSE", function(x = NULL, pos = 1, ...) standardGeneric("SSE"))   
setGeneric("PC", function(x = NULL, pos = 1, ...) standardGeneric("PC"))

setGeneric("boot",
  function(x = NULL,
    pos = 1,
    Bootstrap = "parametric",
    B = 100, 
    n = numeric(), 
    replace = TRUE,
    prob = numeric(), ...)
  standardGeneric("boot"))
  
setGeneric("RCLRMIX",
  function(model = "RCLRMIX",
    x = NULL,
    pos = 1, 
    Zt = factor(), ...)
  standardGeneric("RCLRMIX"))  

setGeneric("RCLSMIX",
  function(model = "RCLSMIX",
    x = list(),
    Dataset = data.frame(), 
    Zt = factor(), ...)
  standardGeneric("RCLSMIX"))

