
### a class for model environments

setClass("ModelEnv",
    representation(
        env = "environment",
        get = "function",
        set = "function",
        hooks = "list"))

### a class for formulae

setClass("FormulaParts",
    representation(
        formula = "list"
    )
)

### model environments given by formulae       

setClass("ModelEnvFormula", contains = c("ModelEnv", "FormulaParts"))

### A prototype for a model class in R

setClass("StatModelCapabilities",
    representation(
        weights = "logical",
        subset  = "logical"),
    prototype(weights = TRUE, subset = TRUE)
)

setClass("StatModel",
    representation(
        name         = "character",
        dpp          = "function",
        fit          = "function",
        predict      = "function",
        capabilities = "StatModelCapabilities")
)

