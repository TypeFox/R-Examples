
####  class for decoupled grids ####

setClass( "doubleEntry", 
          representation( meta = "list",
                          scale = "list",
                          coupled = "logical",
                          elements = "list",
                          constructs = "list",
                          elicitation = "list",
                          ratings = "array",
                          calcs = "list",
                          plotdata = "data.frame"))

