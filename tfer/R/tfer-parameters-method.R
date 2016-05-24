setGeneric("parameters",
           function(object)
           standardGeneric("parameters")
           )


setMethod("parameters", signature(object = "transfer"),
          function(object){
              paras = para(object)
              cat("\n", "Simulation N = ", paras[1],"\n",
                  "Distance = ", paras[2],"\n",
                  "Average transferred = ", paras[4], "\n",
                  "% High persistence = ", paras[5], "\n",
                  paras[6], "<= % loss in first hour <=", paras[7],"\n",
                  paras[8], "<= % high persistence loss in first hour <=", paras[9],"\n",
                  paras[10], "<= % loss in j'th hour <=", paras[11],"\n",
                  paras[12], "<= % high persistence loss in j'th hour <=", paras[13],"\n",
                  paras[14], "<= % detected in lab <=", paras[15],"\n",
                  "Time = ", paras[16],"\n", "\n")
          }
          )
