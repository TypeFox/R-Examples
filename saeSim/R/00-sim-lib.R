################################################################################

setClass("sim_fun",
         slots = c(order = "numeric", call = "call"),
         contains = "function")

################################################################################

setClass("sim_setup", 
         slots = c(base = "data.frame", simName = "character"), 
         contains = "list")

setClass("summary.sim_setup",
         slots = c(
           sim_setup = "sim_setup",
           duration = "table",
           expression = "expression",
           dim = "numeric"))