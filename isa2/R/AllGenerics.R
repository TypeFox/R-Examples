
setGeneric("isa", function(data, ...) standardGeneric("isa"))
setGeneric("isa.normalize",
           function(data, ...) standardGeneric("isa.normalize"))
setGeneric("isa.iterate",
           function(normed.data, ...) standardGeneric("isa.iterate"))
setGeneric("isa.unique",
           function(normed.data, isaresult, ...) standardGeneric("isa.unique"))

setGeneric("ppa", function(data, ...) standardGeneric("ppa"))
setGeneric("ppa.normalize",
           function(data, ...) standardGeneric("ppa.normalize"))
setGeneric("ppa.iterate",
           function(normed.data, ...) standardGeneric("ppa.iterate"))
setGeneric("ppa.unique",
           function(normed.data, pparesult, ...) standardGeneric("ppa.unique"))

setGeneric("robustness",
           function(normed.data, ...) standardGeneric("robustness"))
setGeneric("isa.filter.robust",
           function(data, ...) standardGeneric("isa.filter.robust"))
setGeneric("ppa.filter.robust",
           function(data, ...) standardGeneric("ppa.filter.robust"))

setGeneric("isa.sweep",
           function(data, ...) standardGeneric("isa.sweep"))
setGeneric("sweep.graph",
           function(sweep.result, ...) standardGeneric("sweep.graph"))

setGeneric("plotModules",
           function(modules, ...) standardGeneric("plotModules"))

