setClass("wt.filter", representation(L="integer", level="integer",
                                     h="numeric", g="numeric",
                                     wt.class="character", wt.name="character",
                                     transform="character"))

setClass("dwt", representation(W="list", V="list", filter="wt.filter",
                               level="integer", n.boundary="numeric",
                               boundary="character", series="matrix",
                               class.X="character", attr.X="list",
                               aligned="logical", coe = "logical"))

setClass("modwt", representation(W="list", V="list", filter="wt.filter",
                                 level="integer", n.boundary="numeric",
                                 boundary="character", series="matrix",
                                 class.X="character", attr.X="list",
                                 aligned="logical", coe = "logical"))

setClass("mra", representation(D="list", S="list", filter="wt.filter",
                               level="integer", boundary="character",
                               series="matrix", class.X="character",
                               attr.X="list", method="character"))

#setMethod("plot", signature=(x="x", y="missing"), plot.dwt)

#setMethod("print", signature(x="dwt"), function(x,...) print.dwt(x,...))
#UseMethod("print", dwt)
