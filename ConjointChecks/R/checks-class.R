setClassUnion("list.null", c("list","NULL"))
setClass("checks", representation(N="matrix", n="matrix", Checks="list.null",tab="matrix",means="list",check.counts="matrix"))

