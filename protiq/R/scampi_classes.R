## class for input data, mostly used internally
scampi <- setClass("scampi",
                   representation(peptides = "data.frame",
                                  proteins = "data.frame",
                                  edgespp = "data.frame"))
scampi <- function(peptides, proteins, edgespp) {
        new("scampi", peptides=peptides, proteins=proteins, edgespp=edgespp)
}

## class for returning results (quantified proteins/peptides)
scampiVal <- setClass("scampiVal",
                      representation(call = "call",
                                     parameters = "list",
                                     ppGraph = "graphNEL",
                                     ccList = "list"),
                      contains = c("scampi"))
scampiVal <- function(peptides, proteins, edgespp, call, parameters,
                      ppGraph, ccList) {
  new("scampiVal", peptides=peptides, proteins=proteins, edgespp=edgespp,
      call=call, parameters=parameters, ppGraph=ppGraph, ccList=ccList)
}
