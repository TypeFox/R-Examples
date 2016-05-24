################################################
#### ECOGEN CLASS DEFINITION
################################################

#' ecogen class
#' @name ecogen-class
#' @keywords internal
#' @slot XY P data frame
#' @slot P P data frame
#' @slot G G data frame
#' @slot A A data frame
#' @slot E E data frame
#' @slot S S data frame
#' @slot C C data frame
#' @slot OUT results
#' @slot INT int.gendata slot
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @aliases ecogen-class


setClass("ecogen",
         
         representation(XY = "data.frame",
                        P = "data.frame",
                        G = "data.frame",
                        A = "data.frame",
                        E = "data.frame",
                        S = "data.frame",
                        C = "data.frame",
                        OUT = "list",
                        INT = "int.gendata"),
         
         prototype(XY = data.frame(), 
                   P = data.frame(),
                   G = data.frame(),
                   A = data.frame(),
                   E = data.frame(),
                   S = data.frame(), 
                   C = data.frame(),
                   OUT = list(),
                   INT = new("int.gendata"))
         )
