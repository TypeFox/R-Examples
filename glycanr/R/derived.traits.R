#' Derived traits for Glycan peaks in IgG for UPLC
#'
#' Calcuates values of derived traits for Glycan peaks in IgG for UPLC
#'
#' @author Ivo Ugrina, Frano Vučković
#' @export iudt
#' @param data data frame that holds columns representing Glycans.
#'   These column names should start with 'GP'.
#' @param method year of the derived traits definition. By default 2014.
#' @param print.exp.names If \code{TRUE} return expected column names 
#'   representing glycans. 
#' @return Returns the data frame with derived traits
#' @details
#' Calculates derived traits from basic glycan peaks. User can choose
#' which definition of the derived traits he will use
#' (see references for different versions/definitions of derived traits).
#' 
#' @references
#' Jennifer E. Huffman et al. (2014)
#' "Comparative Performance of Four Methods for High-throughput Glycosylation Analysis of Immunoglobulin G in Genetic and Epidemiological Research*"
#' \url{http://dx.doi.org/10.1074/mcp.M113.037465}
iudt <- function(data=NULL, method="2014", print.exp.names=FALSE) {
    x <- NULL

    if(method == "2014"){
        if(!print.exp.names & !is.data.frame(data)){
            warning("Either use print.exp.names=TRUE or
                    set the paramater data to be a data frame")
            return(x)
        }
        x <- igg.uplc.derived.traits.2014(data, print.exp.names)
    }

    x
}


#' Derived traits for Glycan peaks in PLASMA for HPLC
#'
#' Calcuates values of derived traits for Glycan peaks in Plasma for HPLC
#'
#' @author Ivo Ugrina, Lucija Klarić
#' @export phdt
#' @param data data frame that holds columns representing Glycans.
#'   These column names should start with 'GP'.
#' @param method year of the derived traits definition. By default 2011.
#' @param print.exp.names If \code{TRUE} return expected column names 
#'   representing glycans. 
#' @return Returns the data frame with derived traits
#' @details
#' Calculates derived traits from basic glycan peaks. User can choose
#' which definition of the derived traits he will use
#' (see references for different versions/definitions of derived traits).
#' 
#' @references
#' Lu et al. (2011)
#' "Screening Novel Biomarkers for Metabolic Syndrome by Profiling 
#'  Human Plasma N-Glycans in Chinese Han and Croatian Populations"
#' \url{http://dx.doi.org/10.1021/pr2004067}
phdt <- function(data=NULL, method="2011", print.exp.names=FALSE) {
    x <- NULL

    if(method == "2011"){
        if(!print.exp.names & !is.data.frame(data)){
            warning("Either use print.exp.names=TRUE or
                    set the paramater data to be a data frame")
            return(x)
        }
        x <- plasma.hplc.derived.traits.2011(data, print.exp.names)
    }

    x
}


#' Derived traits for Glycan peaks in IgG for LCMS
#'
#' Calcuates values of derived traits for Glycan peaks in IgG for LCMS
#'
#' @author Ivo Ugrina
#' @export ildt
#' @param data data frame that holds columns representing Glycans.
#' @param method year of the derived traits definition. By default 2014.
#' @param print.exp.names If \code{TRUE} return expected column names 
#'   representing glycans. 
#' @return Returns the data frame with derived traits
#' @details
#' Calculates derived traits from basic glycan peaks. User can choose
#' which definition of the derived traits he will use
#' (see references for different versions/definitions of derived traits).
#' 
#' @references
#' Jennifer E. Huffman et al. (2014)
#' "Comparative Performance of Four Methods for High-throughput Glycosylation Analysis of Immunoglobulin G in Genetic and Epidemiological Research*"
#' \url{http://dx.doi.org/10.1074/mcp.M113.037465}
ildt <- function(data=NULL, method="2014", print.exp.names=FALSE) {
    x <- NULL
 
    if(method == "2014"){
        if(!print.exp.names & !is.data.frame(data)){
            warning("Either use print.exp.names=TRUE or
                    set the paramater data to be a data frame")
            return(x)
        }
        x <- igg.lcms.derived.traits.2014(data, print.exp.names)
    }

    x
}

#' Translate names between computer readable and human readable
#' for derived traits of IgG with LCMS
#'
#' Translates names between computer readable and human readable
#' for derived traits of IgG with LCMS
#'
#' @author Ivo Ugrina
#' @export ildt.translate
#' @importFrom utils read.delim
#' @param orignames vector; type string
#' @param method year of the derived traits definition. By default 2014.
#' @param to type of translation. If \code{inverse} is used everything will be
#'   translated. For \code{computer} names will be translated to computer
#'   readable, and for \code{human} names will be translated to human readable.
#' @return Returns a character vector with original and translated names
#' @details
#' User can choose which definition of the derived traits he will use
#' (see references for different versions/definitions of derived traits).
#' 
#' @references
#' Jennifer E. Huffman et al. (2014)
#' "Comparative Performance of Four Methods for High-throughput Glycosylation Analysis of Immunoglobulin G in Genetic and Epidemiological Research*"
#' \url{http://dx.doi.org/10.1074/mcp.M113.037465}
ildt.translate <- function(orignames, to="inverse", method="2014") {
    x <- NULL
  
    if(method == "2014"){
       x <- ildt.translate.2014(orignames, to)
    }

    x
}

#' Translate names between computer readable and human readable
#' for derived traits of IgG with UPLC
#'
#' Translates names between computer readable and human readable
#' for derived traits of IgG with UPLC
#'
#' @author Ivo Ugrina
#' @export iudt.translate
#' @importFrom utils read.delim
#' @param orignames vector; type string
#' @param method year of the derived traits definition. By default 2014.
#' @param to type of translation. If \code{inverse} is used everything will be
#'   translated. For \code{computer} names will be translated to computer
#'   readable, and for \code{human} names will be translated to human readable.
#' @return Returns a character vector with original and translated names
#' @details
#' User can choose which definition of the derived traits he will use
#' (see references for different versions/definitions of derived traits).
#' 
#' @references
#' Jennifer E. Huffman et al. (2014)
#' "Comparative Performance of Four Methods for High-throughput Glycosylation Analysis of Immunoglobulin G in Genetic and Epidemiological Research*"
#' \url{http://dx.doi.org/10.1074/mcp.M113.037465}
iudt.translate <- function(orignames, to="inverse", method="2014") {
    x <- NULL
  
    if(method == "2014"){
       x <- iudt.translate.2014(orignames, to)
    }

    x
}
