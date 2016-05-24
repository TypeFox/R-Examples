# Fix the Non-standard evaluation usage for check()
if(getRversion() >= "2.15.1"){
    utils::globalVariables(c("value", "variable", "glycan", "groups", ".",
                             "transpose", "gid", "mxxx", "s"))
}

#' Total Area Normalization of glycan data
#'
#' Returns glycans normalized with Total Area Normalization approach.
#'
#' @author Ivo Ugrina
#' @export tanorm
#' @param d data frame in long format containing glycan measurements
#' @param grouping should data be normalized per groups
#' @return Returns a data.frame with original glycan values substituted by normalized ones
#' @details
#' Input data frame should have at least the following three columns: \cr
#'   - gid - representing a unique name of a sample \cr
#'   - glycan - representing glycan names \cr
#'   - value - representing measured values \cr
#' and if the grouping argument is \code{TRUE} it should also have column: \cr
#'   - groups - representing groupings (e.g. IgG1, IgG2 and IgG4)
#' @examples
#' data(mpiu)
#' mpiun <- tanorm(mpiu)
#' head(mpiun)
tanorm <- function(d, grouping=FALSE){
    warning("Version 0.3 of glycanr introduces a change in tanorm function. It expects data frame input in long format now.",
            call. = FALSE)

    if(grouping==FALSE){
        return(tanorm_basic(d))
    }else{
        return(tanorm_groups(d))
    }
}

tanorm_basic <- function(d){
    d <- d %>% 
        dplyr::group_by(gid) %>%
		dplyr::mutate(value = value/sum(value, na.rm = TRUE)*100) %>%
		dplyr::ungroup()
	d
}

tanorm_groups <- function(d){	
	d <- d %>%
        dplyr::group_by(groups, gid) %>%
		dplyr::mutate(value = value/sum(value, na.rm = TRUE)*100) %>%
		dplyr::ungroup()
	d
}

#' Reference Peak Normalization of glycan data
#'
#' Returns glycans normalized with Reference Peak Normalization approach.
#'
#' @author Ivo Ugrina, Lucija Klarić
#' @export refpeaknorm
#' @param d data frame in long format containing glycan measurements
#' @param grouping should data be normalized per groups
#' @param peak glycan name to use as the reference peak. If \code{NULL}
#'   peak with maximal value (summed through all samples) is used
#' @return Returns a data.frame with original glycan values substituted by normalized ones
#' @details
#' Input data frame should have at least the following three columns: \cr
#'   - gid - representing a unique name of a sample \cr
#'   - glycan - representing glycan names \cr
#'   - value - representing measured values \cr
#' and if the grouping argument is \code{TRUE} it should also have column: \cr
#'   - groups - representing groupings (e.g. IgG1, IgG2 and IgG4)
#' @examples
#' data(mpiu)
#' mpiun <- refpeaknorm(mpiu)
#' head(mpiun)
refpeaknorm <- function(d, grouping=FALSE, peak=NULL){
    if(grouping==FALSE){
        return(refpeaknorm_basic(d, peak))
    }else{
        return(refpeaknorm_groups(d, peak))
    }
}

refpeaknorm_basic <- function(d, peak=NULL){
    if(is.null(peak)){
        tmp <- d %>% 
            dplyr::group_by(glycan) %>% 
            dplyr::summarise(s = sum(value, na.rm=TRUE)) %>% 
            dplyr::ungroup() %>% 
            dplyr::arrange(desc(s))

        ref_peak <- as.character(tmp[1,]$glycan)
    } else {
        ref_peak <- peak
    }

    d <- d %>%
        dplyr::group_by(gid) %>% 
        dplyr::mutate(value=value/value[glycan==ref_peak]) %>% 
        dplyr::ungroup()
	return(d)
}

refpeaknorm_groups <- function(d, peak=NULL){
    d <- d %>% 
        dplyr::group_by(groups) %>% 
        dplyr::do(refpeaknorm_basic(., peak)) %>% 
        dplyr::ungroup()

    d
}

#' Median Normalization of glycan data
#'
#' Returns glycans normalized with Median Normalization approach.
#'
#' @author Ivo Ugrina, Lucija Klarić
#' @export mediannorm
#' @param d data frame in long format containing glycan measurements
#' @param grouping should data be normalized per groups
#' @return Returns a data.frame with original glycan values substituted by normalized ones
#' @details
#' Input data frame should have at least the following three columns: \cr
#'   - gid - representing a unique name of a sample \cr
#'   - glycan - representing glycan names \cr
#'   - value - representing measured values \cr
#' and if the grouping argument is \code{TRUE} it should also have column: \cr
#'   - groups - representing groupings (e.g. IgG1, IgG2 and IgG4)
#' @examples
#' data(mpiu)
#' mpiun <- mediannorm(mpiu)
#' head(mpiun)
mediannorm <- function(d, grouping=FALSE){
    if(grouping==FALSE){
        return(mediannorm_basic(d))
    }else{
        return(mediannorm_groups(d))
    }
}

mediannorm_basic <- function(d) {
	d <- d %>%
		dplyr::group_by(gid) %>%
		dplyr::mutate(value = (value - stats::median(value, na.rm = TRUE))/IQR(value, na.rm = TRUE)) %>%
		dplyr::ungroup() 
    d
}

mediannorm_groups = function(d) {
	d <- d %>%
		dplyr::group_by(groups, gid) %>%
		dplyr::mutate(value = (value - stats::median(value, na.rm = TRUE))/IQR(value, na.rm = TRUE)) %>%
		dplyr::ungroup() 
    d
}

#' Median Quotient Normalization of glycan data
#'
#' Returns glycans normalized with Median Quotient Normalization approach.
#'
#' @author Ivo Ugrina, Lucija Klarić
#' @export medianquotientnorm
#' @param d data frame in long format containing glycan measurements
#' @param grouping should data be normalized per groups
#' @return Returns a data.frame with original glycan values substituted by normalized ones
#' @details
#' Input data frame should have at least the following three columns: \cr
#'   - gid - representing a unique name of a sample \cr
#'   - glycan - representing glycan names \cr
#'   - value - representing measured values \cr
#' and if the grouping argument is \code{TRUE} it should also have column: \cr
#'   - groups - representing groupings (e.g. IgG1, IgG2 and IgG4)
#' @references
#' Dieterle F,Ross A, Schlotterbeck G, Senn H.: \cr
#' Probabilistic Quotient Normalization as Robust Method to Account for
#' Diluition of Complex Biological Mixtures. Application in 1H NMR Metabolomics. \cr
#' Anal Chem 2006;78:4281-90. \cr
#' \url{http://dx.doi.org/10.1021/ac051632c}
#' @examples
#' data(mpiu)
#' mpiun <- medianquotientnorm(mpiu)
#' head(mpiun)
medianquotientnorm <- function(d, grouping=FALSE){
    if(grouping==FALSE){
        return(medianquotientnorm_basic(d))
    }else{
        return(medianquotientnorm_groups(d))
    }
}

medianquotientnorm_basic <- function(d){
    ref_chromx <- d %>% 
        dplyr::group_by(glycan) %>% 
        dplyr::mutate(mxxx=value/stats::median(value, na.rm=TRUE)) %>% 
        dplyr::ungroup()

    d <- ref_chromx %>% 
        dplyr::group_by(gid) %>% 
        dplyr::mutate(value=value/stats::median(mxxx, na.rm=TRUE)) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(-mxxx)

    d
}

medianquotientnorm_groups <- function(d){
    d <- d %>% 
        dplyr::group_by(groups) %>% 
        dplyr::do(medianquotientnorm_basic(.)) %>% 
        dplyr::ungroup()

    d
}

#' Quantile Normalization of glycan data
#'
#' Returns glycans normalized with Quantile Normalization approach.
#'
#' @author Ivo Ugrina, Lucija Klarić
#' @export quantilenorm
#' @param d data frame in long format containing glycan measurements
#' @param transpose transpose the data prior to normalization
#' @param grouping should data be normalized per groups
#' @return Returns a data.frame with original glycan values substituted by normalized ones
#' @details
#' Input data frame should have at least the following three columns: \cr
#'   - gid - representing a unique name of a sample \cr
#'   - glycan - representing glycan names \cr
#'   - value - representing measured values \cr
#' and if the grouping argument is \code{TRUE} it should also have column: \cr
#'   - groups - representing groupings (e.g. IgG1, IgG2 and IgG4)
#' @references
#' Bolstad, B. M., Irizarry R. A., Astrand, M, and Speed, T. P.: \cr
#' A Comparison of Normalization Methods for High Density Oligonucleotide
#' Array Data Based on Bias and Variance.\cr
#' Bioinformatics 19(2), p. 185-193, 2003. \cr
#' \url{http://dx.doi.org/10.1093/bioinformatics/19.2.185}
#' @examples
#' data(mpiu)
#' mpiun <- quantilenorm(mpiu)
#' head(mpiun)
#' 
#' # transpose (change) subjects and measurements
#' mpiunt <- quantilenorm(mpiu, transpose=TRUE)
#' head(mpiunt)
quantilenorm <- function(d, grouping=FALSE, transpose=FALSE){
    if(!requireNamespace("preprocessCore", quietly=TRUE)){
        stop("Unable to proceed since package preprocessCore from
        BioConductor is not available on this system. This
        package is a prerequisite to use the quantilenorm function!")
    }

    if(grouping==FALSE){
        return(quantilenorm_basic(d, transpose))
    }else{
        return(quantilenorm_groups(d, transpose))
    }
}


quantilenorm_basic <- function(d, transpose=FALSE){
    if(!requireNamespace("preprocessCore", quietly=TRUE)){
        stop("Unable to proceed since package preprocessCore from
        BioConductor is not available on this system. This
        package is a prerequisite to use the quantilenorm function!")
    }

    tmp <- d %>% 
        dplyr::select(gid, glycan, value) %>% 
        tidyr::spread(glycan, value)

	glycans <- as.character(unique(d$glycan))

    if(transpose){
        tempD <- preprocessCore::normalize.quantiles(as.matrix(tmp[, glycans]))
    }else{
        tempD <- preprocessCore::normalize.quantiles(t(as.matrix(tmp[, glycans])))
    }

    if(transpose){
        tmp[,glycans] <- tempD
    }else{
        tmp[,glycans] <- t(tempD)
    }

    tmp <- tmp %>% 
        tidyr::gather_("glycan", "value", glycans)

    d <- d %>% 
        dplyr::select(-value)
    d <- merge(d, tmp, by=c("gid", "glycan"))

	return(d)
}

quantilenorm_groups <- function(d, transpose=FALSE){
    d <- d %>% 
        dplyr::group_by(groups) %>% 
        dplyr::do(quantilenorm_basic(., transpose)) %>% 
        dplyr::ungroup()

    d
}

