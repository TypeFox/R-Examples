## last modified June 2002

mixconstr <- function(conpi = "NONE", conmu = "NONE", consigma = "NONE", 
    fixpi = NULL, fixmu = NULL, fixsigma = NULL, cov = NULL, 
    size = NULL) 
{
    if (is.na(match(conpi, c("NONE", "PFX")))) 
        stop(paste("Unknown conpi ", conpi, ".", sep = ""))
    if (is.na(match(conmu, c("NONE", "MFX", "MEQ", "MES", "MGC")))) 
        stop(paste("Unknown conmu ", conmu, ".", sep = ""))
    if (is.na(match(consigma, c("NONE", "SFX", "SEQ", "FCV", 
        "CCV", "BINOM", "NBINOM", "POIS")))) 
        stop(paste("Unknown consigma ", consigma, ".", sep = ""))
    if (!(is.null(fixpi) | is.vector(fixpi))) 
        stop("fixpi is invalid.")
    if (conpi == "NONE" & !is.null(fixpi)) {
        fixpi <- NULL
        warning("fixpi was not used")
    }
    if (!(is.null(fixmu) | is.vector(fixmu))) 
        stop("fixmu is invalid.")
    if (conmu == "NONE" & !is.null(fixmu)) {
        fixmu <- NULL
        warning("fixmu was not used")
    }
    if (!(is.null(fixsigma) | is.vector(fixsigma))) 
        stop("fixsigma is invalid.")
    if (consigma == "NONE" & !is.null(fixsigma)) {
        fixsigma <- NULL
        warning("fixsigma was not used")
    }
    if (!(is.null(cov) | is.numeric(cov))) 
        stop("cov is invalid.")
    if (!(is.null(size) | is.numeric(size))) 
        stop("size is invalid.")
    if (!is.null(size) & consigma == "BINOM") 
        size <- ceiling(size)
    if (conpi == "PFX" & !is.vector(fixpi)) 
        stop("PFX needs fixpi.")
    else if (conmu == "MFX" & !is.vector(fixmu)) 
        stop("MFX needs fixmu.")
    else if (consigma == "SFX" & !is.vector(fixsigma)) 
        stop("SFX needs fixsigma.")
    else if (consigma == "FCV" & length(cov) != 1) 
        stop("FCV needs a scalar cov.")
    else if (consigma == "BINOM" & !is.vector(size)) 
        stop("BINOM needs size.")
    else if (consigma == "NBINOM" & !is.vector(size)) 
        stop("NBINOM needs size.")
    else if (conmu == "MEQ" & consigma == "SEQ") 
        stop("MEQ cannot have SEQ.")
    else if (conmu == "MEQ" & consigma == "FCV") 
        stop("MEQ cannot have FCV.")
    else if (conmu == "MEQ" & consigma == "CCV") 
        stop("MEQ cannot have CCV.")
    else if (conmu == "MEQ" & consigma == "POIS") 
        stop("MEQ cannot have POIS.")
    else if (conmu == "MFX" & consigma == "SFX" & length(fixmu) != 
        length(fixsigma)) 
        stop("Invalid constraints.")
    else if (conmu == "MFX" & consigma == "BINOM" & length(fixmu) != 
        length(size)) 
        stop("Invalid constraints.")
    else if (conmu == "MFX" & consigma == "NBINOM" & length(fixmu) != 
        length(size)) 
        stop("Invalid constraints.")
    else list(conpi = conpi, conmu = conmu, consigma = consigma, 
        fixpi = fixpi, fixmu = fixmu, fixsigma = fixsigma, cov = cov, 
        size = size)
}
