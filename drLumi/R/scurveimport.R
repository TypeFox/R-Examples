scurveimport <- function(dtblock, scurve_vars,name_batch){
    noncommon <- scurve_vars[scurve_vars%nin%names(dtblock)]
    if(length(noncommon)>=1){
        stop("Following variables not found '",paste(noncommon, 
                collapse=", "), "'")
    } 
    scurves_block <- lapply(scurve_vars, function(x) dtblock[[x]][])
    names(scurves_block) <- scurve_vars
    scurves_block <- ldply(scurves_block)
    type_var <- apply(scurves_block[,c(1,2)], 1, 
                function(x) paste(x, collapse="_"))
    scurves_block[,1] <- NULL
    scurves_block[,1] <- NULL
    scurves_block <- as.data.frame(t(scurves_block))
    names(scurves_block) <- type_var
    scurves_block$analyte <- rownames(scurves_block)
    scurves_block$batch <- name_batch
    scurves_data <- scurves_block
    scurves_data$batch_analyte <- apply(scurves_data[,c("batch","analyte")],1,
                function(x) paste(x, collapse="*") )
    names(scurves_data) <- gsub(" ","_",(tolower(names(scurves_data))))
    naux <- names(scurves_data)
    selec_names <- naux[naux%nin%c("batch_analyte")]
    vars <- c("batch_analyte",selec_names)
    scurves_data <- scurves_data[, vars]
    scurves_data
}
