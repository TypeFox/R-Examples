regionimport <- function(dtblock, region_vars, name_batch){
    noncommon <- region_vars[region_vars%nin%names(dtblock)]
    if(length(noncommon)>=1){
    stop("Following variables not found '", 
        paste(noncommon, collapse=", ") , "'")
    } 
    region_block <- lapply(region_vars, function(x) dtblock[[x]][])
    region_data <- as.data.frame(t(as.data.frame(region_block)))
    region_data <- region_data[-1,]
    new_frame <- data.frame(region_data )
    names(new_frame) <- c("rid","units")
    new_frame$analyte <- rownames(region_data)
    new_frame$batch <- name_batch
    region_data <- new_frame
    region_data$batch_analyte <- apply(region_data[,c("batch","analyte")],1,
                                        function(x) paste(x, collapse="*") )
    vars <- c("batch_analyte","batch","analyte","rid","units")
    region_data <- region_data[,vars]
    region_data
}



