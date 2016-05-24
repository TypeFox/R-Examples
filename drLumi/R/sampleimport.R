sampleimport <- function(dtblock, sample_vars, name_batch, 
                backname, stanname, samplename, 
                dilutionname, batchname){
    noncommon <- sample_vars[sample_vars%nin%names(dtblock)]
    if(length(noncommon)>=1){
        stop("Following variables not found '", 
        paste(noncommon, collapse=", ") ,"'")
    } 
    sample_block <- lapply(sample_vars, function(x) dtblock[[x]][])[[1]]
    back <- unique(grep(backname,sample_block[,samplename], value=TRUE))
    standard <- unique(grep(stanname,sample_block[,samplename], value=TRUE))
    sample_selection <- sample_block[,samplename]%in%c(back, standard)
    sample_block <- sample_block[sample_selection, c(samplename,dilutionname)]
    sample_block <- unique(sample_block)
    sample_data <- sample_block
    names(sample_data) <- c("sample","dilution")
    if (nrow(sample_data)>0) { 
        sample_data$batch <- name_batch
    } else {
        sample_data <- data.frame( sample = numeric(0), 
                        dilution = numeric(0), 
                        batch = character(0))
    }
    sample_data$batch_sample <- apply(sample_data[,c(batchname,"sample")],1,
                                    function(x) paste(x, collapse="*") )
    sample_data$batch <- sample_data[, batchname]
    vars <- c("batch_sample","batch","sample","dilution")
    sample_data <- sample_data[,vars]
    sample_data
}

