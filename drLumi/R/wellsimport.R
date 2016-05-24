#' @importFrom reshape melt
wellsimport <- function(dtblock, well_vars, name_batch){
    noncommon <- well_vars[well_vars%nin%names(dtblock)]
    if(length(noncommon) >= 1){ 
        stop(paste("Following variables not found '",noncommon, "'"))
    }  
    well_block <- lapply(well_vars, function(x) dtblock[[x]][])
    well_block <- lapply(well_block, 
                    function(x) melt(x, id.vars = c("Location","Sample")))
    for(i in 1:length(well_vars)){  
        names(well_block[[i]]) <- c("well","sample","analyte",well_vars[i])
    }    
    well_data <- well_block[[1]]
    if(length(well_vars)>1){
        for(i in 2:length(well_vars)){
            well_data <- merge(well_data, well_block[[i]], 
                by=c("well","sample","analyte"), 
                all.x=TRUE, all.y=TRUE, sort=FALSE)
        }
    }
    
    well_aux <- unlist(lapply(as.character(well_data$well), 
                        function(x) strsplit(x,",")[[1]][2]))
    
    well_aux <- unlist(lapply(well_aux, function(x) substr(x,1,nchar(x)-1)))
    
    plate_aux <- unlist(lapply(as.character(well_data$well), 
                        function(x) strsplit(x,",")[[1]][1]))
    plate_aux <- unlist(lapply(as.character(plate_aux), 
                        function(x)strsplit(x,"(",fixed=TRUE)[[1]][2]))
    
    well <- paste("P",plate_aux,"_",well_aux, sep="")
    
    well_data$well <- well
    
    well_data$batch <- name_batch
    vars <- c("batch","well","analyte")
    well_data$batch_well_analyte <- apply(well_data[,vars],1,function(x){
        paste(x, collapse="*") })
    vars <- c("batch_well_analyte",vars,"sample",well_vars)
    well_data <- well_data[,vars]

    names(well_data) <- gsub(" ","_",(tolower(names(well_data))))
    names(well_data) <- gsub("%","pc",names(well_data))
    return(well_data)
}

