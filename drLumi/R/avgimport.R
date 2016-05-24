#' @importFrom reshape melt
avgimport <- function(dtblock, average_vars, name_batch){
    noncommon <- average_vars[average_vars%nin%names(dtblock)]
    if(length(noncommon)>=1){
        stop("Following variables not found '",
            paste(noncommon, collapse=", "), "'")
    } 
    avg_block <- lapply(average_vars, function(x) dtblock[[x]][])
    var_names <- lapply(avg_block, 
                    function(x) names(x) <- c("Sample",names(x)[-1]))
    for(i in 1:length(average_vars)){      
        names(avg_block[[i]]) <- var_names[[i]]    
    }
    avg_block <- lapply(avg_block, function(x) melt(x, id.vars = c("Sample")))
    for(i in 1:length(average_vars)){  
        names(avg_block[[i]]) <- c("sample","analyte",average_vars[i])
    }
    avg_data <- avg_block[[1]]
    if(length(average_vars)>1){
        for(i in 2:length(average_vars)){
            avg_data <- merge(avg_data, avg_block[[i]], 
                        by=c("sample","analyte"), 
                        all.x=TRUE, all.y=TRUE, sort=FALSE)
        }
    }
    avg_data$batch <- name_batch
    vars <- c("batch","sample","analyte")
    avg_data$batch_sample_analyte <- apply(avg_data[,vars],1,
                                    function(x) paste(x, collapse="*") )
    vars <- c("batch_sample_analyte","batch", "sample", 
                "analyte", average_vars)
    avg_data <- avg_data[,vars]
    names(avg_data) <- gsub(" ","_",(tolower(names(avg_data))))
    names(avg_data) <- gsub("%","pc_",names(avg_data))
    if( "standard_expected_concentration"%in%names(avg_data) ){
        old.name <- "standard_expected_concentration"
        new.name <- "st_exp_conc"
        avg_data <- rename.vars(avg_data, from=old.name, 
                        to=new.name, info=FALSE)
    }
    return(avg_data)
}

