checkDat <-
function(dataset, treatment, outcome, match.on, keep.vars){
    keep.columns <- unique(c(treatment, outcome, match.on, keep.vars))
   
    # Check if all the variables are in the data
    if(sum(!(keep.columns %in% colnames(dataset))) > 0){ 
        missing.cols <- keep.columns[!(keep.columns %in% colnames(dataset))]
        error.msg <- paste('the following columns are not in the data: ',
                           paste(missing.cols, collapse = '\n'), sep = '\n'
                           )
        customStop(error.msg, 'makeFrontier()')
    }
    
    # Make sure user isn't trying to match on the treatment or the outcome
    if(treatment %in% match.on){
        customStop("the treatment is in 'match.on'. You shouldn't match on the treatment, that's bad.", 'makeFrontier()')
    }
    if(outcome %in% match.on){
        customStop("the outcome is in 'match.on'. You shouldn't match on the treatment, that's bad.", 'makeFrontier()')
    }

    # Check treatment
    if(sum(!(dataset[,treatment] %in% c(0,1))) != 0){
        customStop('the treatment must be either 0/1 (integers) or "TRUE"/"FALSE" (logical).', 'makeFrontier()')
    }
    
    # Trim the dataset to the stuff we need
    dataset <- dataset[, keep.columns]

    # Check for missing values
    if(sum(is.na(dataset)) != 0){
        customStop("missing values in the data; remove them (or impute) and try again.", 'makeFrontier()')
    }
    
    rownames(dataset) <- 1:nrow(dataset)
    return(dataset)

}
