makeFrontier <-
function(dataset, treatment, outcome, match.on, 
         keep.vars = NULL, QOI = 'FSATT', metric = 'Mahal',
         ratio = 'fixed', breaks = NULL){

    # Check the frontier arguments 
    checkArgs(QOI, metric, ratio)
    
    # Check data and trim to suff we need
    dataset <- checkDat(dataset, treatment, outcome, match.on, keep.vars)
    
    if(QOI == 'FSATT' & metric == 'Mahal' & ratio == 'variable'){
        frontier <- MahalFrontierFSATT(treatment = treatment,
                                       outcome = outcome,
                                       dataset = dataset,
                                       ratio = ratio,
                                       match.on = match.on)
        class(frontier) <- "MahalFSATTClass"
        return(frontier)        
    }
    if(QOI == 'FSATT' & metric == 'Mahal' & ratio == 'fixed'){
        frontier <- MahalFrontierFSATT(treatment = treatment,
                                       outcome = outcome,
                                       dataset = dataset,
                                       ratio = ratio,
                                       match.on = match.on)
        class(frontier) <- "MahalFSATTClass"
        return(frontier)
    }
    if(QOI == 'SATT' & metric == 'L1' & ratio == 'fixed'){
        frontier <- L1FrontierSATT(treatment = treatment,
                                   outcome = outcome,
                                   dataset = dataset,
                                   breaks = breaks,
                                   match.on = match.on)
        class(frontier) <- "L1SATTClass"
        return(frontier)
    }
    if(QOI == 'FSATT' & metric == 'L1' & ratio == 'variable'){
        print("See the 'cem' package")
        #return(L1FrontierCEM(treatment = treatment, outcome = outcome, dataset = dataset, breaks = breaks))
    }
    else{
        msg <- paste('the ', ratio, '-ratio ', metric, ' theoretical frontier is not presently calculable.', sep = '')
        customStop(msg, 'makeFrontier()')
    }
}
