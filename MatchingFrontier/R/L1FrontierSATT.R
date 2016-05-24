L1FrontierSATT <-
function(treatment, outcome, dataset, breaks, match.on){
    cat("Calculating L1 binnings...\n")
    binnings <- getBins(dataset, treatment, match.on, breaks)
    cat("Calculating L1 frontier... This may take a few minutes...\n")
    frontier <- binsToFrontier(binnings)
    out <- list(
        frontier = frontier,
        cuts = binnings$cuts,
        treatment = treatment,
        outcome = outcome,
        QOI = 'SATT',
        metric = 'L1',
        ratio = 'fixed',
        dataset = dataset,
        match.on = match.on
        )
    return(out)
}
