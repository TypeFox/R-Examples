#' Partition de.lorean object by cells
#'
#' @param dl de.lorean object
#' @param pieces How many pieces to partition into
#'
partition.de.lorean <- function(
    dl,
    pieces = 2
) {
    partition <- (1:(dim(dl)[2]) %% pieces) + 1
    cells <- sample(colnames(dl$expr))
    get.piece <- function(p) {
        partition.cells <- cells[partition == p]
        cell.filter <- function(cells) cells %in% partition.cells
        filter.cells(dl, cell.filter)
    }
    lapply(1:pieces, get.piece)
}

#' Test robustness of pseudotime estimation on subsets of de.lorean
#' object
#'
#' @param dl de.lorean object
#' @param pieces How many pieces to partition into
#'
test.robustness.de.lorean <- function(
    dl,
    pieces = 2
) {
    # Partition de.lorean object into several pieces
    partition <- partition.de.lorean(dl, pieces)
    # Define function to fit each piece
    run.model <- function(piece) {
        piece <- prepare.for.stan(piece)
        piece <- compile.model(piece)
        piece <- find.best.tau(piece)
        piece <- fit.model(piece)
        process.posterior(piece)
    }
    # Fit full de.lorean object
    full.model <- run.model(dl)
    # Get tau posterior
    full.tau <- (full.model$samples.l$tau
                    %>% dplyr::select(-c)
                    %>% mutate(fit="part"))
    # Fit each piece
    piece.models <- lapply(partition, run.model)
    # Get tau posterior from each piece and adjust mean to match full fit
    # mean
    get.tau.posterior <- function(piece.model) {
        posterior <- (piece.model$samples.l$tau
            %>% dplyr::select(-c)
            %>% mutate(fit="full"))
        cells <- unique(posterior$cell)
        full.mean <- mean((full.tau %>% filter(cell %in% cells))$tau)
        piece.mean <- mean(posterior$tau)
        posterior$tau <- posterior$tau - piece.mean + full.mean
        posterior
    }
    # Gather tau posteriors from each piece
    pieces.tau <- do.call(rbind,
                          lapply(piece.models, get.tau.posterior))
    # Combine fit pieces with fit full model
    tau.posterior <- rbind(full.tau, pieces.tau)
    # Sort by median tau
    cells <- (tau.posterior
        %>% group_by(cell)
        %>% dplyr::summarise(median.tau=median(tau))
        %>% arrange(median.tau)
    )
    # Create result list with plot
    list(tau.posterior = tau.posterior,
         gp = (ggplot(tau.posterior
                      %>% mutate(cell=factor(cell, levels=cells$cell)),
                      aes(x=cell, y=tau, fill=fit))
                  + geom_boxplot()
              ))
}


