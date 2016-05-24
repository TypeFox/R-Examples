#' Calculate the roughness of the vector. The roughness is the RMS
#' of the differences between consecutive points.
#'
#' @param x Values
#'
calc.roughness <- function(x) {
    N <- length(x)
    stopifnot(N > 0)
    if (0 == N) return(0)
    S <- sd(x)
    if (S == 0) return(0)
    sqrt(sum((x[1:(N-1)] - x[2:N])**2) / (N-1)) / S
}


#' Permute a data frame, x. If group.col is given it should name an ordered
#' factor that the order of the permutation should respect.
#'
#' @param .df Data frame
#' @param group.col Name of an ordered factor that the permutation
#'    should respect.
#'
permute.df <- function(.df, group.col=NULL) {
    if (is.null(group.col)) {
        return(.df[sample(nrow(.df)),])
    } else {
        stopifnot(is.ordered(.df[,group.col]))
        return(
            .df
            %>% group_by_(as.symbol(group.col))
            %>% do(permute.df(.)))
    }
}


#' Calculate the roughness of the held out genes given the sample.
#'
#' @param dl de.lorean object
#' @param expr.held.out The expression matrix including the held out genes
#' @param sample.iter Which sample to use
#'
roughness.of.sample <- function(
    dl,
    expr.held.out=dl$expr.held.out,
    sample.iter=dl$best.sample)
with(dl, {
    tau <- tau.for.sample(dl, sample.iter=sample.iter)
    mean(apply(expr.held.out[,order(tau)], 1, calc.roughness))
})


#' Permute cells and test roughness of expression.
#'
#' @param dl de.lorean object
#' @param expr.held.out The expression matrix of the held out genes
#'
permuted.roughness <- function(
    dl,
    expr.held.out=dl$expr.held.out)
{
    permuted <- permute.df(dl$cell.map, "capture")
    apply(expr.held.out[,permuted$c], 1, calc.roughness)
}


#' Apply permutation based roughness test to held out genes
#'
#' @param dl de.lorean object
#' @param expr.held.out The expression matrix including the held out genes
#' @param num.perms Number of permutations to test
#'
roughness.of.permutations <- function(
    dl,
    expr.held.out=dl$expr.held.out,
    num.perms=1000)
{
    colMeans(sapply(1:num.perms, function(i) permuted.roughness(dl, expr.held.out)))
}


#' Calculate roughnesses under fit samples and also under random
#' permutations
#'
#' @param dl de.lorean object
#' @param expr.held.out The expression matrix including the held out genes
#' @param num.perms Number of permutations to test
#'
#' @export
#'
roughness.test <- function(
    dl,
    expr.held.out=dl$expr.held.out,
    num.perms=1000)
within(dl, {
    message('Performing roughness test on ', nrow(expr.held.out),
            ' genes and ', ncol(expr.held.out),
            ' cells. Using ', num.perms, ' permutations.')
    # Combine both types into a dataframe
    roughnesses <- rbind(
        data.frame(type="pseudotime", sample.iter=sample.iters(dl))
        %>% mutate(roughness=Vectorize(
                                roughness.of.sample,
                                vectorize.args='sample.iter'
                             )(
                                dl,
                                expr.held.out=expr.held.out,
                                sample.iter
                             )),
        data.frame(type="permutation",
                   sample.iter=NA,
                   roughness=roughness.of.permutations(dl, held.out.expr)))
    roughness.test <- t.test(
        x=filter(roughnesses, "permutation" == type)$roughness,
        y=filter(roughnesses, "pseudotime"  == type)$roughness,
        alternative="greater")
})


#' Plot results of roughness test
#'
#' @param dl de.lorean object
#'
#' @export
#'
roughnesses.plot <- function(dl) with(dl, (
    ggplot(roughnesses, aes(x=roughness,
                            fill=type, color=type))
    + geom_histogram(aes(y=..density..), position='dodge')
    + geom_rug()
    + geom_vline(
        xintercept=filter(roughnesses, dl$best.sample==sample.iter)$roughness,
        linetype='dashed', color='blue')
))
