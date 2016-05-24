### The plot method for the `prinsimp` class produces a scree plot: a
### bar or line plot of the variances of each simple basis, same as
### the plot method for classes `princomp` and `prcomp`.
###
### See: screeplot.default
plot.prinsimp <- function(x, main = deparse(substitute(x)), ...) {
    screeplot(x, main = main, ...)
}

plot.simpart <- function(x, display=list(model=TRUE, simple=TRUE), layout, ...) {
    if (missing(layout)) {
        model_basis <- min(length(x$variance$model),
                             length(x$variance$model[display$model]))
        simple_basis <- min(length(x$variance$simple),
                             length(x$variance$simple[display$simple]))
        basis_to_plot <- model_basis + simple_basis
        if (basis_to_plot > 6) {
            stop("For default plot, the total number of basis vector plots
cannot exceed 6. Use 'display=list(model=...., simple=...)' to
reduce the number of basis vectors plotted or specify the plot
layout manually, using the 'layout' argument")
        }

        layout <- matrix(0, nrow=3, ncol=4)

        ## the first column consists of the model basis vectors plus
        ## the wrap-around tail of the simple basis, if it exists
        left_col_model <- min(3, model_basis)
        left_col_simple <- max(0, simple_basis - 3)
        layout[seq_len(left_col_model)] <- seq_len(left_col_model)
        layout[4 - seq_len(left_col_simple)] <- model_basis + seq(to = simple_basis,
                                                                  len = left_col_simple)
        
        ## the second column consists of the simple basis vectors plus
        ## the wrap-around tail of the model basis, if it exists
        right_col_simple <- min(3, simple_basis)
        right_col_model <- max(0, model_basis - 3)
        layout[3+seq_len(right_col_simple)] <- model_basis + seq_len(right_col_simple)
        layout[7 - seq_len(right_col_model)] <- seq(to = model_basis,
                                                    len = right_col_model)

        if (all(layout[ , 1] == 0)) layout[ , 1] <- layout[ , 2]
        else if (all(layout[ , 2] == 0)) layout[ , 2] <- layout[ , 1]
        
        layout[, 3:4] <- c(basis_to_plot+1, basis_to_plot+1, basis_to_plot+2,
                           basis_to_plot+1, basis_to_plot+1, basis_to_plot+2)
    }
    
    def.par <- par(no.readonly = TRUE) # save default, for resetting...
    on.exit(par(def.par))
    layout(layout)
    
    basisplot(x, display, ...)
    varsimp(x, display, ...)
    varperc(x, ...)
}


varsimp <- function(x, display=list(model=TRUE, simple=TRUE), full.simple=TRUE, ...) {
    model_simplicity <- x$simplicity$model
    simple_simplicity <- x$simplicity$simple
    full_simplicity <- x$simplicity$full

    d <- length(full_simplicity)
    nullnum <- length(simple_simplicity)
    
    vpcrelval <- x$variance$model / sum(x$variance$full)
    vpcrelval <- vpcrelval[seq_len(d - nullnum)]
    
    vsmoothrelval <- x$variance$simple / sum(x$variance$full)

    plot(NULL, xlim=c(0, max(vpcrelval, vsmoothrelval)),
         ylim=range(0, unlist(x$simplicity), finite = TRUE),
         main='Variance-simplicity View',
         xlab='Proportion of Var Explained',
         ylab='Simplicity score (complex to simple)')

    if (length(vsmoothrelval) > 0) {
        to_plot <- display$simple
        if (length(vsmoothrelval[to_plot]) > 0) {
            points(vsmoothrelval[to_plot], simple_simplicity[to_plot], col='red')
            text(vsmoothrelval[to_plot], simple_simplicity[to_plot],
                 labels=as.character(seq(vsmoothrelval)[to_plot]),
                 pos=4, col='red')
        }
    }
    
    if (length(vpcrelval) > 0) {
        to_plot <- display$model
        if (length(vpcrelval[to_plot]) > 0) {
            points(vpcrelval[to_plot], model_simplicity[to_plot], pch='*', col='blue')
            text(vpcrelval[to_plot], model_simplicity[to_plot],
                 labels=as.character(seq(vpcrelval)[to_plot]),
                 pos=4, col='blue')
        }
    }

    if (full.simple) {
        abline(h=full_simplicity, lty='dashed')
    }
}


basisplot <- function(x, display=list(model=TRUE, simple=TRUE), ...) {
    model_basis <- x$model
    simple_basis <- x$simple

    d <- nrow(model_basis)
    npcout <- ncol(model_basis)
    nullnum <- ncol(simple_basis)

    model_var_explained <- x$variance$model / sum(x$variance$full) * 100
    simple_var_explained <- x$variance$simple / sum(x$variance$full) * 100

    y_range <- range(model_basis[, seq_len(npcout)[display$model]],
                     simple_basis[, seq_len(nullnum)[display$simple]],
                     0, finite = TRUE)
    
    for (i in seq_len(npcout)) {
        if (!(i %in% seq_len(npcout)[display$model])) next;
        plot(rownames(model_basis), model_basis[,i],
             type='l', lty='solid', col='blue',
             xlab = 'x', ylab = 'Loadings',
             ylim = y_range)
        abline(h=0, lty='dashed')

        ## probe index in the top-left corner
        text(par('usr')[1], par('usr')[4], as.character(i),
             col='blue', adj=c(-.5, 1.5))

        ## %-variance explained in the top-right corner
        var_explained <- format_varperc(model_var_explained[i], 2)
        
        text(par('usr')[2], par('usr')[4],
             paste0('var: ', var_explained, '%'),
             adj=c(1, 1.5))
    }
    for (i in seq_len(nullnum)) {
        if (!(i %in% seq_len(nullnum)[display$simple])) next;
        plot(rownames(simple_basis), simple_basis[,i],
             type='l', lty='dashed', col='red',
             xlab = 'x', ylab = 'Loadings',
             ylim = y_range)
        abline(h=0, lty='dashed')

        ## probe index in the top-left corner
        text(par('usr')[1], par('usr')[4], as.character(i),
             col='red', adj=c(-.5, 1.5))

        ## %-variance explained in the top-right corner
        var_explained <- format_varperc(simple_var_explained[i], 2)
        text(par('usr')[2], par('usr')[4],
             paste0('var: ', var_explained, '%'),
             adj=c(1, 1.5))
    }
}


varperc <- function(x, ...) {
    simple_proportion <- sum(x$variance$simple) / sum(x$variance$full)

    plot(NULL, xlim=c(0, 1), ylim=c(-1, 1),
         main='% Var Explained')
    lines(c(0, simple_proportion), c(0, 0), lwd=30, col='red', lend='butt')
    lines(c(simple_proportion, 1), c(0, 0), lwd=30, col='blue', lend='butt')
    lines(rep(simple_proportion, 2), c(-.5, .5), lwd=2)

    model_label <- paste('Model% =', sprintf("%.0f", (1-simple_proportion)*100))
    null_label <- paste('Null% =', sprintf("%.0f", simple_proportion*100))
    text(simple_proportion, .5, model_label,
         col='blue', pos=position_text(model_label, simple_proportion, 3))
    ## calculate position using the model label because it's longer,
    ## ensuring that both labels are aligned in the same direction
    text(simple_proportion, -.5, null_label,
         col='red', pos=position_text(model_label, simple_proportion, 1))
}


## Returns the position specifier that can be used as arguments to the
## `text` function to place a string centered (pos=2 or 4) at the
## specified location if it will fit within the plot boundaries, or
## anchor it to the left (pos=1) or right (pos=3) if the location is
## too close to the right or left boundary of the plot, respectively.
position_text <- function(label, x, desired_pos) {
    plot_margins <- par('usr')
    left_margin <- plot_margins[1]
    right_margin <- plot_margins[2]
    
    width <- strwidth(label)
    midpoint <- width/2
    
    if (x - midpoint < left_margin)
        4
    else if (midpoint + x > right_margin)
        2
    else desired_pos
}
