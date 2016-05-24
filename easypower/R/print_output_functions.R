get.highest <- function(sort.frame) {
    largest.value <- sort.frame[1, ]

    output <- list(treatment = largest.value$variable,
                   es = largest.value$eta.sq,
                   total = largest.value$sample.size,
                   group.n = ceiling(largest.value$sample.size/max(sort.frame$lev)))
    class(output) <- "highestClass"

    return(output)
}

#' @export
print.highestClass <- function(x, ...) {
    string.out <- sprintf("\nTreatment: %s\nEffect Size: %s\nTotal N: %d\nn per cell: %d\n\n",
                          x$treatment, x$es, x$total, x$group.n)
    cat(string.out)
}

# Gets and prints all of the required sample sizes
get.all <- function(final.frame) {
    group.n <- sapply(1:length(final.frame$variable), FUN = function(x) {
        final.frame$sample.size[x]/max(final.frame$lev)
    }
    )
    out.frame <- data.frame(final.frame, "n.group" = ceiling(group.n))
    output <- data.frame("Treatment" = out.frame$variable,
                         "...Effect Size" = out.frame$eta.sq,
                         "...Total" = out.frame$sample.size,
                         "...n.per.cell" = out.frame$n.group)

    return(output)
}

# Gets the data for result="select"
# The highest sample size and the results where the user
# input a numeric value for the effect size
get.select <- function(sort.frame) {
    largest.value <- sort.frame[1, ]
    es.list <- list("small", "med", "large")

    frame.select <- data.frame()
    for(i in 1:length(sort.frame$variable)) {
        if(!(sort.frame$eta.sq[i] %in% es.list) & sort.frame$variable[i] != largest.value$variable) {
            frame.select <- rbind.data.frame(frame.select, sort.frame[i, ], make.row.names = FALSE)
        }
    }

    update.frame <- data.frame(rbind(frame.select, largest.value))
    group.n <- sapply(1:length(update.frame$variable), FUN = function(x) {
        update.frame$sample.size[x]/max(update.frame$lev)
    }
    )

    out.frame <- data.frame("Treatment" = update.frame$variable,
                            "...Effect Size" = update.frame$eta.sq,
                            "...Total N" = update.frame$sample.size,
                            "...n.per.cell" = ceiling(group.n))


    return(out.frame)
}

# Gets the data from the sample.oneway function
# Formats the data for output
get.oneway <- function(sample.size, name.iv) {
    output <- list(treatment = name.iv,
                   es = sample.size$f,
                   total.n = ceiling(sample.size$n * sample.size$k),
                   group.n = ceiling(sample.size$n)
    )

    class(output) <- "oneway"
    return(output)

}

#' @export
print.oneway <- function(x, ...) {
    cat(sprintf("Effect size used in calculation: Cohen's f\n"))
    cat(sprintf("Cutoffs: small = 0.10, med = 0.25, large = 0.40\n\n"))
    cat(sprintf("Treatment: %s\nEffect Size (f): %1.3f\nTotal N: %d\nn per cell: %d\n\n",
                x$treatment, x$es, x$total.n, x$group.n))

}


