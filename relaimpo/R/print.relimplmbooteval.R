"print.relimplmbooteval" <-
function (x, ...) 
{
    if (!(is(x, "relimplmbooteval") || is(x, "relimplmbootMI"))) 
        stop("x must be the output from function booteval.relimp or mianalyze.relimp")
    p <- length(slot(x, "namen")) - 1
    nlev <- length(x@level)

    ## all information on options etc. of metric calculation included in print.relimplm(x)
    print.relimplm(x)
    cat("\n", "\n", "Confidence interval information (", x@nboot, 
        "bootstrap replicates, bty=", x@bty, "):", "\n")
    offset <- max(nchar(rownames(x@mark))) + 12
    if (x@sort) 
        cat("Sorted Relative Contributions with confidence intervals:", 
            "\n", "\n")
    else cat("Relative Contributions with confidence intervals:", 
        "\n", "\n")
    if (x@rank) {
        cat(paste(rep(" ", offset + nlev * (p + 1)), collapse = ""), 
            "Lower", paste(rep(" ", 7 * nlev - 5), collapse = ""), 
            "Upper", "\n", sep = "")
    }
    else {
        cat(paste(rep(" ", offset), collapse = ""), "Lower", 
            paste(rep(" ", 7 * nlev - 5), collapse = ""), "Upper", 
            "\n", sep = "")
    }
    print(noquote(x@mark))
    cat("\n")
    if (x@rank) {
        cat("Letters indicate the ranks covered by bootstrap CIs.", 
            "\n")
        cat("(Rank bootstrap confidence intervals always obtained by percentile method)", 
            "\n")
    }
    cat("CAUTION: Bootstrap confidence intervals can be somewhat liberal.", 
        "\n")
    if (x@fixed) cat("NOTE: X-matrix has been considered as fixed for bootstrapping.", 
        "\n")

    # differences
    if (x@diff) {
        if (x@sort) 
            cat("\n", "\n", "Differences between Relative Contributions (sorted by absolute value):", 
                "\n", "\n")
        else cat("\n", "\n", "Differences between Relative Contributions:", 
            "\n", "\n")
        offset <- max(nchar(rownames(x@markdiff))) + 12
        cat(paste(rep(" ", offset + nchar(paste(x@level,collapse=" ")) + 1), collapse = ""), 
            "Lower", paste(rep(" ", 8 * nlev - 5), collapse = ""), 
            "Upper", "\n", sep = "")
        print(noquote(x@markdiff))
        cat("\n")
        cat("* indicates that CI for difference does not include 0.", 
            "\n")
        cat("CAUTION: Bootstrap confidence intervals can be somewhat liberal.", 
            "\n")
        if (x@fixed) cat("NOTE: X-matrix has been considered as fixed for bootstrapping.", 
            "\n")
    }
}

