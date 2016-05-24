gx.hypergeom <-
function (tt, aa, kk, xx) 
{
    pp <- 100 * round(dhyper(xx, aa, tt - aa, kk), 3)
    cat("  Traverse length is", tt, "sites with", aa, "sites 'expected' a priori to be anomalous,", 
        "\n  number of >Threshold sites is", kk, "with", xx, 
        "coinciding at 'expected' sites.", paste("\n  Probability that this is due to 'chance' is ", 
            pp, "%\n", sep = ""))
    invisible()
}
