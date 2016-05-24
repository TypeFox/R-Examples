
# layout(num.rows, num.cols, mat, num.figures, col.widths, row.heights,
#        cm.widths, cm.heights, respect, respect.mat)

C_layout <- function(x) {
    dev.set(recordDev())
    # Mimic call on off-screen device (so get the right answer when
    # query off-screen device in drawing functions)
    do.call("layout", x[-1])
    dev.set(playDev())    
}
