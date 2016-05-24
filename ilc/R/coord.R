coord <-
function(pos){
    # find current plotting area coordinates:
    x <- par("usr")  
    # coordonates of the legends
    if(mode(pos) == "character") {
        test.arg <- c("UR", "LR", "UL", "LL", "UC", "LC")
        # U - upper, L - lower, R - right , L - left, C - centre
        pos <- check.arg(pos, test.arg, n = 1, value = F)
        pos <- switch(pos,
                      c(0.67, 0.98),  # "UR"
                      c(0.67, 0.2),   # "LR"
                      c(0.015, 0.98), # "UL"
                      c(0.015, 0.2),  # "LL"
                      c(0.35, 0.98),  # "UC"
                      c(0.35, 0.2))   # "LC"
    }
    if(mode(pos) == "list") pos <- c(pos$x, pos$y)
    list(x = x[1] + pos[1] * (x[2] - x[1]), y = x[3] + pos[2] * (x[4] - x[3]))
}
