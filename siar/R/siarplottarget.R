siarplottarget <-
function(siardata, isox, isoy, a, grps) {
        pchseq <- c(1:2, 4:20)
        for (j in 1:nrow(siardata$targets)) {
            if (siardata$numgroups != 1 & !is.null(grps)) {
                if (any(siardata$targets[j, 1] == grps)) {
                  dx <- siardata$targets[j, isox + a]
                  dex <- 2 * siardata$corrections[isox, 3]
                  dy <- siardata$targets[j, isoy + a]
                  dey <- 2 * siardata$corrections[isoy, 3]
                  siaraddcross(x = dx, y = dy, clr = "grey50",
                    upch = pchseq[siardata$targets[j, 1]])
                }
            }
            else {
                dx <- siardata$targets[j, isox + a]
                dex <- 2 * siardata$corrections[isox, 3]
                dy <- siardata$targets[j, isoy + a]
                dey <- 2 * siardata$corrections[isoy, 3]
                siaraddcross(x = dx, y = dy, clr = "grey50")
            }
        }
    }
