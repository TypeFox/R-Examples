panelCorr <- function(data,...){
                x<-data$x
                y<-data$y
                points(x, y, pch = 16)
                chh <- par()$cxy[2]
                x1 <- min(x)
                y1 <- max(y) - chh/4
                r1 <- cor(x, y)
                text(x1, y1, paste(round(r1, 3)), cex = 0.8, adj = 0)
        }

