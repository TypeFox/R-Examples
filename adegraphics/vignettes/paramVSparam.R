Adegpar <- adegpar()
names <- names(Adegpar)
xx <- as.null()
yy <- as.null()

for(i in 1:length(Adegpar)) {
  if(is.list(Adegpar[[i]])) {
    for(j in 1:length(Adegpar[[i]])) {
      if(is.list(Adegpar[[i]][[j]])) {  ## sublist of list
        xx <- c(xx, paste(names(Adegpar)[i], '.', names(Adegpar[[i]])[j], sep = ""))
        yy <- c(yy, names(Adegpar[[i]][[j]]))
      } else {
        yy <- c(yy, names(Adegpar[[i]])[j])
        xx <- c(xx, names[i])
      }
    }
  }
  else
    xx <- c(xx, names[i])
}

yy <- unique(yy)
xx <- unique(xx)
paramVSsub <- data.frame(matrix(0, nrow = length(yy), ncol = length(xx)))
row.names(paramVSsub) <- yy
colnames(paramVSsub) <- xx

## filling
for(i in 1:length(Adegpar)) {
  if(is.list(Adegpar[[i]])) {
    for(j in 1:length(Adegpar[[i]])) {
      if(is.list(Adegpar[[i]][[j]])) ## sublistof list
        paramVSsub[names(Adegpar[[i]][[j]]), paste(names(Adegpar)[i], '.', names(Adegpar[[i]])[j], sep = "")] <- 100
      else
        paramVSsub[names(Adegpar[[i]])[j], names[i]] <- 100
  	}
  }
}

table.value(t(paramVSsub), axis.text = list(cex = 0.9), symbol = "circle", plegend.drawKey = FALSE, ppoints.cex = 0.3, 
            ptable.y = list(srt = 60, pos = "left"),
            ptable.margin = list(bottom = 2, left = 17, top = 17, right = 2))
            



