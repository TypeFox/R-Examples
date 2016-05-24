# A. Multiple pages on a single screen device 
# A1. Pause for multiple pages
windows(); bringToTop(stay = TRUE)
devAskNewPage(ask = TRUE); plot(1:3, col = 'red')
devAskNewPage(ask = TRUE); plot(5:7, col = 'black')

# A2. Save and replay on a screen
windows(); bringToTop(stay = TRUE)
plot(1:10, col = 'red'); aa <- recordPlot() 
plot(20:30, col = 'black'); bb <- recordPlot()
replayPlot(aa); replayPlot(bb)
replayPlot(bb); replayPlot(aa)

# A3. Replay and save as a file
setwd("C:/aErer")
pdf(file = "test3.pdf")  # save two graphs on a single pdf file
replayPlot(aa); replayPlot(bb)
dev.off()

# A4. Manual recording and replaying
windows(); bringToTop(stay = TRUE) 
# use mouse on pull-down menus in the window: History => Recording
plot(1:10, col = 'red')
plot(rnorm(30), col = 'black')
# use mouse on pull-down menus in the window: History => Previous, Next

# B. Multiple pages on multiple screen devices for comparison
# B1. Initiate multiple windows manually
windows(); bringToTop(stay = TRUE); plot(1:3, col = 'red')
windows(); bringToTop(stay = TRUE); plot(5:7, col = 'purple')
windows(); bringToTop(stay = TRUE); plot(7:9, col = 'green')

# B2. Initiate multiple windows by loop
data <- list(da = 1:3, db = 5:7, dc = 7:9)
color <- c('red', 'purple', 'green')
for (i in 1:3) {
  windows(); plot(data[[i]], col = color[i])
}

# C. Saving multiple pages as a single file
pdf(file = "testMpage.pdf", onefile = TRUE)  # 3 graphs in 1 file
plot(rnorm(500)); plot(1:100); plot(3:20)
dev.off()

# D. Saving multiple pages as multiple files
pdf(file = "testMpage%01d.pdf", onefile = FALSE)  # 3 graphs in 3 files
plot(rnorm(500)); plot(1:100); plot(3:20)
dev.off()

png(file = "testNpage.png")  # 1 file for the last graph; usually an error
plot(rnorm(500)); plot(1:100); plot(3:20)
dev.off()

png(file = "testNpage%01d.png")  # 3 graphs in 3 files
plot(rnorm(500)); plot(1:100); plot(3:20)
dev.off()

# using a loop statement
dnor <- list(da = rnorm(30), db = rnorm(300), dc = rnorm(1000))
color <- c('red', 'purple', 'green')
for (i in 1:3) {
  png(file = paste("testGpage", i, ".png", sep = ""))  
  plot(dnor[[i]], col = color[i])
  dev.off()  
}  