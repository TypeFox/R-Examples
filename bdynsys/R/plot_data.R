# plot_data creates a phase plane plotting the data for two variables in dynamic interaction
# argument entidxn defines the selected cases to be highlightened in the data_plot
# it corresponds to the line number in the wide format 

# to call: plot_data(datap, datap$logGDP, datap$DemocrH, seq(0, 12, by = 0.5), seq(0, 1, by = 0.1), 1, 2, 4, 5, 7, 9)

plot_data <- function(dataset, xvar, yvar, rangeX, rangeY, entidx1, entidx2, entidx3, entidx4, entidx5, entidx6) 
{ 
  procdata <- preprocess_data(2, xvar, yvar)
  xwide <- procdata$xwide
  ywide <- procdata$ywide
  
  tmpx <- rangeX
  tmpy <- rangeY
  xmin=min(tmpx)
  xmax=max(tmpx)
  ymin=min(tmpy)
  ymax=max(tmpy)
  
  # phase portrait with highlighted data trajectories
  dev.set(1)
  postscript("dataplot.eps", horizontal=FALSE, width=5, height=5,
             onefile=FALSE, paper="special", family="ComputerModern")
  
 # title(xlab=strcat(xvar, "(x)", ylab=strcat(yvar, "(y)"), pointsize=18))
  plot(xwide[1:1, 1:20], ywide[1:1, 1:20], type='l', lwd=1.5, col='gray', xlab = "X-Variable", ylab = "Y-Variable", 
       xlim=c(xmin, xmax), ylim=c(ymin, ymax))
  #points(xwide[1:20, 10], ywide[1:20, 10], pch = 20, col='darkgray')
  
  # plotting 15 random cases as trajectories and highligh selected entities' tranjectories  
  matplot(xwide[10,], ywide[10,], type='l', ldw=2, col = 'gray', add=TRUE)
  matplot(xwide[11,], ywide[11,], type='l', ldw=2, col = 'gray', add=TRUE)
  matplot(xwide[12,], ywide[12,], type='l', ldw=2, col = 'gray', add=TRUE)
  matplot(xwide[13,], ywide[13,], type='l', ldw=2, col = 'gray', add=TRUE)
  matplot(xwide[14,], ywide[14,], type='l', ldw=2, col = 'gray', add=TRUE)
  matplot(xwide[15,], ywide[15,], type='l', ldw=2, col = 'gray', add=TRUE)
  matplot(xwide[16,], ywide[16,], type='l', ldw=2, col = 'gray', add=TRUE)
  matplot(xwide[17,], ywide[17,], type='l', ldw=2, col = 'gray', add=TRUE)
  matplot(xwide[18,], ywide[18,], type='l', ldw=2, col = 'gray', add=TRUE)
  matplot(xwide[19,], ywide[19,], type='l', ldw=2, col = 'gray', add=TRUE)
  matplot(xwide[20,], ywide[20,], type='l', ldw=2, col = 'gray', add=TRUE)
  matplot(xwide[21,], ywide[21,], type='l', ldw=2, col = 'gray', add=TRUE)
  matplot(xwide[22,], ywide[22,], type='l', ldw=2, col = 'gray', add=TRUE)
  matplot(xwide[23,], ywide[23,], type='l', ldw=2, col = 'gray', add=TRUE)
  matplot(xwide[24,], ywide[24,], type='l', ldw=2, col = 'gray', add=TRUE)
  matplot(xwide[25,], ywide[25,], type='l', ldw=2, col = 'gray', add=TRUE)

  matplot(xwide[entidx1,], ywide[entidx1,], type='l', ldw=2, col = 'blue', add=TRUE)
  matplot(xwide[entidx2,], ywide[entidx2,], type='l', ldw=2, col = 'darkgreen', add=TRUE)
  matplot(xwide[entidx3,], ywide[entidx3,], type='l', ldw=2, col = 'red', add=TRUE)
  matplot(xwide[entidx4,], ywide[entidx4,], type='l', ldw=2, col = 'cyan', add=TRUE)
  matplot(xwide[entidx5,], ywide[entidx5,], type='l', ldw=2, col = 'magenta', add=TRUE)
  matplot(xwide[entidx6,], ywide[entidx6,], type='l', ldw=2, col = 'black', add=TRUE)
  
  # adding markers (circles) for marking initial conditions (default)
  points(xwide[10,1], ywide[10,1], pch = 20, col = 'darkgray')
  points(xwide[11,1], ywide[11,1], pch = 20, col = 'darkgray')
  points(xwide[12,1], ywide[12,1], pch = 20, col = 'darkgray')
  points(xwide[13,1], ywide[13,1], pch = 20, col = 'darkgray')
  points(xwide[14,1], ywide[14,1], pch = 20, col = 'darkgray')
  points(xwide[15,1], ywide[15,1], pch = 20, col = 'darkgray')
  points(xwide[16,1], ywide[16,1], pch = 20, col = 'darkgray')
  points(xwide[17,1], ywide[17,1], pch = 20, col = 'darkgray')
  points(xwide[18,1], ywide[18,1], pch = 20, col = 'darkgray')
  points(xwide[19,1], ywide[19,1], pch = 20, col = 'darkgray')
  points(xwide[20,1], ywide[20,1], pch = 20, col = 'darkgray')
  points(xwide[21,1], ywide[21,1], pch = 20, col = 'darkgray')
  points(xwide[22,1], ywide[22,1], pch = 20, col = 'darkgray')
  points(xwide[23,1], ywide[23,1], pch = 20, col = 'darkgray')
  points(xwide[24,1], ywide[24,1], pch = 20, col = 'darkgray')
  points(xwide[25,1], ywide[25,1], pch = 20, col = 'darkgray')

  points(xwide[entidx1,1], ywide[entidx1,1], pch = 20, col = 'blue')
  points(xwide[entidx2,1], ywide[entidx2,1], pch = 20, col = 'darkgreen')
  points(xwide[entidx3,1], ywide[entidx3,1], pch = 20, col = 'red')
  points(xwide[entidx4,1], ywide[entidx4,1], pch = 20, col = 'cyan')
  points(xwide[entidx5,1], ywide[entidx5,1], pch = 20, col = 'magenta')
  points(xwide[entidx6,1], ywide[entidx6,1], pch = 20, col = 'black')
  
  legend("topleft", bg="white", legend=c('Entity1  ', 'Entity2  ', 'Entity3', 'Entity4', 'Entity5', 'Entity6'), 
         lwd = 2, col=c('blue', 'darkgreen', 'red', 'cyan', 'magenta', 'black'), bty='n')
        
  dev.off()
}
