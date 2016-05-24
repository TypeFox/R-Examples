plotTrimmOutl <- function(data,trimmOutl,nsizes,bustVariable,variable,col,xlim,ylim,main){

 if(variable == "chest"){
  plot(data[, bustVariable], data[, variable], pch = "*", col = "thistle1", xlab = bustVariable, ylab = variable, 
       main = main, xlim = xlim, ylim = ylim, xaxt = "n", yaxt = "n")
  axis(1, at = seq(xlim[1], xlim[2], 10), labels = seq(xlim[1], xlim[2], 10))
  axis(2, at = seq(ylim[1], ylim[2], 10), labels = seq(ylim[1], ylim[2], 10))

  for(i in 1 : nsizes){
   points(data[as.character(trimmOutl[[i]]), bustVariable], data[as.character(trimmOutl[[i]]), variable], pch = i, 
          col = col[i])
  }
 }

 if(variable == "hip"){
  plot(data[, bustVariable], data[, variable], pch = "*", col = "thistle1", xlab = bustVariable, ylab = variable, 
       main = main, xlim = xlim, ylim = ylim, xaxt = "n", yaxt = "n")
  axis(1, at = seq(xlim[1], xlim[2], 10), labels = seq(xlim[1], xlim[2], 10))
  axis(2, at = seq(ylim[1], ylim[2], 10), labels = seq(ylim[1], ylim[2], 10))

  for(i in 1 : nsizes){
   points(data[as.character(trimmOutl[[i]]), bustVariable], data[as.character(trimmOutl[[i]]), variable], pch = i, 
          col = col[i])
  }
 }

 if(variable == "necktoground"){
  plot(data[, bustVariable], data[, variable], pch = "*", col = "thistle1", xlab = bustVariable, ylab = variable, 
       main = main, xlim = xlim, ylim = ylim, xaxt = "n", yaxt = "n")
  axis(1, at = seq(xlim[1], xlim[2], 10), labels = seq(xlim[1], xlim[2], 10))
  axis(2, at = seq(ylim[1], ylim[2], 10), labels = seq(ylim[1], ylim[2], 10))

  for(i in 1 : nsizes){
   points(data[as.character(trimmOutl[[i]]), bustVariable], data[as.character(trimmOutl[[i]]), variable], pch = i, 
          col = col[i])
  }
 }

 if(variable == "waist"){
  plot(data[, bustVariable], data[, variable], pch = "*", col = "thistle1", xlab = bustVariable, ylab = variable, 
       main = main, xlim = xlim, ylim = ylim, xaxt = "n", yaxt = "n")
  axis(1, at = seq(xlim[1], xlim[2], 10), labels = seq(xlim[1], xlim[2], 10))
  axis(2, at = seq(ylim[1], ylim[2], 10), labels = seq(ylim[1], ylim[2], 10))

  for(i in 1 : nsizes){
   points(data[as.character(trimmOutl[[i]]), bustVariable], data[as.character(trimmOutl[[i]]), variable], pch = i, 
          col = col[i])
  }
 }

 if((variable != "chest") & (variable != "hip") & (variable != "necktoground") & (variable != "waist")){ 
   stop("This variable doesn't belong to the database")
 }

}
