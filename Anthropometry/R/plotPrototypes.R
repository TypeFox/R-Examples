plotPrototypes <- function(data,prototypes,nsizes,bustVariable,variable,col,xlim,ylim,main,EN){

 #Bust intervals defined by the European standard to sizing system:
 Bust_4 <- seq(76, 104, 4)  ; Bust_6 <- seq(110, 128, 6)  ; BustVec <- c(Bust_4, Bust_6) 

 if(variable == "chest"){
  plot(data[, bustVariable], data[, variable], pch = "*", col = "thistle1", xlab = bustVariable, ylab = variable, 
       main = main, xlim = xlim, ylim = ylim, xaxt = "n", yaxt = "n")
  axis(1, at = seq(xlim[1], xlim[2], 10), labels = seq(xlim[1], xlim[2], 10))
  axis(2, at = seq(ylim[1], ylim[2], 10), labels = seq(ylim[1], ylim[2], 10))

  for(i in 1 : nsizes){
  #To locate correctly the rows of the prototypes in the whole database, the prototypes labels must be 
  #converted into a character.
   points(data[as.character(prototypes[[i]]), bustVariable], data[as.character(prototypes[[i]]), variable], pch = i, 
          col = col[i])
  }

  if(EN){
   #The European standard to sizing system does not fix the chest standard measurements. In order to 
   #overcome this limitation, we round the values obtained by means of a linear regression 
   #(see Ibanez et al. (2012)), to the nearest integer:
   Chest_4 <- seq(79, 107, 4) ; Chest_6 <- seq(112, 130, 6) ; Chest <- c(Chest_4, Chest_6) 

   for(i in 1:length(Bust_4)){
    symbols(Bust_4[i], Chest_4[i], rectangles = matrix(c(4,4), 1, 2), add = TRUE, inches = FALSE)
    points(Bust_4[i], Chest_4[i], pch = ".", cex = 2.3)
   }
   for(i in 1:length(Bust_6)){
    symbols(Bust_6[i] - 1, Chest_6[i] - 1, rectangles = matrix(c(6,6),1,2), add = TRUE, inches = FALSE)
    points(Bust_6[i] - 1, Chest_6[i] - 0.5, pch = ".", cex = 2.3)
   }
  } 

 }


 if(variable == "hip"){
  plot(data[, bustVariable], data[, variable], pch = "*", col = "thistle1", xlab = bustVariable, ylab = variable, 
       main = main, xlim = xlim, ylim = ylim, xaxt = "n", yaxt = "n")
  axis(1, at = seq(xlim[1], xlim[2], 10), labels = seq(xlim[1], xlim[2], 10))
  axis(2, at = seq(ylim[1], ylim[2], 10), labels = seq(ylim[1], ylim[2], 10))

  for(i in 1 : nsizes){
   points(data[as.character(prototypes[[i]]), bustVariable], data[as.character(prototypes[[i]]), variable], pch = i, 
          col = col[i])
  }

  if(EN){ 
   #Hip intervals defined by the European standard to sizing system.
   Hip_4 <- seq(84,112,4) ; Hip_5 <- seq(117,132,5) ; Hip = c(Hip_4,Hip_5) 

   for(i in 1:length(Bust_4)){
    symbols(Bust_4[i], Hip_4[i], rectangles = matrix(c(4,4), 1, 2), add = TRUE, inches = FALSE)
    points(Bust_4[i], Hip_4[i], pch = ".", cex = 2.3)
   }

    for(i in 1:length(Bust_6)){
     symbols(Bust_6[i] - 1, Hip_5[i] - 0.5, rectangles = matrix(c(6,5),1,2), add = TRUE, inches = FALSE)
     points(Bust_6[i] - 1, Hip_5[i] - 0.5, pch = ".", cex = 2.3)
    }
  }

 }


 if(variable == "necktoground"){
  plot(data[, bustVariable], data[, variable], pch = "*", col = "thistle1", xlab = bustVariable, ylab = variable, 
       main = main, xlim = xlim, ylim = ylim, xaxt = "n", yaxt = "n")
  axis(1, at = seq(xlim[1], xlim[2], 10), labels = seq(xlim[1], xlim[2], 10))
  axis(2, at = seq(ylim[1], ylim[2], 10), labels = seq(ylim[1], ylim[2], 10))

  for(i in 1 : nsizes){
   points(data[as.character(prototypes[[i]]), bustVariable], data[as.character(prototypes[[i]]), variable], pch = i, 
          col = col[i])
  }

  if(EN){
   #We take as neck to ground measures for the standard sizing system, the values 132, 136 and 140 cm 
   #because those are the most repeated measurements:
   vec <- seq(130, 142, 4)
   Bust_ng_4 <- seq(74, 102, 4)  ; Bust_ng_6 <- seq(107, 131, 6)  ; Bust_ng <- c(Bust_ng_4, Bust_ng_6) 

    for(i in 1:length(Bust_ng)){
     segments(Bust_ng[i], vec[1], Bust_ng[i], vec[length(vec)])
   }

    for(i in 1:length(vec)){
     segments(Bust_ng[1], vec[i], Bust_ng[length(Bust_ng)], vec[i])
    }
  }

 }

 
 if(variable == "waist"){
  plot(data[, bustVariable], data[, variable], pch = "*", col = "thistle1", xlab = bustVariable, ylab = variable, 
       main = main, xlim = xlim, ylim = ylim, xaxt = "n", yaxt = "n")
  axis(1, at = seq(xlim[1], xlim[2], 10), labels = seq(xlim[1], xlim[2], 10))
  axis(2, at = seq(ylim[1], ylim[2], 10), labels = seq(ylim[1], ylim[2], 10))

  for(i in 1 : nsizes){
   points(data[as.character(prototypes[[i]]), bustVariable], data[as.character(prototypes[[i]]), variable], pch = i, 
          col = col[i])
  }

  if(EN){
   #Waist intervals defined by the European standard to sizing system:
   Waist_4 <- seq(60,88,4) ; Waist_6 <- seq(94,112,6) ; Waist = c(Waist_4,Waist_6) 

   for(i in 1:length(Bust_4)){
    symbols(Bust_4[i], Waist_4[i], rectangles = matrix(c(4,4), 1, 2), add = TRUE, inches = FALSE)
    points(Bust_4[i], Waist_4[i], pch = ".", cex = 2.3)
   }

   for(i in 1:length(Bust_6)){
    symbols(Bust_6[i] - 1, Waist_6[i] - 1, rectangles = matrix(c(6,6),1,2), add = TRUE, inches = FALSE)
    points(Bust_6[i] - 1, Waist_6[i] - 0.5, pch = ".", cex = 2.3)
   }
  }

 }

 if((variable != "chest") & (variable != "hip") & (variable != "necktoground") & (variable != "waist")){ 
   stop("This variable doesn't belong to the database")
 }
}
