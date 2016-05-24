singleplotter <-
function(data=ClassPatternData$result1, img=ClassPatternData$demoimage1, metrics=c(1,5,10), rows=1, cols=3, addactual=TRUE, colour=TRUE) {

  #--------------------------------------------------------------
  # 
  # TITLE:     singleplotter()
  # AUTHOR:    TARMO REMMEL
  # DATE:      16 JULY 2013
  # CALLS:     NA
  # CALLED BY: NA
  # NEEDS:     SDMTools LIBRARY
  # NOTES:     USED TO DRAW A BOXPLOT OF THE RANGE OF EXPECTED
  #            CLASS METRIC VALUES GIVEN A PROPORTION AND A RHO
  #            THAT IS CORRECTED FOR BIAS.  IT IS POSSIBLE TO
  #            ADD THE LANDSCAPE METRIC VALUE TO THE BOXPLOT
  #            FROM THE ORIGINAL IMAGE FROM WHICH THE PARAMETERS
  #            WERE ESTIMATED TO ILLUSTRATE HOW EXPECTED OR
  #            THAT LANDSCAPE IS, GIVEN SIMULATED RESULTS.
  #            SETTING colour=FALSE RESULTS IN A BW PLOT
  #
  # FUTURE:    FIX PLOTTING YLIM TO FORCE INCLUSION OF RED DOT
  #            THIS ONLY HAPPENS WHEN THE RED DOT (ACTUAL LANDSCAPE)
  #            VALUE IS FAR BEYOND THE EXPECTATION.
  #
  #--------------------------------------------------------------

  # EXAMPLE USAGE:
  # singleplotter(data=result1, metrics=c(2,7,18,20,21,22), rows=2, cols=3, addactual=TRUE)

  # USED FOR PRODUCING BOXPLOTS OF CLASS-LEVEL LPI RESULTS
  # [1] "LOW.class"                    "LOW.n.patches"               
  # [3] "LOW.total.area"               "LOW.prop.landscape"          
  # [5] "LOW.patch.density"            "LOW.total.edge"              
  # [7] "LOW.edge.density"             "LOW.landscape.shape.index"   
  # [9] "LOW.largest.patch.index"      "LOW.mean.patch.area"         
  # [11] "LOW.sd.patch.area"            "LOW.min.patch.area"          
  # [13] "LOW.max.patch.area"           "LOW.perimeter.area.frac.dim" 
  # [15] "LOW.mean.perim.area.ratio"    "LOW.sd.perim.area.ratio"     
  # [17] "LOW.min.perim.area.ratio"     "LOW.max.perim.area.ratio"    
  # [19] "LOW.mean.shape.index"         "LOW.sd.shape.index"          
  # [21] "LOW.min.shape.index"          "LOW.max.shape.index"         
  # [23] "LOW.mean.frac.dim.index"      "LOW.sd.frac.dim.index"       
  # [25] "LOW.min.frac.dim.index"       "LOW.max.frac.dim.index"      
  # [27] "LOW.total.core.area"          "LOW.prop.landscape.core"     
  # [29] "LOW.mean.patch.core.area"     "LOW.sd.patch.core.area"      
  # [31] "LOW.min.patch.core.area"      "LOW.max.patch.core.area"     
  # [33] "LOW.prop.like.adjacencies"    "LOW.aggregation.index"       
  # [35] "LOW.lanscape.division.index"  "LOW.splitting.index"         
  # [37] "LOW.effective.mesh.size"      "LOW.patch.cohesion.index"    
  # [39] "HIGH.class"                   "HIGH.n.patches"              
  # [41] "HIGH.total.area"              "HIGH.prop.landscape"         
  # [43] "HIGH.patch.density"           "HIGH.total.edge"             
  # [45] "HIGH.edge.density"            "HIGH.landscape.shape.index"  
  # [47] "HIGH.largest.patch.index"     "HIGH.mean.patch.area"        
  # [49] "HIGH.sd.patch.area"           "HIGH.min.patch.area"         
  # [51] "HIGH.max.patch.area"          "HIGH.perimeter.area.frac.dim"
  # [53] "HIGH.mean.perim.area.ratio"   "HIGH.sd.perim.area.ratio"    
  # [55] "HIGH.min.perim.area.ratio"    "HIGH.max.perim.area.ratio"   
  # [57] "HIGH.mean.shape.index"        "HIGH.sd.shape.index"         
  # [59] "HIGH.min.shape.index"         "HIGH.max.shape.index"        
  # [61] "HIGH.mean.frac.dim.index"     "HIGH.sd.frac.dim.index"      
  # [63] "HIGH.min.frac.dim.index"      "HIGH.max.frac.dim.index"     
  # [65] "HIGH.total.core.area"         "HIGH.prop.landscape.core"    
  # [67] "HIGH.mean.patch.core.area"    "HIGH.sd.patch.core.area"     
  # [69] "HIGH.min.patch.core.area"     "HIGH.max.patch.core.area"    
  # [71] "HIGH.prop.like.adjacencies"   "HIGH.aggregation.index"      
  # [73] "HIGH.lanscape.division.index" "HIGH.splitting.index"        
  # [75] "HIGH.effective.mesh.size"     "HIGH.patch.cohesion.index"   
  

  # COMPUTE ACTUAL CLASS METRICS
  actual <- ClassStat(img)
  
  # SETUP GRAPHIC ENVIRONMENT   
  par(mfrow=c(rows,cols), pty="s")
 
  for(num in 1:length(metrics)) {
    boxplot(data[metrics[num]])
    title(names(data[metrics[num]]))
    
    # ADD DATA POINT FOR ACTUAL LPI VALUE IF DESIRED
    if(addactual == TRUE) {
       if(metrics[num] < 39) {
         if(colour) {
           # THIS IS THE NUMBER OF SIMULATED MAPS
           numsim <- length(data[[metrics[num]]])
           points(actual[1,metrics[num]], pch=21, col="red", bg="red", cex=1.75)
           actualval <- actual[1,metrics[num]]
           cat("Actual Metric Value (", names(data[metrics[num]]), "): ", actualval, "\n")
           higherthan <- sum(data[[metrics[num]]] > actualval)
           if(is.na(higherthan)) { higherthan <- 0 }
           lowerthan <- sum(data[[metrics[num]]] < actualval)
           if(is.na(lowerthan)) { lowerthan <- 0 }
           sameas <- numsim - lowerthan - higherthan
           probhighereq <- (higherthan + sameas) / numsim
           problowereq <- (lowerthan + sameas) / numsim
           cat(higherthan, " higher values, ", lowerthan, " lower values, and ", sameas, " identical values as the map.\n", sep="")
           cat("Probability of map having a value <= to expectation: P=", formatC(probhighereq, digits=4, format="f"), "\n", sep="")
           cat("Probability of map having a value >= to expectation: P=", formatC(problowereq, digits=4, format="f"), "\n\n", sep="")
         }
         else {
           # THIS IS THE NUMBER OF SIMULATED MAPS
           numsim <- length(data[[metrics[num]]])
           points(actual[1,metrics[num]], pch=21, cex=1.75)
           actualval <- actual[1,metrics[num]]
           cat("Actual Metric Value (", names(data[metrics[num]]), "): ", actualval, "\n")
           higherthan <- sum(data[[metrics[num]]] > actualval)
           if(is.na(higherthan)) { higherthan <- 0 }
           lowerthan <- sum(data[[metrics[num]]] < actualval)
           if(is.na(lowerthan)) { lowerthan <- 0 }
           sameas <- numsim - lowerthan - higherthan
           probhighereq <- (higherthan + sameas) / numsim
           problowereq <- (lowerthan + sameas) / numsim
           cat(higherthan, " higher values, ", lowerthan, " lower values, and ", sameas, " identical values as the map. \n", sep="")
           cat("Probability of map having a value <= to expectation: P=", formatC(probhighereq, digits=4, format="f"), "\n", sep="")
           cat("Probability of map having a value >= to expectation: P=", formatC(problowereq, digits=4, format="f"), "\n\n", sep="")

         }
       }
       else {
         if(colour) {
           # THIS IS THE NUMBER OF SIMULATED MAPS
           numsim <- length(data[[metrics[num]]])
           points(actual[2,metrics[num] - 38], pch=21, col="red", bg="red", cex=1.75)
           actualval <- actual[2,metrics[num] - 38]
           cat("Actual Metric Value (", names(data[metrics[num]]), "): ", actualval, "\n")
           higherthan <- sum(data[[metrics[num]]] > actualval)
           if(is.na(higherthan)) { higherthan <- 0 }
           lowerthan <- sum(data[[metrics[num]]] < actualval)
           if(is.na(lowerthan)) { lowerthan <- 0 }
           sameas <- numsim - lowerthan - higherthan
           probhighereq <- (higherthan + sameas) / numsim
           problowereq <- (lowerthan + sameas) / numsim
           cat(higherthan, " higher values, ", lowerthan, " lower values, and ", sameas, " identical values as the map. \n", sep="")
           cat("Probability of map having a value <= to expectation: P=", formatC(probhighereq, digits=4, format="f"), "\n", sep="")
           cat("Probability of map having a value >= to expectation: P=", formatC(problowereq, digits=4, format="f"), "\n\n", sep="")
         }
         else {
           # THIS IS THE NUMBER OF SIMULATED MAPS
           numsim <- length(data[[metrics[num]]])
           points(actual[2,metrics[num] - 38], pch=21, cex=1.75)
           actualval <- actual[2,metrics[num] - 38]
           cat("Actual Metric Value (", names(data[metrics[num]]), "): ", actualval, "\n")
           higherthan <- sum(data[[metrics[num]]] > actualval)
           if(is.na(higherthan)) { higherthan <- 0 }
           lowerthan <- sum(data[[metrics[num]]] < actualval)
           if(is.na(lowerthan)) { lowerthan <- 0 }
           sameas <- numsim - lowerthan - higherthan
           probhighereq <- (higherthan + sameas) / numsim
           problowereq <- (lowerthan + sameas) / numsim
           cat(higherthan, " higher values, ", lowerthan, " lower values, and ", sameas, " identical values as the map. \n", sep="")
           cat("Probability of map having a value <= to expectation: P=", formatC(probhighereq, digits=4, format="f"), "\n", sep="")
           cat("Probability of map having a value >= to expectation: P=", formatC(problowereq, digits=4, format="f"), "\n\n", sep="")
         }
       }
    } # END IF
    
  } # END FOR
  
}
