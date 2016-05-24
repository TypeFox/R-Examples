doubleplotter <-
function(data1=ClassPatternData$result1, data2=ClassPatternData$result2, img1=ClassPatternData$demoimage1, img2=ClassPatternData$demoimage2, metric=5) {

  #--------------------------------------------------------------
  # 
  # TITLE:     doubleplotter()
  # AUTHOR:    TARMO REMMEL
  # DATE:      24 MAY 2011
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
  #--------------------------------------------------------------

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

  

  # DRAW BOXPLOTS FOR THE TWO SPCIFIED METRICS
  map1 <- as.vector(unlist(data1[metric]))
  map2 <- as.vector(unlist(data2[metric]))
  tempobj <- cbind(map1, map2)
  dimnames(tempobj)[[2]] <- c("Map 1", "Map 2")
  boxplot(tempobj, notch=TRUE, ylab=names(data1)[metric])
  
}
