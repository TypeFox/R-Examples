##' Plot semi-gSEM principle 1 result.
##' 
##' plot.sgSEMp1 plots a structural equation network model diagram based on the best functional form for each selected pairwise variable. 
##' @title Plotting of Principle 1 of Semi-gSEM
##' @param x The returned list from sgSEMp1. Plotting uses the first element of this list (res.print) in which the first column of it is the response, the second column is variable and the other columns are the corresponding best functional form, r-squared, adj-r-squared, P-value1, P-value2 and P-value3.
##' @param cutoff A threshold value for the adjusted R-squared. Solid lines represent a relationship with the adjusted R-sqr greater than the cutoff and dotted lines with less than the cutoff. The default is 0.2.
##' @param width A numeric describing the width of the plot output in pixels.
##' @param height A numeric describing the height of the plot output in pixels.
##' @param filename A character string naming a file to save as an html file.
##' @param ... Other parameters. Currently not used. 
##' @return An html style plot of the pairwise relationship pathway diagram between stressors and responses. Arrows show relationships between each variable with given statistical relations along the connection lines.
##'   
##' @export
##' 
##' @examples
##' # Load acrylic data set
##' data(acrylic)
##' # Build a semi-gSEM model with principle 1
##' ans <- sgSEMp1(acrylic)
##' # Plot the network model with adjusted-R-squred of 0.1
##' plot(ans, cutoff = 0.1)

plot.sgSEMp1 <- function(x, ...,
                         cutoff = 0.2,
                         width = NULL,
                         height = NULL,
                         filename = NULL){
  
  rtp1 <- x$table
  rtp1[, -c(1:3)] <- round(rtp1[, -c(1:3)], 2)
  rtp1.a <- rtp1[rtp1[, "adj-R-Sqr"] < cutoff, ]
  rtp1.b <- rtp1[rtp1[, "adj-R-Sqr"] >= cutoff, ]
  
  # This generates syntax for connections between variables and responses. 
  # <br/> is for linebreak between AIC values.
  
  # node styling options:
  # [] for rectanguler, () for rounded edges in rectangle, (( )) for circle, {} for rhombus
  # Details in "http://knsv.github.io/mermaid/flowchart.html"
  
  if(dim(rtp1.a)[1] > 0) {  
      conp1.a <- sapply(1 :nrow(rtp1.a),
                        function(i){  
                            paste0(rtp1.a[i,2],
                                   "(", rtp1.a[i,2],
                                   ")", "-.->|", 
                                   paste0(colnames(rtp1.a[,3:5]),
                                          ":", rtp1.a[i,3:5],
                                          collapse="<br/>"),
                                   "|", rtp1.a[i,1], "(",
                                   rtp1.a[i,1], ")")
                        }
                        )
  }
  
  if(dim(rtp1.b)[1] > 0) {
      conp1.b <- sapply(1:nrow(rtp1.b),
                        function(i){  
                            paste0(rtp1.b[i,2],
                                   "(", rtp1.b[i,2],
                                   ")", "==>|", 
                                   paste0(colnames(rtp1.b[,3:5]),
                                          ":", rtp1.b[i,3:5],
                                          collapse="<br/>"),
                                   "|", rtp1.b[i,1], "(", rtp1.b[i,1], ")")
                        }
                        ) 
  }
  
  ## This generates syntax to run "mermaid" for plotting using above syntax
  ## "LR" is left to right flow
  
  ## For "fill" and "stroke", CSS style coloring can be used. 
  
  if(exists("conp1.a") == TRUE & exists("conp1.b") == TRUE) {
      conp1.plot <- paste0(
          "graph LR;", "\n", 
          paste(conp1.a, collapse = "\n"), "\n", paste(conp1.b, collapse="\n"), "\n", 
          "classDef default fill:#FFFF99, stroke:#000000, stroke-width:3px;")
          }
  
  if(exists("conp1.a") == FALSE) { 
      conp1.plot <- paste0(
          "graph LR;", "\n",
          paste(conp1.b, collapse="\n"), "\n",  
          "classDef default fill:#FFFF99, stroke:#000000, stroke-width:3px;")
              cat("The cutoff value is lower than all of the adjusted R-sqr values: Only solid lines")
  }
  
  
  
  if(exists("conp1.b") == FALSE) { 
      conp1.plot <- paste0(
          "graph LR;", "\n",
          paste(conp1.a, collapse="\n"), "\n",  
          "classDef default fill:#FFFF99, stroke:#000000, stroke-width:3px;")
              cat("The cutoff value is higher than all of the adjusted R-sqr values: Only dotted lines\n")
  }
  
  
  p1 <- DiagrammeR::mermaid(conp1.plot, width=width, height=height) 
  
  if(!is.null(filename))
      saveWidget(p1, file = filename, selfcontained = TRUE)
  
  return(p1)
}
