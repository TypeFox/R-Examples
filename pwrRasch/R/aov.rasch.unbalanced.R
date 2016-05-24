##########################################################################################################
#
# pwrRasch: Statistical Power Simulation for Testing the Rasch Model
#
# Internal function: ANOVA unbalanced 
#
# Authors: Takuya Yanagida <takuya.yanagida@univie.ac.at>
#  		     Jan Steinfeld <jan.steinfeld@univie.ac.at>
#
##########################################################################################################

aov.rasch.unbalanced <- function(data, group = "group", person = "person", item = "item", response = "response",
                               output = TRUE) {
  
  eval(parse(text = paste0("data$", group, " <- as.factor(data$", group, ")")))
  eval(parse(text = paste0("data$", person, " <- as.factor(data$", person, ")")))
  eval(parse(text = paste0("data$", item, " <- as.factor(data$", item, ")")))
  
  formula <- paste(response, "~", paste(group, item, paste(group, item, sep = ":"),
                                        paste0("Error(", person, " + ", person, ":", item,")"),sep = " + "))
    
  data.aov <- summary(aov(eval(parse(text = formula)), data = data))
  
  #------------------------------------------------------------------------------------------------------#
  # Output
  
  if (output == TRUE) {
    
    cat("Three-way analysis of variance with mixed classification \n\n")
    
    print(unclass(data.aov[[2]])[[1]][2, ])
      
    if (data.aov[[1]][[1]][1, "Pr(>F)"] < .05) {
      
      warning(paste0("Main effect A (group) is statistically significant, F(1, ",
                     data.aov[[1]][[1]][2, "Df"], ") = ", formatC(data.aov[[1]][[1]][1, "F value"], format = "f", digits = 3),
                     ", p = ", formatC(data.aov[[1]][[1]][1, "Pr(>F)"], format = "f", digits = 3), 
                     ", i.e. results may not be trustworthy."))                     
                     
    }
  
  }
  
  #------------------------------------------------------------------------------------------------------#
  # Return   

  F.B <- data.aov[[1]][[1]][2, "Mean Sq"] / data.aov[[2]][[1]][3, "Mean Sq"]
  p.B <- pf(F.B, data.aov[[1]][[1]][2, "Df"], data.aov[[2]][[1]][3, "Df"], lower.tail = FALSE) 
  
  restab <- matrix(unlist(c(data.aov[[1]][[1]][1, ],
                            data.aov[[1]][[1]][2, 1:3], F.B, p.B,
                            data.aov[[2]][[1]][1, ],
                            data.aov[[2]][[1]][2, ],
                            data.aov[[2]][[1]][3, ])), ncol = 5, byrow = TRUE,
                   dimnames = list(c("group", "person(group)", "item", "group:item", "Residuals"),
                                   c("DF", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")))
  
  return(invisible(restab))
  
}