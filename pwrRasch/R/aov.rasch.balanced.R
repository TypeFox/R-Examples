##########################################################################################################
#
# pwrRasch: Statistical Power Simulation for Testing the Rasch Model
#
# Internal function: ANOVA balanced 
#
# Authors: Takuya Yanagida <takuya.yanagida@univie.ac.at>
#  		     Jan Steinfeld <jan.steinfeld@univie.ac.at>
#
##########################################################################################################

aov.rasch.balanced <- function(data, group = "group", person = "person", item = "item", response = "response",
                               output = TRUE) {
  
  eval(parse(text = paste0(group, " <- data$", group)))
  eval(parse(text = paste0(person, " <- data$", person)))
  eval(parse(text = paste0(item, " <- data$", item)))
  eval(parse(text = paste0(response, " <- data$", response)))
  
  a <- length(unique(group))
  b <- length(unique(person))
  c <- length(unique(item))
  
  ###
  
  SS.T <- (1/(a*b/a*c))*sum(response)^2

  SS.A <- sum(rowsum(response, group)^2)*(1/(b/a*c)) - SS.T
  
  SS.B <- sum(rowsum(response, person)^2)*(1/c) - SS.T - SS.A
  
  SS.C <- sum(rowsum(response, item)^2)*(1/(a*b/a)) - SS.T

  SS.interm. <- c(rowsum(response[group == 1], item[group == 1]),rowsum(response[group == 2], item[group == 2]))
  
  SS.AC <- sum(SS.interm.^2)*(1/(b/a)) - SS.A - SS.C - SS.T
  
  SS.BC <- sum(SS.interm.) - SS.T - SS.B - SS.C - SS.AC - SS.A  
  
  ###
  
  df.B <- a*(b/a - 1)
  
  df.C <- c - 1
  
  df.BC <- df.B*df.C
    
  ###
  
  MS.A <- SS.A 
  
  MS.B <- SS.B / df.B
  
  MS.C <- SS.C / df.C
  
  MS.BC <- SS.BC / df.BC
 
  MS.AC <- SS.AC / df.C
  
  ###
  
  F.A <- MS.A / MS.B
  
  F.B <- MS.B / MS.BC
  
  F.C <- MS.C / MS.BC
  
  F.AC <- MS.AC / MS.BC
  
  ###
    
  p.A <- pf(F.A, 1, df.B, lower.tail = FALSE)
  
  p.B <- pf(F.B, df.B, df.BC, lower.tail = FALSE)
  
  p.C <- pf(F.C, df.C, df.BC, lower.tail = FALSE)
  
  p.AC <- pf(F.AC, df.C, df.B*df.C, lower.tail = FALSE)
  
  
  #------------------------------------------------------------------------------------------------------#
  # Output
  
  if (output == TRUE) {
    
    cat("Three-way analysis of variance with mixed classification \n\n")
        
    printCoefmat(matrix(c(df.C, SS.AC, MS.AC, F.AC, p.AC), ncol = 5,
                        dimnames = list("group:item", c("DF", "Sum Sq", "Mean Sq", "F value", "Pr(>F)"))),
                 has.Pvalue = TRUE, cs.ind = 0)
    
    # Warning: Statistically significant Main Effect A
    if (p.A < .05) {
      
      warning(paste0("Main effect A (group) is statistically significant, ", 
                     "F(1, ", df.B, ") = ", formatC(F.A, format = "f", digits = 3),  
                     ", p = ", formatC(p.A, format = "f", digits = 3),  ", i.e. results may not be trustworthy."))
      
    }
  
  }
  
  #------------------------------------------------------------------------------------------------------#
  # Return 
  
  restab <- matrix(c(1, SS.A, MS.A, F.A, p.A,
                     df.B, SS.B, MS.B, F.B, p.B,
                     df.C, SS.C, MS.C, F.C, p.C,
                     df.C, SS.AC, MS.AC, F.AC, p.AC,
                     df.BC, SS.BC, MS.BC, NA, NA), ncol = 5, byrow = TRUE,
                   dimnames = list(c("group", "person(group)", "item", "group:item", "Residuals"),
                                   c("DF", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")))

  return(invisible(restab))
  
}