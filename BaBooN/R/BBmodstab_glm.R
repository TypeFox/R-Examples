# Bootstrap model stabilizer for linear models
# Version:            0.2
# Date:        2015-05-17
# Author: F.M., ctb: T.S.

BB.mod.stab.glm <- function(data, BB.data, s.model, model="linear", maxit.glm=25)
  {
    if (model == "binomial") {
      
      options(warn=-1) 
	  
      regmod <- glm(s.model, data=BB.data, family = binomial(link="logit"),
                    na.action=na.exclude, control=glm.control(maxit = maxit.glm))
					
      options(warn=0) 
	  
	  c.regmod <- glm(s.model, data=data, family = binomial(link="logit"),
                      na.action=na.exclude, control=glm.control(maxit = maxit.glm))
					  
    } else if (model == "linear") {
      regmod <- lm(s.model, data=BB.data, na.action=na.exclude)
      c.regmod <- lm(s.model, data=data, na.action=na.exclude)
    }
    c.namen <- names(c.regmod$coefficients)
    BB.namen <- names(regmod$coefficients)
    mislevpos <- !is.element(c.namen, BB.namen)
    if (any(mislevpos == T)) {
      help.coeff <- regmod$coefficients
      regmod$coefficients <- c.regmod$coefficients
      regmod$coefficients[mislevpos==T] <- 0
      regmod$coefficients[mislevpos==F] <- help.coeff
      regmod$xlevels <- c.regmod$xlevels
      regmod$rank <- c.regmod$rank
      regmod$assign <- c.regmod$assign
      regmod$qr$pivot <- c.regmod$qr$pivot
      regmod$qr$rank <- c.regmod$qr$rank
    }
    x <- list(model=regmod, c.model=c.regmod, mislevpos=mislevpos)
    return(x)
  }
