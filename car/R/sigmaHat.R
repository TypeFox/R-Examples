#-------------------------------------------------------------------------------
# Revision history:
# 2010-07-01: moved from alr3 and renamed.  S. Weisberg
## method to return sigmaHat

sigmaHat <- function(object){UseMethod("sigmaHat")}
sigmaHat.default <- function(object){summary(object)$sigma}
sigmaHat.lm <- function(object) sigmaHat.default(object)
sigmaHat.glm <- function(object){sqrt(summary(object)$dispersion)}