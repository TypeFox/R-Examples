#
#  Copyright (C) 2004-2012 Friedrich Leisch and Bettina Gruen
#  $Id: z.R 4834 2012-08-02 10:17:09Z gruen $
#

###**********************************************************
## Backward compatibility

## component model driver
FLXglm <- FLXMRglm
FLXglmFix <- FLXMRglmfix
FLXmclust <- FLXMCmvnorm
FLXbclust <- FLXMCmvbinary

## concomitant model driver
FLXmultinom <- FLXPmultinom
FLXconstant <- FLXPconstant
