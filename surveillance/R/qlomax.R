################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Simple implementation of the quantile function of the Lomax distribution
### (we could also use VGAM::qlomax, but this would be slightly slower)
###
### Copyright (C) 2012-2013 Sebastian Meyer
### $Revision: 691 $
### $Date: 2013-12-18 14:09:41 +0100 (Mit, 18 Dez 2013) $
################################################################################

qlomax <- function (p, scale, shape)
    scale * ((1-p)^(-1/shape) - 1)

