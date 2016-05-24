# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

"xp.check.scope" <-
  function(object,
           covnam=xvardef("covariates", object),
           nmods=object@Prefs@Gam.prefs$nmods,
           smoother1=object@Prefs@Gam.prefs$smoother1,
           smoother2=object@Prefs@Gam.prefs$smoother2,
           smoother3=object@Prefs@Gam.prefs$smoother3,
           smoother4=object@Prefs@Gam.prefs$smoother4,
           arg1=object@Prefs@Gam.prefs$arg1,
           arg2=object@Prefs@Gam.prefs$arg2,
           arg3=object@Prefs@Gam.prefs$arg3,
           arg4=object@Prefs@Gam.prefs$arg4,
           excl1=object@Prefs@Gam.prefs$excl1,
           excl2=object@Prefs@Gam.prefs$excl2,
           excl3=object@Prefs@Gam.prefs$excl3,
           excl4=object@Prefs@Gam.prefs$excl4,
           extra=object@Prefs@Gam.prefs$extra,
           ...)
{
  scp <- xp.scope3(object,
                   covnam=covnam,
                   nmods=nmods,
                   smoother1=smoother1,
                   smoother2=smoother2,
                   smoother3=smoother3,
                   smoother4=smoother4,
                   arg1=arg1,
                   arg2=arg2,
                   arg3=arg3,
                   arg4=arg4,
                   excl1=excl1,
                   excl2=excl2,
                   excl3=excl3,
                   excl4=excl4,
                   extra=extra)
  return(scp)
}
