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

"xp.get.disp"<-
  function(gamdata,
           parnam,
           covnams,
           family="gaussian",
           ...)
{
  ##
  ## Run the null gam
  ##
  form <- as.formula(paste(parnam,"~1"))
  gam.null <- gam(form,data=gamdata)


  ## check that categorical covariates have more than one factor
  for(cov in covnams){
    if(is.factor(gamdata[, cov]) && nlevels(gamdata[,cov])==1){
      covnams=covnams[covnams!=cov]
    }
  }

  sel=rep(FALSE,times=length(covnams))
  names(sel)=covnams
  sel2=rep(FALSE,times=length(covnams))
  names(sel2)=covnams
  
  for(i in covnams) {
    ##
    ## Run the gam with the covariate entering linearly
    ##
    form <- as.formula(paste(parnam,"~",i))
    gam1 <- gam(form,data=gamdata)

    ##
    ## If we are dealing with a continous covariate
    ## run the gam on the non-linear function as well
    ## and comapre it to the linear fit
    if(!is.factor(gamdata[, i])) {
      form <- as.formula(paste(parnam," ~ ns(",i,", df=2)"))
      gam2 <- gam(form,data=gamdata)
      p    <- anova(gam1, gam2, test = "F")$Pr[2]
      
      ##
      ## If the non-linear was the better one, compare it to the null
      ## gam else compare the ninear to the null gam
      ##
      if(p < 0.05) {
	 p <- anova(gam.null, gam2, test = "F")$Pr[2]	 
	 if(p < 0.05){
	   sel2[i] <- T
       	 }
       } else {
         p <- anova(gam.null, gam1, test = "F")$Pr[2]
         if(p < 0.05){
	    sel[i] <- T
          }
       }
      ##
      ## If we are delaing with a factor, comapre it to the null gam
      ##
    } else {
      p <- anova(gam.null, gam1, test = "F")$Pr[2]
      if(p < 0.05){
        sel[i] <- T
      }
    }
  }
  
    
  ##
  ## Assemble the formula to use in the dispersion getting gam
  ##
  
  form <- NULL
  if(any(sel)) {
    form <- paste(names(sel)[sel], collapse = "+")
  }
  if(any(sel2)) {
    ncov <- names(sel2)[sel2]
    for(i in 1:length(ncov))
      ncov[i] <- paste("ns(",ncov[i],",df=2)",sep="")

    if(!any(is.null(form)))
      form <- paste(form,"+",paste(ncov,collapse="+"))
    else
      form <- paste(ncov,collapse="+")
  }

  if(is.null(form))
    gamform <- as.formula(paste(parnam,"~ 1"))
  else
    gamform <- as.formula(paste(parnam,"~",form))

  ##
  ## Run the dispersion getting GAM
  ##

  if(family == "gaussian") {
    gam3 <- gam(gamform,data=gamdata)
  } else {
    gam3 <- gam(gamform,data=gamdata,family=quasi(link=identity,variance="mu^2"))
  }

  disp <- summary(gam3)$dispersion

  ret.list <- list(covs = form,
		   formula = gamform,
		   dispersion = disp)
  return(ret.list)
}
