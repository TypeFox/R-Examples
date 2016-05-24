#  R package CosmoPhotoz file R/plotDiagPhotoZ.R
#  Copyright (C) 2014 COIN
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License version 3 as published by
#the Free Software Foundation.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#' @title Plot diagnostics for photometric redshift  estimations
#'
#' @description \code{plotDiagPhotoZ} returns diagnostic plots from the results of 
#' photometric redshifts. Different types of plots are available: a density plot of the 
#' error distribution (\code{errordist}), a predicted versus observed contour plot 
#' (\code{predobs}), violin plot showing the error distribution at different redshift 
#' bins (\code{errorviolins}) and a box plot showing the errors at each different 
#' redshift bin (\code{box}). The produced plots are returned as ggplot2 objects.
#' 
#' @import ggplot2 ggthemes gridExtra
#' @param  photoz vector containing photoz data
#' @param  specz vector containing spectroscopic redshift data
#' @param  type a string with one of the following values: \code{errordist}, \code{predobs}, \code{errorviolins} or \code{box}
#' @param  npoints an integer indicating how many points should be used to create the \code{predobs} plot (if 0, all points will be used)
#' @return ggplot object 
#' @examples
#' # First, generate some mock data
#' ppo <- runif(1000, min=0.1, max=2)
#' ppo_ph <- rnorm(length(ppo), mean=ppo, sd=0.05)
#' 
#' # Then generate the plots
#' plotDiagPhotoZ(ppo_ph, ppo, type="errordist")
#' #plotDiagPhotoZ(ppo_ph, ppo, type="predobs")
#' #plotDiagPhotoZ(ppo_ph, ppo, type="errorviolins")
#' #plotDiagPhotoZ(ppo_ph, ppo, type="box")
#' 
#' @usage plotDiagPhotoZ(photoz, specz, type, npoints)
#' 
#' @author Rafael S. de Souza, Alberto Krone-Martins
#' 
#' @keywords hplot
#' @export
plotDiagPhotoZ <- function(photoz, specz, type=c("errordist", "predobs", "errorviolins", "box"), npoints=0) {
  
  # First some basic error control
  if( ! (type %in% c("errordist", "predobs", "errorviolins","box"))) {
    stop("Error in plotDiagPhotoZ :: the chosen plot type is not implemented.")
  } 
  if( ! is.vector(photoz) ) {
    stop("Error in plotDiagPhotoZ :: photoz is not a vector, and the code expects a vector.")
  }
  if( ! is.vector(specz) ) {
    stop("Error in plotDiagPhotoZ :: specz is not a vector, and the code expects a vector.")
  }
  
  # To prevent CRAN check notes regarding no visible bindings
  sigma=zspec=zphot=z_spec=z_photo=..level..=NULL

  # Now, for the real work
  # If the user wants to plot the error distributions
  if(type=="errordist") {
    sigm <- (photoz-specz)/(1+specz)
    sigmT <- sigm
    # Just to make sure that very big outliers will not be used for the density estimation
    sigm <- sigm[-which(abs(sigm) > abs(median(sigm) + 60*mad(sigm)))]
    # if the rejection was too strong, get back...
    if (length(sigm) == 0) {
      sigm <- sigmT    
    }
    sig <- data.frame(sigma=sigm)

    g1 <- ggplot(sig) + geom_density(aes(x=sigma), fill="#31a354", alpha=0.6) + coord_cartesian(c(-1, 1)) +
      xlab(expression((z[phot]-z[spec])/(1+z[spec]))) +
      theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans") +
      theme(plot.title = element_text(hjust=0.5),
          axis.title.y=element_text(vjust=0.75),
          axis.title.x=element_text(vjust=-0.25),
          text = element_text(size=20))
    return(g1)
  }
  
  # If the user wants to plot the predicted versus the reference values
  if(type=="predobs") {
    comb <- data.frame(zspec=specz, zphot=photoz)
    # Just to make sure that very big outliers will not be used for the density estimation
    pppp <- abs(comb$zphot - comb$zspec)
    comb <- comb[-which(pppp > abs(median(pppp) + 60*mad(pppp))),]
    # if the rejection was too strong, get back...
    if (length(comb$zspec)/length(pppp) < 0.95) {
      comb <- data.frame(zspec=specz, zphot=photoz)
    }
    
    if(npoints != 0) {
      idx <- sample(1:length(comb$zspec), replace=F, size = npoints)
      comb <- comb[idx,]
    }
    
    p1 <- ggplot(comb, aes(x=zspec, y=zphot))
    p2 <- p1 + stat_density2d(bins=200, geom="polygon", aes(fill =..level.., alpha=..level..), na.rm = TRUE, trans="log", n = 250,contour = TRUE) +
#      p2 <- p1 + stat_density2d(bins=200, geom="polygon", aes(), na.rm = TRUE, trans="log", n = 250,contour = TRUE) +
      coord_cartesian(xlim=c(0, max(specz)), ylim=c(0, max(specz)))+xlab(expression(z[spec]))+ylab(expression(z[phot])) +
      scale_fill_gradient2(guide="none", low = "#c7e9c0", mid="#41ab5d", high = "#00441b",space = "rgb") +
      geom_abline(intercept = 0) + theme(legend.text = element_text(colour="gray40"), legend.title=element_blank(), text = element_text(size=20), legend.position=c(0.1,0.75), axis.line = element_line(color = 'black')) +
      geom_density2d(colour="gray60", alpha=0.3, breaks = c(1, 5,10,25,50,100,200,250)) + theme_gdocs() +
      scale_alpha(guide="none")
#      scale_fill_gradient2(guide="none",low = "red",mid="cyan",high = "blue2", space = "Lab") 
    return(p2)
  }

  # If the user wants to plot the error distribution as violin plots within predetermined bins
  if(type=="errorviolins") {
    b2 <- factor(floor(specz * 10)/10)
    error_photoZ <- (specz-photoz)/(1+specz)
    dfd <- data.frame(z_photo=error_photoZ, z_spec=b2)  
    p <- ggplot(dfd) + xlab(expression(z[spec])) + ylab(expression((z[photo]-z[spec])/(1+z[spec]))) + ylim(-0.5, 0.5)
    p <- p + theme(legend.position = "none", axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))
    p <- p + geom_violin(aes(z_spec, z_photo), fill="#31a354", alpha=0.8)+
      theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans") +
      theme(plot.title = element_text(hjust=0.5),
            axis.title.y=element_text(vjust=0.75),
            axis.title.x=element_text(vjust=-0.25),
            text = element_text(size=20))  
    return(p)
    }

  # If the user wants to plot the error distribution as box plots within predetermined bins
  if(type=="box") {
    b2 <- factor(floor(specz * 10)/10)
    error_photoZ <- (specz-photoz)/(1+specz)
    dfd <- data.frame(z_photo=error_photoZ, z_spec=b2)  
    p <- ggplot(dfd) + xlab(expression(z[spec])) + ylab(expression((z[photo]-z[spec])/(1+z[spec]))) + ylim(-0.5, 0.5)
    p <- p + theme(legend.position = "none", axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))
    p <- p + geom_boxplot(aes(z_spec, z_photo), notch=F,fill="#31a354", alpha=0.8,outlier.colour = "gray")+
    theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans") +
    theme(plot.title = element_text(hjust=0.5),
          axis.title.y=element_text(vjust=0.75),
          axis.title.x=element_text(vjust=-0.25),
          text = element_text(size=20))  
    return(p)
  }
}
