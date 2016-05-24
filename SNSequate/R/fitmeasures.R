### fitmeasures.R                   
### A collection of measures to investigate goodness of fit
###
### Copyright: Jorge Gonzalez, 2012.
### Last modification: 25-05-2012.
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or (at
### your option) any later version.
###
### This program is distributed in the hope that it will be useful, but
### WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
### General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
###
### Author contact information:
###
###      Jorge Gonzalez B.
###      Department of Statistics
###      Facultad de Matematicas
###      Pontificia Universidad Catolica de Chile
###      Casilla 306, Correo 22 
###      Santiago
###      Chile
###      Voice: +56-2-3545467  URL  : http://www.mat.puc.cl/~jgonzale
###      Fax  : +56-2-3547729  Email: jgonzale@mat.puc.cl
###

gof <- function(obs, fit, methods=c("FT"), p.out=TRUE){
  retvals <- list()
  cat("\n")
  
  if("FT" %in% methods){
    cat("Freeman-Tukey Residuals:\n")
    cat("------------------------\n")
    ftres <- ft.res(obs, fit, p.out=p.out)
    print(ftres)
    retvals <- c(retvals, list(ft.res=ftres))
    cat("\n")
  }
  
  if("KL" %in% methods){
    cat("Symmetrised Kullback-Leibler divergence:\n")
    cat("----------------------------------------\n")
    kl <- kl.divergence(obs, fit) + kl.divergence(fit, obs)
    print(kl)
    retvals <- c(retvals, list(kl.div=kl))
    cat("\n")
  }
  
  if("Chisq" %in% methods){
    cat("Pearson's Chi-squared Test:\n")
    cat("---------------------------\n")
    chisq <- chisq.test(obs, fit)
    chisq <- list(statistic=chisq$stat, df=chisq$parameter, p.value=chisq$p.value)
    cat(paste0("X-squared = ",chisq$stat, ", df = ",chisq$df, ", p-value = ",chisq$p.value))
    retvals <- c(retvals, list(chisq=chisq))
    cat("\n")
  }
  
  cat("---------------------------\n\n")
  
  return(retvals)
  
}


ft.res <- function(obs, fit, p.out=TRUE){
  FT <- as.numeric(sqrt(obs) + sqrt(obs+1) - sqrt(4*fit+1))
  
  if(p.out){
    par(mfrow=c(1,1)) # Reset the plot grid
    qqnorm(FT, main="Freeman-Tukey Residuals QQPlot")
    qqline(FT, col=2)
  }
  
  return(FT)
}


kl.divergence <- function(p, q){
  pn <- p / sum(p)
  qn <- q / sum(q)
  return( sum(pn*(log(pn) - log(qn)))  )
}
