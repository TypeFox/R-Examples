#-------------------------------------------------------------------------------------
# SHASHo
demo.SHASHo <- function()
{ if (interactive())  
   {   mu <- sigma <- nu <- tau <- NULL
      density.draw <- function(panel) 
         {
           y <- seq(-10, 10, length=200)
         plot(y, dSHASHo(x=y, mu=panel$mu, sigma=panel$sigma, nu=panel$nu, tau=panel$tau), ylab="f(y)", main="SHASHo probability density functions", type="l")
          panel   
         }
         SHASHpanel <- rp.control('SHASHo family', sigma = 1, mu=0, nu=0, tau=1)
        rp.slider(SHASHpanel, variable=mu, from=-4, to=4, resolution=0.01,  action = density.draw, title="mu",  showvalue = TRUE)  
        rp.slider(SHASHpanel, variable=sigma, from=0.1, to=5, resolution=0.01,  action = density.draw, title="sigma",  showvalue = TRUE) 
        rp.slider(SHASHpanel, variable=nu,  from=-4, to=4,  resolution=0.01,  action = density.draw, title="nu",  showvalue = TRUE)
        rp.slider(SHASHpanel, variable=tau, from=0.01, to=2,  resolution=0.01,  action = density.draw, title="tau",  showvalue = TRUE)  
   }
}
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# SHASHo2
demo.SHASHo2 <- function()
{ if (interactive())  
   {   mu <- sigma <- nu <- tau <- NULL
      density.draw <- function(panel) 
         {
           y <- seq(-10, 10, length=200)
         plot(y, dSHASHo2(x=y, mu=panel$mu, sigma=panel$sigma, nu=panel$nu, tau=panel$tau), ylab="f(y)", main="SHASHo2 probability density functions", type="l")
          panel   
         }
         SHASHpanel <- rp.control('SHASHo2 family', sigma = 1, mu=0, nu=0, tau=1)
        rp.slider(SHASHpanel, variable=mu, from=-4, to=4, resolution=0.01,  action = density.draw, title="mu",  showvalue = TRUE)  
        rp.slider(SHASHpanel, variable=sigma, from=0.1, to=5, resolution=0.01,  action = density.draw, title="sigma",  showvalue = TRUE) 
        rp.slider(SHASHpanel, variable=nu,  from=-3, to=3,  resolution=0.01,  action = density.draw, title="nu",  showvalue = TRUE)
        rp.slider(SHASHpanel, variable=tau, from=0.01, to=4,  resolution=0.01,  action = density.draw, title="tau",  showvalue = TRUE)  
   }
}
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# count data
# PARETO2
demo.PARETO2 <- function()
{ if (interactive())  
   {   mu <- sigma <- NULL
      density.draw <- function(panel) 
         {
           y <- seq(0.01, 10, length=200)
         plot(y, dPARETO2(x=y, mu=panel$mu, sigma=panel$sigma), ylab="f(y)", main="PARETO2 probability density functions", type="l")
          panel   
         }
         PARETOpanel <- rp.control('PARETO2 family', sigma = 1, mu=1)
        rp.slider(PARETOpanel, variable=mu, from=0.1, to=10, resolution=0.01,  action = density.draw, title="mu",  showvalue = TRUE)  
        rp.slider(PARETOpanel, variable=sigma, from=0.1, to=10, resolution=0.01,  action = density.draw, title="sigma",  showvalue = TRUE) 
   }
}
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# count data
# PARETO2
demo.PARETO2o <- function()
{ if (interactive())  
   {   mu <- sigma <- NULL
      density.draw <- function(panel) 
         {
           y <- seq(0.01, 10, length=200)
         plot(y, dPARETO2o(x=y, mu=panel$mu, sigma=panel$sigma), ylab="f(y)", main="PARETO2o probability density functions", type="l")
          panel   
         }
         PARETOpanel <- rp.control('PARETO2o family', sigma = 1, mu=1)
        rp.slider(PARETOpanel, variable=mu, from=0.1, to=10, resolution=0.01,  action = density.draw, title="mu",  showvalue = TRUE)  
        rp.slider(PARETOpanel, variable=sigma, from=0.1, to=10, resolution=0.01,  action = density.draw, title="sigma",  showvalue = TRUE) 
   }
}
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# count data
# IGAMMA
demo.IGAMMA <- function()
{ if (interactive())  
   {   mu <- sigma <- NULL
      density.draw <- function(panel) 
         {
           y <- seq(0.01, 50, length=200)
         plot(y, dIGAMMA(x=y, mu=panel$mu, sigma=panel$sigma), ylab="f(y)", main="IGAMMA probability density functions", type="l")
          panel   
         }
         IGAMMApanel <- rp.control('IGAMMA family', sigma = 1, mu=1)
        rp.slider(IGAMMApanel, variable=mu, from=0.1, to=10, resolution=0.01,  action = density.draw, title="mu",  showvalue = T)  
        rp.slider(IGAMMApanel, variable=sigma, from=0.1, to=10, resolution=0.01,  action = density.draw, title="sigma",  showvalue = T) 
   }
}
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# count data
# YULE
demo.YULE <- function()
{ if (interactive())  
   {   mu <- NULL
      density.draw <- function(panel) 
         {
         op<- par(mfrow=c(2,1))	
           x <- seq(0, 30, 1)
         plot(x, dYULE(x, mu=panel$mu), type="h", ylab="f(y)", main="probability density functions")
         points(x, dYULE(x, mu=panel$mu), col="blue")
        plot(x, pYULE(x, mu=panel$mu), type="h", ylab="F(y)", main="cumulative probability functions")
        par(op) 
          panel   
         }
         YULEpanel <- rp.control('YULE',  mu=1, ymax=1)
        rp.slider(YULEpanel, variable=mu, from=0.1, to=50, resolution=0.01,  action = density.draw, title="mu",  showvalue = T)  
   }
}
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# count data
# WARING
demo.WARING <- function()
{ if (interactive())  
   {   mu <- NULL
      density.draw <- function(panel) 
         {
         op<- par(mfrow=c(2,1))	
           x <- seq(0, 30, 1)
         plot(x, dWARING(x, mu=panel$mu), type="h", ylab="f(y)", main="probability density functions")
         points(x, dWARING(x, mu=panel$mu), col="blue")
        plot(x, pWARING(x, mu=panel$mu), type="h", ylab="F(y)", main="cumulative probability functions")
        par(op) 
          panel   
         }
         WARINGpanel <- rp.control('WARING',  mu=1, ymax=1)
        rp.slider(WARINGpanel, variable=mu, from=0.1, to=50, resolution=0.01,  action = density.draw, title="mu",  showvalue = T)  
   }
}
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# count data
# GEOM
demo.GEOM <- function()
{ if (interactive())  
   {   mu <- NULL
      density.draw <- function(panel) 
         {
         op<- par(mfrow=c(2,1))	
           x <- seq(0, 30, 1)
         plot(x, dGEOM(x, mu=panel$mu), type="h", ylab="f(y)", main="probability density functions")
         points(x, dGEOM(x, mu=panel$mu), col="blue")
        plot(x, pGEOM(x, mu=panel$mu), type="h", ylab="F(y)", main="cumulative probability functions")
        par(op) 
          panel   
         }
         GEOMpanel <- rp.control('GEOM',  mu=1, ymax=1)
        rp.slider(GEOMpanel, variable=mu, from=0.1, to=50, resolution=0.01,  action = density.draw, title="mu",  showvalue = T)  
   }
}
