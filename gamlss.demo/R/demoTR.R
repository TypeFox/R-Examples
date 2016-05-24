
demo.NOtr <- function()
{if (interactive())     
      {   
        mu <- 0.5 ; sigma <- 0.5 ;  left <- -10 ; right <- 10
       # require(gamlss.tr)
      plotf <- function(panel) 
         {   
              with(panel, {
         mu   <- as.numeric(panel$mu)
      sigma   <- as.numeric(panel$sigma)
         nu   <- as.numeric(panel$nu)
         left <- as.numeric(panel$left)
        right <- as.numeric(panel$right)
         dNOTR <- trun.d(c(left,right),"NO", type="both")
         xgrid <- seq(-6, 6, length=200)
         txgrid <- seq(left, right, length=200) 
         dgrid <- dNO(xgrid,  mu=mu, sigma=sigma)
         tdgrid <- dNOTR(txgrid, mu=mu, sigma=sigma)
         maxy<- max(c(dgrid, tdgrid))
         plot(xgrid, dgrid, type = "l", col = "blue", ylab="pdf", xlab="x", ylim=c(0, maxy), main="Normal truncated distribution" )
        lines(txgrid, tdgrid, col="red") 
        lines(rep(left,30), seq(0, maxy, length=30), col="grey" )
        lines(rep(right,30), seq(0, maxy, length=30), col="grey")
        points(left,  dNOTR(left,  mu=mu,sigma=sigma), col="red")
        points(right, dNOTR(right, mu=mu,sigma=sigma), col="red")
               })
      panel
         }
            }
   bipanel <- rp.control("NO truncated distribution", mu=0, sigma=1, left=-3L, right=3L)
   rp.slider(bipanel, variable = mu, from=-5, to=5, resolution=0.001,  action =  plotf, title="mu",  showvalue = T) 
   rp.slider(bipanel, variable = sigma, from=0.01, to=5, resolution=0.001,  action =  plotf, title="sigma",  showvalue = T)  
   rp.doublebutton(bipanel, variable = left,  action = plotf, initval = -3, step = .1, range = c(-5, 5), showvalue = T, "left point")
   rp.doublebutton(bipanel, variable = right, action = plotf, initval =  3, step = .1, range = c(-5, 5), showvalue = T, "right point")
   
   # rp.slider(bipanel, variable = left, from=-5, to=5, resolution=0.1,  action =  plotf, title="left point",  showvalue = T)  
   #rp.slider(bipanel, variable = right, from=-5, to=5, resolution=0.1,  action =  plotf, title="right point",  showvalue = T)  
   
   #rp.textentry(bipanel, left, plotf, labels = c("left"), initval = c(-5))
   #rp.textentry(bipanel, right, plotf, labels = c("right"), initval = c(5))
   #rp.do(bipanel, plotf)
        #rp.textentry(tables.panel, xobs, tables.redraw, "Observed value",
         #           pos = c(140, 75, 160, 25))
}


demo.GAtr <- function()
{if (interactive())     
      {   
        mu <- 5 ; sigma <- 0.5 ;  left <- 0 ; right <- 10
       # require(gamlss.tr)
      plotf <- function(panel) 
         {   
              with(panel, {
         mu   <- as.numeric(panel$mu)
      sigma   <- as.numeric(panel$sigma)
         nu   <- as.numeric(panel$nu)
         left <- as.numeric(panel$left)
        right <- as.numeric(panel$right)
         dGATR <- trun.d(c(left,right),"GA", type="both")
         xgrid <- seq(0.001, 15, length=200)
         txgrid <- seq(left+0.001, right, length=200) 
         dgrid <- dGA(xgrid,  mu=mu, sigma=sigma)
         tdgrid <- dGATR(txgrid, mu=mu, sigma=sigma)
         maxy<- max(c(dgrid, tdgrid))
         plot(xgrid, dgrid, type = "l", col = "blue", ylab="pdf", xlab="x", ylim=c(0, maxy), main="Gamma truncated distribution" )
        lines(txgrid, tdgrid, col="red") 
        lines(rep(left,30), seq(0, maxy, length=30), col="grey" )
        lines(rep(right,30), seq(0, maxy, length=30), col="grey")
        points(left,  dGATR(left,  mu=mu,sigma=sigma), col="red")
        points(right, dGATR(right, mu=mu,sigma=sigma), col="red")
               })
      panel
         }
            }
   bipanel <- rp.control("NO truncated distribution", mu=5, sigma=.5, left=0, right=10L)
   rp.slider(bipanel, variable = mu,    from=0.02, to=10, resolution=0.01,  action =  plotf, title="mu",  showvalue = T) 
   rp.slider(bipanel, variable = sigma, from=0.02, to=5, resolution=0.01,  action =  plotf, title="sigma",  showvalue = T)  
   rp.doublebutton(bipanel, variable = left,  action = plotf, initval = 0,  step = .1, range = c(0, 5), showvalue = T, "left point")
   rp.doublebutton(bipanel, variable = right, action = plotf, initval = 10, step = .1, range = c(5, 15), showvalue = T, "right point")
   #rp.slider(bipanel, variable = left, from=0, to=5, resolution=0.01,  action =  plotf, title="left point",  showvalue = T)  
   #rp.slider(bipanel, variable = right, from=5, to=15, resolution=0.01,  action =  plotf, title="right point",  showvalue = T)  
   
   #rp.textentry(bipanel, left, plotf, labels = c("left"), initval = c(-5))
   #rp.textentry(bipanel, right, plotf, labels = c("right"), initval = c(5))
   #rp.do(bipanel, plotf)
        #rp.textentry(tables.panel, xobs, tables.redraw, "Observed value",
         #           pos = c(140, 75, 160, 25))
}
