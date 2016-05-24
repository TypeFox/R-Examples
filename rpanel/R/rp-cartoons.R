rp.cartoons <- function(hscale = 1) {

panel.launch <- function(menu.panel) {
   if (menu.panel$demo == "q-q plot") {
      bc.fn <- function(y, lambda) {
         if (abs(lambda) < 0.001) z <- log(y)
         else z <- (y^lambda - 1)/ lambda
      }
      qq.draw <- function(panel) {
         z <- bc.fn(panel$y, panel$lambda)
         qqnorm(z, main = paste("lambda =", round(panel$lambda, 2)))
         panel
      }
      panel <- rp.control(y = exp(rnorm(50)), lambda = 1)
      rp.slider(panel, lambda, -2, 2, qq.draw)
      rp.do(panel, qq.draw)
   }
   else if (menu.panel$demo == "bubbleplot")
      rp.bubbleplot(log(gdp), log(co2.emissions), 1960:2007, size = population, 
         col = life.expectancy, interpolate = TRUE, hscale = hscale)
   else if (menu.panel$demo == "Binomial distribution") {
      plot.binomial <- function(panel) {
         with(panel, {
            n <- as.numeric(n)
            probs <- dbinom(0:n, n, prob)
            plot(c(0,n), c(0,1), type = "n", xlab = "x", ylab = "Probability")
            segments(0:n, rep(0, n+1), 0:n, probs)
            title(paste("Binomial:  n =", n, "  p =", round(prob, 3)))
         })
         panel
      }
      rp.binomial <- function() {
         pname <- rp.control("Binomial probabilities", n = 20, prob = 0.5)
         rp.slider(pname, prob, 0, 1, initval = 0.5, title = "Binomial proby, p:", 
                   action = plot.binomial)
         rp.textentry(pname, n, plot.binomial, "Sample size, n:")
         rp.do(pname, plot.binomial)
      }
      rp.binomial()
   }
   else if (menu.panel$demo == "Tables") {
      rp.tables(hscale = hscale)
   }
   else if (menu.panel$demo == "Normal fitting") {
      y <- rnorm(50, mean = 10, sd = 0.5)
      rp.normal(y, hscale = hscale)
      rm(y)
   }
   else if (menu.panel$demo == "Confidence intervals") {
      rp.ci(hscale = hscale)
   }
   else if (menu.panel$demo == "Regression - CofE (Attend)") {
      with(CofE, {
         rp.regression(Attend, Giving, hscale = hscale)
      })
   }
   else if (menu.panel$demo == "Regression - rodent") {
      with(rodent, {
         rp.regression(log(Mass), log(Speed), hscale = hscale)
      })
   }
   else if (menu.panel$demo == "Regression - CofE (Attend & Employ)") {
      with(CofE, {
         rp.regression(cbind(Employ, Attend), Giving, hscale = hscale)
      })
   }
   else if (menu.panel$demo == "Anova - two-way") {
      with(poisons, {
         rp.anova(1/stime, treatment, poison, hscale = hscale)
      })
   }
   else if (menu.panel$demo == "Ancova") {
      with(gullweight, {
         rp.ancova(hab, weight, month, hscale = hscale)
      })
   }
   else if (menu.panel$demo == "Logistic regression") {
      with(river, {
         rp.logistic(Temperature, Low, hscale = hscale)
      })
   }
   else if (menu.panel$demo == "Repeated measurements") {
      LH    <- luthor[,2:16]
      gp    <- factor(luthor[,1])
      times <- c(1:5,(5+(1:10)/2))
      rp.rmplot(log(LH), fac = gp, timept = times, hscale = hscale)
   }
   else if (menu.panel$demo == "Likelihood - exponential") {
      rp.likelihood("sum(log(dexp(data, theta)))", aircond, 0.005, 0.03, hscale = hscale)
   }
   else if (menu.panel$demo == "Likelihood - gamma") {
      rp.likelihood("sum(log(dgamma(data, theta[1], theta[2])))",
                 aircond, c(0.3, 0.005), c(3, 0.06))
   }
   else if (menu.panel$demo == "Power") {
      rp.power(hscale = hscale)
   }
   else if (menu.panel$demo == "Density estimation (1d)") {
      sm.density(tephra$Al2O3, panel = TRUE)
   }
   else if (menu.panel$demo == "Density estimation (2d)") {
      with(airpc, {
         y <- cbind(Comp.1, Comp.2)[Period==3,]
         sm.density(y, panel = TRUE)
      })
   }
   else if (menu.panel$demo == "Flexible regression (1d)") {
      with(trawl, {
         sm.regression(Longitude, Score1, panel = TRUE)
      })
   }
   else if (menu.panel$demo == "Flexible regression (2d)") {
      with(trawl, {
         Position <- cbind(Longitude - 143, Latitude)
         sm.regression(Position, Score1, panel = TRUE)
      })
   }
   else if (menu.panel$demo == "Surface uncertainty") {
      with(trawl, {
         location  <- cbind(Longitude, Latitude)
         model     <- sm.regression(location, Score1, ngrid = 15, display = "none")
         longitude <- model$eval.points[ , 1]
         latitude  <- model$eval.points[ , 2]
         xgrid     <- as.matrix(expand.grid(longitude, latitude))
         S         <- sm.weight2(location, xgrid, model$h)
         covar     <- tcrossprod(S) * model$sigma^2
         rp.surface(model$estimate, covar, longitude, latitude, location, Score1, hscale = hscale)
      })
   }
   else if (menu.panel$demo == "Gulls") {
      rp.gulls()
   }
   else if (menu.panel$demo == "Quakes") {
      with(quakes, {
         rp.plot4d(cbind(long, lat), depth, mag, hscale = hscale)
      })
   }
   else if (menu.panel$demo == "Clyde DO") {
      with(Clyde, {
         rp.plot4d(cbind(Doy, DO), Station, hscale = hscale)
      })
   }
   else if (menu.panel$demo == "SO2") {
      with(SO2, {
         location <- cbind(longitude, latitude)
         if (require(mgcv) & require(maps)) {
            location1 <- location[,1]
            location2 <- location[,2]
            model   <- gam(logSO2 ~ s(location1, location2, year))
            loc1    <- seq(min(location1), max(location1), length = 30)
            loc2    <- seq(min(location2), max(location2), length = 30)
            yr      <- seq(min(year), max(year), length = 30)
            newdata <- expand.grid(loc1, loc2, yr)
            names(newdata) <- c("location1", "location2", "year")
            model <- predict(model, newdata)
            model <- list(x = cbind(loc1, loc2), z = yr,
                          y = array(model, dim = rep(30, 3)))
            mapxy <- map('world', plot = FALSE,
                         xlim = range(longitude), ylim = range(latitude))
            rp.plot4d(location, year, logSO2, model, hscale = hscale,
                      col.palette = rev(heat.colors(20)),
                      foreground.plot = function() map(mapxy, add = TRUE))
         }
         else
            rp.plot4d(location, year, logSO2, model, col.palette = rev(heat.colors(20)),
                      hscale = hscale)
      })
   }
   else if (menu.panel$demo == "Spatial simulation") {
      rp.geosim(hscale = hscale)
   }
   else if (menu.panel$demo == "Mururoa") {
      rp.mururoa(hscale = hscale)
   }
   else if (menu.panel$demo == "Firth") {
      rp.firth(hscale = hscale)
   }
   menu.panel
}

menu.panel <- rp.control("Cartoons", homer = FALSE, number.list = list(),
      ss = list(), trans = list(), theta = list())
menu.list  <-  list(list("Introductory",
                         "q-q plot",
                         "bubbleplot",
                         "Binomial distribution",
                         "Tables",
                         "Normal fitting",
                         "Confidence intervals"
                         ),
                    list("Regression",
                         "Regression - CofE (Attend)",
                         "Regression - rodent",
                         "Regression - CofE (Attend & Employ)",
                         "Anova - two-way",
                         "Ancova",
                         "Logistic regression"
                         ),
                    list("Advanced",
                         "Repeated measurements",
                         "Likelihood - exponential",
                         "Likelihood - gamma",
                         "Power"
                          ),
                    list("Smoothing",
                         "Density estimation (1d)",
                         "Density estimation (2d)",
                         "Flexible regression (1d)",
                         "Flexible regression (2d)",
                         "Surface uncertainty"
                          ),
                    list("Applications",
                         "Gulls",
                         "Quakes",
                         "Clyde DO",
                         "SO2",
                         "Spatial simulation",
                         "Mururoa",
                         "Firth"
                          )
                    )
                    
if (!require(sm)) menu.list <- menu.list[-4]
rp.menu(menu.panel, demo, menu.list, action = panel.launch)
image.file <- file.path(system.file(package = "rpanel"), "images", "cartoons.gif")
rp.image(menu.panel, image.file)

invisible()
}

