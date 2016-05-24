#####################################################################
## Individual-based re-implementation in R of:
##
## Rinke, K. & Petzoldt, T. (2003). Effects of temperature and food
## on individual growth and reproduction of Daphnia and their 
## consequences on the population level. Limnologica 33: p293-304.
##
#####################################################################


##===================================================================
## 1. Define the Model as a simecol object
##===================================================================

##-------------------------------------------------------------------
## 1.1 Derive a new class from the simecol base class simObj
##-------------------------------------------------------------------

setClass("indbasedModel",
         representation(
           parms  = "list",
           init   = "data.frame"
         ),
         contains = "simObj"
)

##-------------------------------------------------------------------
## 1.2 Define the model object
##-------------------------------------------------------------------

ibm_daphnia <- new("indbasedModel",
  main = function(time, init, parms) {
     init <- live(init, parms)
     init <- survive(init, parms)
     init <- hatch(init, parms)
     init
  },
  equations = list(
    newdaphnia = function(SON, n) {
      if (n>0) {
        data.frame(age = rep(0, n), size = SON, eggs = 0, eggage = 0)
      } else {
        NULL
      }
    },
    bottrell = function(temp) {
      exp(3.3956 + 0.2193 * log(temp) - 0.3414 * log(temp)^2)
    },
    tefi = function(time, temp, food, parms){
      with(parms, {
        deltaL <- L_0 - L_0_Hall
        k      <- b1 * exp(b2 * temp)
        L_max  <- (a1 * food)/(a2 + food) + a3 - k * a4
        L      <- L_max - (L_max - L_0_Hall) * exp (-k * time) + deltaL
        E      <- (X_max_slope * food)/(K_s_slope + food) * L +
                    beta_min * (1 - exp(-u_c * food))
        as.data.frame(cbind(L, E))
    })},
    live = function(inds, parms){
      with(parms,{
        ninds       <- nrow(inds)
        inds$age    <- inds$age + DELTAT
        inds$eggage <- ifelse(inds$size > SAM & inds$eggs > 0,
                              inds$eggage + DELTAT, 0)
        tefi_out    <- tefi(inds$age, temp, food, parms)
        inds$size   <- tefi_out$L
        neweggs     <- round(tefi_out$E)
        inds$eggs   <- ifelse(inds$size > SAM & inds$eggage==0,
                              neweggs, inds$eggs)
        inds
    })},
    survive  = function(inds, parms) subset(inds, inds$age < parms$maxage),
    hatch = function(inds, parms) {
      newinds <- NULL
      with(parms, {
        have.neo  <- 0
        new.neo   <- 0
        edt       <- bottrell(temp)
        have.neo  <- which(inds$eggs > 0 & inds$eggage > edt)
        eggs      <- inds$eggs[have.neo]
        new.neo   <- sum(eggs)
        inds$eggs[have.neo]   <- 0
        inds$eggage[have.neo] <- 0
        newinds <- newdaphnia(L_0, new.neo)
        rbind(inds, newinds)
      })
    }
  ),
  parms = list(
    # parameters of the somatic growth equation
    a1          = 1.167,    # (mm)
    a2          = 0.573,    # (mg L^-1)
    a3          = 1.420,    # (mm)
    a4          = 2.397,    # (d),
    b1          = 1.089e-2, # (d^-1)
    b2          = 0.122,    # ((deg. C)^-1)
    # parameters of the clutch size equation
    X_max_slope = 23.83,    # (eggs)
    K_s_slope   = 0.65,     # (mg L^-1)
    beta_min    = -29.28,   # (eggs)
    u_c         = 1,        # (L mg^-1) unit conversion factor
    # parameters of the individual-based model
    L_0_Hall    = 0.35,     # (mm) SON (size of neonanates) of the Hall data
    L_0         = 0.65,     # (mm) SON
    SAM         = 1.50,     # (mm) SAM (size at maturity)
    maxage      = 60,       # (d)
    # constant environmental conditions
    temp        = 20,       # (deg C)
    food        = 0.5       # (mg L^-1)
  ),
  init = data.frame(age=0, size=0.65, eggs=0, eggage=0),
  times = c(from=0, to=60, by=1),
  solver = "iteration"
)

##-------------------------------------------------------------------
## 1.3 Optionally define a user provided solver function that stores
##     the data more model specific than the default iteration
##-------------------------------------------------------------------

myiteration <- function(y, times = NULL, func = NULL, parms = NULL,
                        animate = FALSE, ...) {
  observer <- function(res) {
    # eggs, size, age, eggage
    number   <- nrow(res)
    meansize <- mean(res$size)
    meanage  <- mean(res$age)
    meaneggs <- mean(res$eggs)
    c(number=number, meansize=meansize, meanage=meanage, meaneggs=meaneggs)
  }
  init              <- y@init
  times             <- fromtoby(y@times)
  func              <- y@main
  parms             <- y@parms
  inputs            <- y@inputs
  equations         <- y@equations
  equations         <- addtoenv(equations)
  environment(func) <- environment()
  parms$DELTAT <- 0
  res <- observer(init)
  out <- res
  for (i in 2:length(times)) {
    time <- times[i]
    parms$DELTAT <- times[i] - times[i-1]
    init <- func(time, init, parms)
    res  <- observer(init)
    out  <- rbind(out, res)
   }
  row.names(out) <- NULL
  out <- cbind(times, out)
  as.data.frame(out)
}

##-------------------------------------------------------------------
## 1.4 Define a user specified plot function according to the
##     data structure stored by "myiteration"
##-------------------------------------------------------------------
setMethod("plot", c("indbasedModel", "missing"), function(x, y, ...) {
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  par(mfrow=c(2, 2))
  o <- out(x)
  plot(o$times, o$meanage,  type="l", xlab="Time", ylab="Mean age (d)")
  plot(o$times, o$meaneggs, type="l", xlab="Time", ylab="Eggs per individual")
  plot(o$times, o$number,   type="l", xlab="Time", ylab="Abundance")
  plot(o$times, o$number,   type="l", xlab="Time", ylab="Abundance", log="y")
})

##===================================================================
## 2. Simulate the Model
##===================================================================

solver(ibm_daphnia) <- "myiteration"
#ibm_daphnia <- sim(ibm_daphnia)
#plot(ibm_daphnia)
