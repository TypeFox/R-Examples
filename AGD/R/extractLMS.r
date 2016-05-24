# extractLMS, ageGrid



#'Extracts LMS values from a gamlss object.
#'
#'Extract LMS values from a gamlss object for solutions that transform the age
#'axis according to the M-curve.
#'
#'It is crucial that \code{t.age} in \code{data} correspond to exactly the same
#'age transformation as used to fit the \code{gamlss} object. Age grid values
#'beyond the range of \code{data$age} produce \code{NA} in the L, M and S
#'values. Parameter \code{flatAge} should be one of the values of the age grid.
#'
#'@param fit A gamlss object containing the final fit on transformed age,
#'\code{t.age}.
#'@param data A data frame containing the original data, with both \code{age}
#'and \code{t.age}
#'@param sex A character vector indicating whether the fit applied to males
#'\code{sex="M"} or females \code{sex="F"}. The default is \code{sex="M"}.
#'@param grid A character vector indicating the desired age grid. See
#'\code{ageGrid()} for possible options. The default is a
#'\code{grid="classic"}, a grid of 59 age points.
#'@param decimals A numerical vector of length 3 indicating the number of
#'significant digits for rounding of the L, M and S curves, respectively.
#'@param flatAge A scalar indicating the age beyond which the L, M and S values
#'should be constant. The default (NULL) is not to flatten the curves.
#'@return A data frame with rows corresponding to time points, and with the
#'following columns: \code{sex},\code{x},\code{L},\code{M},\code{S}.
#'@author Stef van Buuren, 2010
#'@keywords distribution
#'@examples
#'
#'\dontrun{
#'#
#'library(gamlss)
#'boys <- boys7482
#'
#'# calculate initial M curve
#'data <- na.omit(boys[,1:2])
#'f0154  <- gamlss(hgt~cs(age,df=15,c.spar=c(-1.5,2.5)),
#'                 sigma.formula=~cs(age,df=4,c.spar=c(-1.5,2.5)),
#'                 data=data,family=NO,
#'                 control=gamlss.control(n.cyc=3))                      
#'
#'# calculate transformed age
#'t.age <- fitted(lm(data$age~fitted(f0154)))
#'t.age <- t.age - min(t.age)
#'data.t <- data.frame(data,t.age=t.age)
#'
#'# calculate final solution
#'f0106r <- gamlss(hgt~cs(t.age,df=10,c.spar=c(-1.5,2.5)),
#'                 sigma.formula=~cs(t.age,df=6,c.spar=c(-1.5,2.5)),
#'                 data=data.t,family=NO,
#'                 control=gamlss.control(n.cyc=3))
#'
#'# extract the LMS reference table in the 'classic' age grid
#'nl4.hgt.boys <- extractLMS(fit = f0106r, data=data.t, grid="compact", 
#'                dec = c(0,2,5))
#'nl4.hgt.boys
#'
#'
#'# flatten the reference beyond age 20Y (not very useful in this data)
#'nl4.hgt.boys.flat <- extractLMS(fit = f0106r, data=data.t, flatAge=20)
#'nl4.hgt.boys.flat
#'
#'# use log age transformation
#'data.t <- data.frame(data, t.age = log(data$age))
#'f0106rlog <- gamlss(hgt~cs(t.age,df=10,c.spar=c(-1.5,2.5)),
#'                 sigma.formula=~cs(t.age,df=6,c.spar=c(-1.5,2.5)),
#'                 data=data.t,family=NO,
#'                 control=gamlss.control(n.cyc=1))
#'
#'nl4.hgt.boys.log <- extractLMS(fit = f0106rlog, data=data.t)
#'nl4.hgt.boys.log
#'}
extractLMS <- function(fit, data, sex="M", grid="classic",
            decimals = c(4,4,4), 
            flatAge = NULL)
{
  # Extracts the LMS table after the 'Cole-transformation'
  # or any other transformation in t.age
  # fit  final gamlss object
  # data should contain both age and t.age
  
  check.names(df=data, needed=c("age","t.age"))
  if (!is.gamlss(fit)) stop("fit not a gamlss object.")
  tm <- data$t.age[which.min(data$age)]
  # if (abs(tm)>0.0001) stop("wrong offset of transformed age")
  # if (min(data$t.age) < 0) warning("Negative transformed age found. Results are unpredictable.")

  grd <- ageGrid(grid)
  grid.age <- grd$year

  minage <- min(data$age, na.rm=TRUE)
  maxage <- max(data$age, na.rm=TRUE)
  outside <- grid.age < minage | grid.age > maxage
  grid.age <- grid.age[!outside]

  t.grid.age <- approx(x=data$age, y=data$t.age, xout=grid.age, ties=mean)$y
  if (length(t.grid.age)==0) stop("No overlap between age grid and data.")
  newdata <- data.frame(t.age=t.grid.age)

  # lms <- predictAll(fit, newdata=newdata, data=data)
  # print(lms)
  lms <- predictAll(fit, newdata=newdata, data=data)
  lms$mu    <- round(lms$mu, decimals[2])
  lms$sigma <- round(lms$sigma/lms$mu, decimals[3])  ## THIS IS PROBABLY AN ERROR!! 25/07/2013
  if (length(lms)>2) lms$nu <- round(lms$nu, decimals[1])
  else lms$nu <- 1
  
  lms <- as.data.frame(lms)
  result <- data.frame(sex=sex, x=grd$year, L=NA, M=NA, S=NA)
  result[!outside, "L"] <- lms$nu
  result[!outside, "M"] <- lms$mu
  result[!outside, "S"] <- lms$sigma
  
  if (!is.null(flatAge)) {
    if (!(flatAge %in% result$x)) stop("FlatAge value (', FlatAge,') not found in age grid")
    flatLMS <- result[flatAge==result$x, c("L","M","S")]
    result[!outside & flatAge<result$x, c("L","M","S")] <- flatLMS
  }
  
  return(result)
}



#'Creates an age grid according to a specified format.
#'
#'Creates an age grid according to a specified format.
#'
#'
#'@param grid A character string specifying one of the following:
#'\code{"compact"}, \code{"classic"}, \code{"extensive"}, \code{"0-104w"},
#'\code{"0-24m"}, \code{"0-21y"}, \code{"0-21yd"} or \code{"0-21yc"}.  The
#'default is \code{"compact"}, which produces an age grid between 0 and 21
#'years with 95 points.
#'@return A list with five components: \code{format}, \code{year},
#'\code{month}, \code{week} and \code{day} containing the age grid in different
#'units.
#'@author Stef van Buuren, 2010
#'@keywords distribution
#'@examples
#'
#'
#'age <- ageGrid("classic")$year
#'
#'
ageGrid <- function(grid="compact"){
  formats <- c("compact", "classic", "extensive", "0-104w", "0-24m", "0-21y", "0-21yd", "0-21yc")
  fmi <- pmatch(grid, formats)
  if (is.na(fmi)) stop("Grid format ",grid," unknown.")
  grid <- switch(fmi,
    compact =
      c((0:14)/365.25,
        (3:13)*7/365.25,
        seq(3,11.5,0.5)/12,
        (12:23)/12,
        seq(2, 21, 0.5)
      ),
    classic =
    {
      grid.weeks <- c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,
                      28,32,36,40,44,48,52,56,60,64)
      c(grid.weeks*7/365.25,seq(1.5,21,0.5))
    },
    extensive = (0:(365.25*21+1))/365.25,
    week = (0:104)*7/365.25,
    month = (0:24)/12,
    year = 0:21,
    dyear = seq(0, 21, 0.1),
    cyear = seq(0, 21, 0.01)
  )
  
  year    <- round(grid, 4)
  month   <- round(grid*12, 4)
  week    <- round(grid*365.25/7, 4)
  day     <- round(grid*365.25, 4)
  return(list(format = formats[fmi],
   year = year, month = month, week = week, day = day))
}


