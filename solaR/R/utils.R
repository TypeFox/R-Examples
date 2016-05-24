###Indices temporales

###Quizás sería útil que todas estas funciones admitiesen un zoo directamente
hour <- function(x) 
{
    as.numeric(format(x, "%H"))
}

minute <- function(x)
{
    as.numeric(format(x, "%M"))
}

second <- function(x) 
{
    as.numeric(format(x, "%S"))
}

hms <- function(x)
{
    hour(x)+minute(x)/60+second(x)/3600
}

doy <- function(x){
  as.numeric(format(x, '%j'))
}

dom <- function(x){
  as.numeric(format(x, '%d'))
}

month <- function(x){
  as.numeric(format(x, '%m'))
}

year <- function(x){
  as.numeric(format(x, '%Y'))

}

DoY <- function(x){format(x, '%j')}

DoM <- function(x){format(x, '%d')}

Month <- function(x){format(x, '%m')}

Year <- function(x){format(x, '%Y')}

dst <- function(x)                      #Adelanto horario por verano
   {
     as.POSIXlt(x)$isdst
   }

##Angulos

d2r <- function(x){x*pi/180}

r2d <- function(x){x*180/pi}

h2r <- function(x){x*pi/12}
h2d <- function(x){x*180/12}

r2h <- function(x){x*12/pi}
d2h <- function(x){x*12/180}

r2sec <- function(x){x*12/pi*3600}


##Trunca un POSIXct a días
truncDay  <-  function(x){as.POSIXct(trunc(x, units='days'))}

## Check if daily indexes are equal (used in fCompD and fTemp)
checkIndexD <- function(ix, iy)
{
    dx <- truncDay(ix)
    dy <- truncDay(iy)
    test <- all.equal(dx, dy,  check.attributes = FALSE)
    if (!isTRUE(test))
        stop('daily indexes do not match.')
}


###Husos horarios
lonHH<-function(tz)
    {            #Calcula la longitud (en radianes) de un huso horario
      stopifnot(class(tz)=='character')
      tHH <- as.POSIXct('2000-1-1 12:00:00', tz=tz)
      tUTC <- as.POSIXct(format(tHH, tz='UTC'), tz=tz)
      h2r(as.numeric(tHH-tUTC))
    }

  
local2Solar <- function(x, lon=NULL){	
  tz=attr(x, 'tzone')
  if (tz=='' || is.null(tz)) {tz='UTC'}
  ##Adelanto oficial por verano
  AO=3600*dst(x)
  AOneg=(AO<0)
  if (any(AOneg)) {
    AO[AOneg]=0
    warning('Some Daylight Savings Time unknown. Set to zero.')
  }
  ##Diferencia entre la longitud del lugar y la longitud del huso horario LH
  LH=lonHH(tz)
  if (is.null(lon)) 
    {deltaL=0
   } else
  {deltaL=d2r(lon)-LH
 }
  ##Hora local corregida en UTC
  ##    tt <- format(x-AO+r2sec(deltaL), tz=tz)
  tt <- format(x, tz=tz)
  result <- as.POSIXct(tt, tz='UTC')-AO+r2sec(deltaL)
  ##      result <- as.POSIXct(tt, tz='UTC')
  result
}


##cbind garantizando conservación del index (para tz='UTC', principalmente)
CBIND <- function(..., index=NULL){
  args <- list(...)
  cdata <- lapply(args, coredata)
  result0 <- as.data.frame(cdata)
  if (is.null(index)){
    return(zoo(result0, index(args[[1]])))
  } else {
    return(zoo(result0, index))
  }
}

##Convierte un difftime en un número de horas
diff2Hours  <- function(by){
  if (!inherits(by, 'difftime')) {
    stop('This function is only useful for difftime objects.')
  } else {
    return(as.numeric(by, units='hours'))
  }
}

char2diff <- function(by){
  if (!is.character(by)) {
    stop('This function is only useful for character strings.')
  } else {
    ##Adaptado de seq.POSIXt
    by2 <- strsplit(by, " ", fixed = TRUE)[[1L]]
    if (length(by2) > 2L || length(by2) < 1L) 
      stop("invalid 'by' string")
    units <- c("secs", "mins", "hours")
    valid <- pmatch(by2[length(by2)], units)
    if (is.na(valid)) {
      stop("invalid string for 'by'")
    } else {
      unitValid <- units[valid]
      if (length(by2)==1) {
        by2=1
      } else {
        by2=as.numeric(by2[1])
      }
      result <- as.difftime(by2,units=unitValid)
      return(result)
    }
  }
}

sample2Hours <- function(by){
  if (is.character(by)) {
    y <- char2diff(by)
    return(diff2Hours(y))
  } else if (inherits(by, 'difftime')) {
    return(diff2Hours(by))
  } else {stop('by must be a character or difftime.')}
}
  
P2E <- function(x, by){
  Nm=1/sample2Hours(by)
  sum(x, na.rm=1)/Nm                    #Potencia a Energía
} 
###OJO: no exportadas
solvePac <- function(x, Cinv){
  Vdc=x[1]
  PdcN=x[2]
  V <- c(1, Vdc, Vdc^2)
  Ki=t(colSums(V*t(Cinv)))
  A=Ki[3]
  B=Ki[2]+1
  C=Ki[1]-(PdcN)
  result <- (-B+sqrt(B^2-4*A*C))/(2*A)
  result
}

dailySum <- function(x, by){##x is a time series
  if (missing(by)) {by=DeltaT(x)}
  res <- aggregate(x, by=truncDay, FUN=P2E, by)
  return(res)
  }

monthlySum <- function(x, by){##x is a INTRADAILY time series
  if (missing(by)) {by=DeltaT(x)}
  res <- aggregate(x, by=as.yearmon, FUN=P2E, by)
  return(res)
  }

yearlySum <- function(x, by){##x is a INTRADAILY time series
  if (missing(by)) {by=DeltaT(x)}
  res <- aggregate(x, by=year, FUN=P2E, by)
  return(res)
  }

dailyMean <- function(x){##x is a time series
  res <- aggregate(x, by=doy, FUN=mean, na.rm=1)
  return(res)
  }

monthlyMean <- function(x){##x is a time series
  res <- aggregate(x, by=as.yearmon, FUN=mean, na.rm=1)
  return(res)
  }

yearlyMean <- function(x){##x is a time series
  res <- aggregate(x, by=year, FUN=mean, na.rm=1)
  return(res)
  }


##No exportada
DeltaT <- function(x){
  return(median(diff(index(x))))###spend a long time with large series, ¿mean(x, 0.2)?
  }


factorI<-function(x, index.rep, breaks=3, ...){
  ##x es una variable extraida con $ de un slot de un objeto
  ##index.rep es el índice que relaciona las variables diarias con las instantátneas.
  ##Se obtiene con indexRep(object)
  var.fac<-cut(x, breaks=breaks, ...)
  ## indexI.day<-as.POSIXct(trunc(indexI, 'day'))
  ## mtch<-match(indexI.day, indexD, nomatch=0)
  result<-var.fac[index.rep]
}



  
