################################################################################
# "Working with dynamic models for agriculture"
# R script for pratical work
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2012-05-11
# Model described in the book, Appendix. Models used as illustrative examples: description and R code
################################ FUNCTIONS #####################################
#' @title The Epirice model (Disease model for rice)
#' @description \strong{Model description.} Adapted from Savary et al.(2012)
#' @param param : a vector of parameters
#' @param weather : weather data.frame for one single year
#' @param sdate : date to begin simulation (day of year) (default 1)
#' @param ldate : date to end simulation (day of year) (default 120)
#' @param H0 : initial number of plant's healthy sites(default 600)
#' @return data.frame with daily state variable
#' @seealso \code{\link{epirice.multi.simule}}
#' @export
epirice.model<-function(param, weather,sdate=1,ldate=120, H0=600)
{
# Parameter values of the model, read from the param vector
Site.size=param["Site.size"]
Sx=param["Sx"]
RRG=param["RRG"]
RRS=param["RRS"]
EODate=param["EODate"]
p=param["p"]
i=param["i"]
rl=param["rl"]
RcOpt=param["RcOpt"]
RcT=param["RcT"]
RcW=param["RcW"]
RcA=param["RcA"]
a=param["a"]
k=param["k"]
SenescType=param["SenescType"]

# Correction coefficient for infection - age (epirice)
fRcA<-function(age){graphRcA=matrix(c(0.00, 1.00, 5.00, 1.00, 10.0, 1.00, 15.0, 0.9, 20.0, 0.8, 25.0, 0.7,
 30.0, 0.64, 35.0, 0.59, 40.0, 0.53,45.0, 0.43, 50.0, 0.32, 55.0, 0.22,
 60.0, 0.16, 65.0, 0.09, 70.0, 0.03, 75.0, 0.02, 80.0, 0.02, 85.0, 0.02,
 90.0, 0.01, 95.0, 0.01, 100, 0.01, 105, 0.01, 110, 0.01, 115, 0.01, 120, 0.01),
 ncol=2, byrow=TRUE,  dimnames = list(NULL,c("day", "RcA")))
return(approx(graphRcA[,"day"], graphRcA[,"RcA"], age, method="linear")$y)}
# Correction coefficient for infection - temperature (epirice)
fRcT<-function(T){graphRcT=matrix(c(0.0, 0.0, 10.0, 0.00, 15.0, 0.5, 20.0, 1.00, 25.0, 0.6, 30.0, 0.2, 35.0,
0.05, 40.0, 0.01, 45.0, 0.00, 60, 0), ncol=2, byrow=TRUE,  dimnames = list(NULL,c("T", "RcT")))
return(approx(graphRcT[,"T"], graphRcT[,"RcT"], T, method="linear")$y)
}
# Correction coefficient for infection - Relative humidity and rain (epirice)
fRcW<-function(RH,RAIN){if(RH>=90 | RAIN>=5) {return(1)} else {return(0)}   }

    # Initialize variables
    # states variables, as vectors initialized to NA
    #H		Number of healthy sites	NSites
    H<-rep(NA,ldate)
    #H		Number of senesced sites	NSites
    S<-rep(NA,ldate)
    #L		Number of latent sites	NSites
    L<-rep(NA,ldate)
    #I		Number of infectious sites	NSites
    II<-rep(NA,ldate)
    #P		Number of post-infectious (removed) sites	NSites
    P<-rep(NA,ldate)
    #TS		Number of total sites	NSites
    TS<-rep(NA,ldate)
    #TOTDIS	Number of total diseased sites	NSites
    TOTDIS<-rep(NA,ldate)

    # Initialize state variables when sowing on day "sdate"
    H[sdate]<- H0
    S[sdate]<- 0
    L[sdate]<- 0
    II[sdate]<- 0
    P[sdate]<- 0
    TOTDIS[sdate]<- L[sdate] + II[sdate] + P[sdate]
    TS[sdate]<- H[sdate]+TOTDIS[sdate]

    # Simulation loop
    for (day in sdate:(ldate-1))
        {
        CC= H[day]/TS[day]
        Rc = RcOpt *fRcA(day-sdate)*fRcT((weather$TMIN[day]+weather$TMAX[day])/2)*fRcW(weather$RH2M[day],weather$RAIN[day])
        
        if (day==(sdate+EODate)) {Starter=1} else {Starter=0}
        # Calculate rates of change of state variables
        RG = RRG*H[day]*(1 - (TS[day]/Sx))
        # We propose a simplification of the cohort formalism for RTransfer and RRemoval
        RTransfer = L[day]/p
        RRemoval = II[day]/i
 
        RInf = Rc*II[day]*CC^a + Starter
        RSenesced = RRS*H[day] + RRemoval*SenescType
        
        dH = RG - RInf - RSenesced
        dS = RSenesced
        dL = RInf - RTransfer
        dII = RTransfer - RRemoval
        dP = RRemoval
        dTOTDIS = dL + dII + dP
        dTS = dH + dTOTDIS
        
        # Update state variables
        H[day+1]= H[day] + dH
        S[day+1]= S[day] + dS
        L[day+1]= L[day] + dL
        II[day+1]= II[day] + dII
        P[day+1]= P[day] + dP
        TOTDIS[day+1]= TOTDIS[day] + dTOTDIS
        TS[day+1]= TS[day] + dTS
        #PRI[day+1]=PRI[day]+dPRI
        }
        severity=(TOTDIS-P)/(TS-P) * 100
        # End simulation loop
    return(data.frame(day=sdate:ldate,DACE=((sdate:ldate)-sdate),H=H[sdate:ldate],L=L[sdate:ldate],II=II[sdate:ldate],P=P[sdate:ldate],TS=TS[sdate:ldate],TOTDIS=TOTDIS[sdate:ldate],severity=severity[sdate:ldate]))
}
################################################################################
#' @title Wrapper function to run Epirice multiple times (for multiple sets of inputs)
#' @param param : a vector of parameters
#' @param multi.simule : matrix of n row definition of input variable : site, year and date of transplantation.
#' @param all : if you want a matrix combining multi.simule and output (default = FALSE)
#' @return matrix with AUDPC for each input vector
#' @seealso \code{\link{epirice.model}}
#' @export
epirice.multi.simule <- function(param,multi.simule, all=FALSE){
# output : AUDPC only
AUDPC <- apply(multi.simule,1,function(v) sum(epirice.model(param, epirice.weather(working.year=v["year"], working.site=v["idsite"]),sdate=v["date_transplantation"],ldate=v["date_transplantation"]+120,H0=600)$severity, na.rm=TRUE))
print(AUDPC)
if(all) AUDPC = cbind(multi.simule,AUDPC = AUDPC)
return(AUDPC)
}

################################################################################
#' @title Define values of the parameters for the Epirice model
#' @return matrix with parameter values (nominal, binf, bsup)
#' @export
epirice.define.param <- function()
{
# nominal, binf, bsup
param<-data.frame(
Site.size= c(45,NA,NA),
Sx=c(30000,NA,NA),
RRG=c(0.1,NA,NA),
RRS=c(0.01,NA,NA),
EODate=c(15,NA,NA),
p=c(5,NA,NA),
i=c(20,NA,NA),
rl=c(0.28,NA,NA),
RcOpt=c(1.14,NA,NA),
RcT=c(25,NA,NA),
RcW=c(1,NA,NA),
RcA=c(1,NA,NA),
a=c(1,NA,NA),
k=c(0.025,NA,NA),
SenescType=c(1,NA,NA))
row.names(param)<-c("nominal","binf","bsup")
return(as.matrix(param))
}
################################################################################
#' @title Read weather data for Epirice (southern Asia weather)
#' @param working.year : year for the subset of weather data (default=NA : all the year)
#' @param working.site : site for the subset of weather data (default=NA : all the site)
#' @return data.frame with daily weather data for one or several site(s) and for one or several year(s)
#' @export
# Reading Weather data function
epirice.weather=function(working.year=NA, working.site=NA)
    {
    #day month year R Tmax Tmin rain ETP
    # idsite","GPSlatitude","GPSlongitude","WEYR","WEDAY","TMAX","TMIN","RAIN","RH2M"
    # Tmax : maximum temperature (celsius)
    # Tmin : minimum temperature (celsius)
    weather=weather_SouthAsia
    names(weather)[names(weather)=="WEYR"|names(weather)=="WEDAY"]= c("year","day")
    # if argument working.year/working.site is specified, work on one particular year/site
    if (!is.na(working.year)&!is.na(working.site)) {weather=weather[(weather$year==working.year)&(weather$idsite==working.site),] }
    else{
      if (!is.na(working.year)) {weather=weather[(weather$year==working.year),]}
      if (!is.na(working.site)) {weather=weather[(weather$idsite==working.site),]}}
    return (weather)
    }

# End of file
