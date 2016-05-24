# MS and VV 25-9-11
# Fitting a group of gamlss.family distributions
# TO DO
# add all available distributions
#-------------------------------------------------------------------------------
# the grouping of distributions
#-------------------------------------------------------------------------------
.realline <- c( 
  "NO", "GU", "RG" ,"LO", "NET", # 2 par
  "TF", "PE", "SN1", "SN2",  # 3 par
  "SHASHo", "EGB2", "JSU", "SEP1", "SEP2", "SEP3", "SEP4", # 4 par
   "ST1", "ST2", "ST3", "ST4", "ST5", "GT")  # 4 par 
#--------------------------------------------------------------------------------
.realplus <- c( "EXP", # 1 par
  "GA","IG","LOGNO", "WEI3", "IGAMMA","PARETO2", # 2 par
  "BCCGo", "exGAUS", "GG", "GIG",  # 3 par
  "BCTo", "BCPEo")  # 4 par 
#--------------------------------------------------------------------------------
.real0to1 <- c("BE", "BEINF", "BEINF0", "BEINF1", "BEOI", "BEZI", "GB1")
#--------------------------------------------------------------------------------
.realAll <- c( "EXP", # 1 par
  "GA","IG","LOGNO", "WEI3", "IGAMMA","PARETO2", # 2 par
  "BCCGo", "exGAUS", "GG", "GIG",  # 3 par
  "BCTo", "BCPEo",
  "GU", "RG" ,"LO", "NET", # 2 par
    "TF", "PE", "SN1", "SN2",    # 3 par
  "SHASHo", "EGB2", "JSU", "SEP1", "SEP2", "SEP3", "SEP4", # 4 par
   "ST1", "ST2", "ST3", "ST4", "ST5", "GT"             
              )  # 4 par
#--------------------------------------------------------------------------------               
#--------------------------------------------------------------------------------
.counts <- c("PO", "GEOM", "LG", "YULE", "WARING", "NBI", "NBII", "PIG", "DEL", "SICHEL", "ZIP","ZIP2", "ZAP", "ZALG", "ZINBI", "ZANBI", "ZIPIG")
#--------------------------------------------------------------------------------
.binom <- c("BI", "BB", "ZIBI", "ZIBB", "ZABI", "ZABB" )
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
fitDist <- function(y,
                    k = 2, # for the AIC
                 type = c("realAll", "realline", "realplus","real0to1","counts", "binom" ), 
           try.gamlss = FALSE,  # whether to try the gamlss() if gamlssML() fails
                extra = NULL,
                 data = NULL, ...) # for extra distributions to include 
{
  # if (!is.null(data)) {attach(data); on.exit(detach(data))}
  #if (!is.null(data)) {attach(data, name="TheDatA"); on.exit(detach(TheDatA))}
  y <- if (!is.null(data)) get(deparse(substitute(y)), envir=as.environment(data)) else y
  type <- match.arg(type)
    DIST <- switch(type, "realAll"=.realAll, 
                        "realline"=.realline, 
                        "realplus"=.realplus,
                         "real0to1"=.real0to1,
                          "counts"=.counts,
                           "binom"=.binom 
                  )
if  (!is.null(extra)) DIST <- c(DIST, extra)
    m0 <- switch(type,  "realAll"= gamlssML(y, family=NO),
                       "realline"= gamlssML(y, family=NO), 
                       "realplus"= gamlssML(y, family=EXP),
                       "real0to1"= gamlssML(y, family=BE),
                         "counts"= gamlssML(y, family=PO),
                          "binom"= gamlssML(y, family=BI) 
                 ) 
  failed <- list() 
    fits <- list()
for (i in 1:length(DIST)) 
{
    m1 <- try(gamlssML(y,family=DIST[i], ...), silent=TRUE)
        if (any(class(m1)%in%"try-error")&&try.gamlss==TRUE) 
        { 
         m1 <-  try(gamlss(y~1,family=DIST[i], trace=FALSE, ...),  silent=TRUE)
        }
    
       if (any(class(m1)%in%"try-error"))
       {
            failed <- c(failed, DIST[i]) 
       }
       else
      {
        aic <- AIC(m1, k=k)
        names(aic) <- DIST[i]
        fits <- c(fits, aic)
        if (AIC(m1, k=k) < AIC(m0, k=k)) 
        {
          m0<-m1 
        }
      }
}
m0$failed <- failed
fits <- unlist(fits)
m0$fits <- fits[order(fits)]         
m0  
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
