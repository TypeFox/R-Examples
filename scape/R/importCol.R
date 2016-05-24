importCol <- function(res.file, Dev=FALSE, CPUE=FALSE, Survey=FALSE, CAc=FALSE, CAs=FALSE, CLc=FALSE, CLs=FALSE,
                      LA=FALSE, quiet=TRUE, ...)
{
  ## Implementation notes
  ## Generic read* functions: Vector, Matrix
  ## Specific get* functions: N, B, Sel, Dev, CPUE, Survey, CAc, CAs, CLc, CLs, LA
  ## The only global objects (used inside get* functions) are *.file, *.vector, and quiet

  ## 1  Define functions
  readVector <- function(keyword, same.line=TRUE, file=res.file, vector=res.vector)
  {
    ## Extract white-space delimited numeric vector followed by keyword
    line <- match(keyword, substring(vector,1,nchar(keyword)))
    v <- if(same.line)
      as.numeric(scan(file, what="", skip=line-1, nlines=1, quiet=TRUE)[-1])
    else
      as.numeric(scan(file, what="", skip=line, nlines=1, quiet=TRUE))
    if(!quiet) cat("vector...")
    return(v)
  }

  readMatrix <- function(keyword, nrow, header=FALSE, stripe=c("no","left","right","upper","lower"), file=res.file,
                         vector=res.vector)
  {
    ## Extract white-space delimited numeric matrix followed by keyword
    stripe <- match.arg(stripe)  # for striped matrices, like N@A in Coleraine: Female|Male|Female|Male|...
    line <- match(keyword,substring(vector,1,nchar(keyword))) + as.numeric(header)
    m <- scan(file, skip=line, nlines=nrow, quiet=TRUE)
    m <- matrix(m, byrow=TRUE, nrow=nrow)
    m <- switch(stripe, left=m[,seq(1,ncol(m)/2)], right=m[,seq(ncol(m)/2+1,ncol(m))], upper=m[seq(1,nrow(m)-1,by=2),],
                lower=m[seq(2,nrow(m),by=2),], m)
    if(!quiet) cat("matrix...")
    return(m)
  }

  getN <- function(sexes, years, ages)
  {
    ## Look for "\"Numbers_at_age_by_Year,sex_and_age\" in colera31, or "Numbers_at_age_by_Year,sex_and_age" in colera32
    if(!quiet) cat("N         ")
    nsexes <- length(sexes)
    nyears <- length(years)
    nages <- length(ages)
    if(nsexes == 1)
    {
      Nu <- readMatrix("Numbers_at_age_by_Year,sex_and_age", nrow=nyears*nsexes)
      N <- data.frame(Sex=rep(sexes,nyears*nages), Year=rep(years,each=nages), Age=rep(ages,nyears), N=as.vector(t(Nu)),
                      ...)
    }
    if(nsexes == 2)
    {
      Nf <- readMatrix("Numbers_at_age_by_Year,sex_and_age", nrow=nyears*nsexes, stripe="upper")
      Nm <- readMatrix("Numbers_at_age_by_Year,sex_and_age", nrow=nyears*nsexes, stripe="lower")
      N <- data.frame(Sex=rep(sexes,each=nyears*nages), Year=rep(rep(years,each=nages),2), Age=rep(ages,2*nyears),
                      N=as.vector(t(rbind(Nf,Nm))), ...)
    }
    if(!quiet) cat("OK\n")
    return(N)
  }

  getB <- function(years, gears)
  {
    ngears <- length(gears)
    if(!quiet) cat("B         ")
    vb <- readMatrix("Vulnerable_Biomass_by_Method_and_Year", nrow=ngears)
    sb <- readVector("Spawning_Biomass_by_Year", same.line=FALSE)
    y  <- c(readVector("Total_Catch_by_Method_and_Year", same.line=FALSE), NA)
    B <- data.frame(years=years, vb=t(vb), sb=sb, y=y, ...)
    names(B) <- if(ngears==1) c("Year","VB","SB","Y") else c("Year",paste("VB",gears,sep="."),"SB","Y")
    if(!quiet) cat("OK\n")
    return(B)
  }

  getSel <- function(gears, surveys, years, sexes, ages)
  {
    if(!quiet) cat("Sel       ")
    ngears <- length(gears)
    nsurveys <- length(surveys)
    nyears <- length(years)
    nsexes <- length(sexes)
    nages <- length(ages)
    com <- readMatrix("Commercial_age-specific_selectivity_by_method,Year,sex_and_age", nrow=ngears*nyears*nsexes)
    com <- com[seq(1, to=ngears*nyears*nsexes, by=nyears),]                # first year: selectivity
    srv <- readMatrix("Survey_age-specific_selectivity_by_survey,Year,sex_and_age", nrow=nsurveys*nsexes)
    fecundity <- readVector("Fecundity_by_year_and_age", same.line=FALSE)  # fecundity
    weight <- readVector("Weight_by_year,sex_and_age", same.line=FALSE)    # weight
    mat <- rep(ifelse(weight>0,fecundity/weight,0), nsexes)  # sexes have same maturity
    if(is.numeric(gears))
      gears <- paste("Gear", gears)
    if(is.numeric(surveys))
      surveys <- paste("Survey", surveys)
    Sel <- data.frame(Series=
                      c(rep(gears,each=nsexes*nages),rep(surveys,each=nsexes*nages),rep("Maturity",nsexes*nages)),
                      Sex=rep(rep(sexes,each=nages),ngears+nsurveys+1), Age=rep(ages,(ngears+nsurveys+1)*nsexes),
                      P=c(t(com),t(srv),mat), ...)
    if(!quiet) cat("OK\n")
    return(Sel)
  }

  getDev <- function(ages, years)
  {
    if(!quiet) cat("Dev       ")
    Dev <- list()
    Dev$sigmaR <- c(readVector("p_log_InitialDev",same.line=TRUE)[6], readVector("p_log_RecDev",same.line=TRUE)[6])
    names(Dev$sigmaR) <- c("Initial", "Annual")
    Dev$Initial <- readVector("log_InitialDev", same.line=TRUE)
    names(Dev$Initial) <- ages[-c(1,length(ages))]  # exclude first and last age
    Dev$Annual <- readVector("log_RecDev", same.line=TRUE)
    names(Dev$Annual) <- years[-length(years)]      # exclude last year
    if(!quiet) cat("OK\n")
    return(Dev)
  }

  getCPUE <- function(gears, years)
  {
    if(!quiet) cat("CPUE      ")
    nseries <- readVector("NCPUEindex")
    ngears <- length(gears)
    nyears <- length(years)
    obs <- readMatrix("indexmethodyearvaluecv", nrow=readVector("Number_of_CPUE_data",same.line=FALSE))
    obs <- data.frame(Series=obs[,1], Gear=obs[,2], Year=obs[,3], Obs=obs[,4], CV=obs[,5], ...)
    fit <- readMatrix("CPUE_Index_Trajectories", nrow=nseries)
    fit <- data.frame(Series=rep(1:nseries,each=nyears), Year=rep(years,nseries), Fit=as.vector(t(fit)), ...)
    CPUE <- merge(obs[,names(obs)!="Gear"], fit, all=TRUE)  # merge without looking at gear
    sgkey <- unique(obs[,c("Series","Gear")])
    CPUE <- merge(sgkey, CPUE)  # add gear column
    CPUE <- data.frame(Series=paste("Series ",CPUE$Series,"-",CPUE$Gear,sep=""), Year=as.integer(CPUE$Year),
                       Obs=CPUE$Obs, CV=CPUE$CV, Fit=CPUE$Fit, ...)
    if(!quiet) cat("OK\n")
    return(CPUE)
  }

  getSurvey <- function(years)
  {
    if(!quiet) cat("Survey    ")
    nyears <- length(years)
    nseries <- readVector("Nsurveyindex")
    obs <- readMatrix("indexyearvaluecv", nrow=readVector("Number_of_survey_data",same.line=FALSE))
    obs <- data.frame(Series=obs[,1], Year=obs[,2], Obs=obs[,3], CV=obs[,4], ...)
    fit <- readMatrix("Survey_Index_Trajectories", nrow=nseries)
    fit <- data.frame(Series=rep(1:nseries,each=nyears), Year=rep(years,nseries), Fit=as.vector(t(fit)), ...)
    Survey <- merge(obs, fit, all=TRUE)
    Survey$Series <- as.integer(Survey$Series)
    Survey$Year <- as.integer(Survey$Year)
    if(!quiet) cat("OK\n")
    return(Survey)
  }

  getCAc <- function(sexes, ages)
  {
    if(!quiet) cat("CAc       ")
    nsexes <- length(sexes)
    nages <- length(ages)
    nobs <- readVector("Number_of_Commercial_C@A", same.line=FALSE)
    obs <- readMatrix("methodyearsamplesizesex1a1sex1a2sex1a3", nrow=nobs)                     # "Observed_C@A"  !unique
    fit <- readMatrix("methodyearsamplesizesex1a1sex1a2sex1a3", nrow=nobs, header=2*(nobs+1))  # "Predicted_C@A" !unique
    CAc <- data.frame(Series=rep(obs[,1],each=nsexes*nages), Year=rep(obs[,2],each=nsexes*nages),
                      SS=rep(obs[,3],each=nsexes*nages), Sex=rep(rep(sexes,each=nages),nobs),
                      Age=rep(ages,nsexes*nobs), Obs=as.vector(t(obs[,-(1:3)])), Fit=as.vector(t(fit)), ...)
    CAc$Series <- as.integer(CAc$Series)
    CAc$Year <- as.integer(CAc$Year)
    CAc$Age <- as.integer(CAc$Age)
    if(!quiet) cat("OK\n")
    return(CAc)
  }

  getCAs <- function(sexes, ages)
  {
    if(!quiet) cat("CAs       ")
    nsexes <- length(sexes)
    nages <- length(ages)
    nobs <- readVector("Number_of_survey_C@A",same.line=FALSE)
    obs <- readMatrix("surveyyearsamplesizesex1a1sex1a2sex1a3", nrow=nobs)                     # "Observed_C@A"  !unique
    fit <- readMatrix("surveyyearsamplesizesex1a1sex1a2sex1a3", nrow=nobs, header=2*(nobs+1))  # "Predicted_C@A" !unique
    CAs <- data.frame(Series=rep(obs[,1],each=nsexes*nages), Year=rep(obs[,2],each=nsexes*nages),
                      SS=rep(obs[,3],each=nsexes*nages), Sex=rep(rep(sexes,each=nages),nobs),
                      Age=rep(ages,nsexes*nobs), Obs=as.vector(t(obs[,-(1:3)])), Fit=as.vector(t(fit)), ...)
    CAs$Series <- as.integer(CAs$Series)
    CAs$Year <- as.integer(CAs$Year)
    CAs$Age <- as.integer(CAs$Age)
    if(!quiet) cat("OK\n")
    return(CAs)
  }

  getCLc <- function(sexes, lengths)
  {
    if(!quiet) cat("CLc       ")
    nsexes <- length(sexes)
    nlengths <- length(lengths)
    nobs <- readVector("Number_of_Commercial_C@L", same.line=FALSE)
    obs <- readMatrix("methodyearsamplesizesex1l1sex1l2sex1l3", nrow=nobs)                 # "Observed_C@L"  !unique
    fit <- readMatrix("methodyearsamplesizesex1l1sex1l2sex1l3", nrow=nobs, header=nobs+1)  # "Predicted_C@L" !unique
    CLc <- data.frame(Series=rep(obs[,1],each=nsexes*nlengths), Year=rep(obs[,2],each=nsexes*nlengths),
                      SS=rep(obs[,3],each=nsexes*nlengths), Sex=rep(rep(sexes,each=nlengths),nobs),
                      Length=rep(lengths,nsexes*nobs), Obs=as.vector(t(obs[,-(1:3)])), Fit=as.vector(t(fit)), ...)
    CLc$Series <- as.integer(CLc$Series)
    CLc$Year <- as.integer(CLc$Year)
    CLc$Length <- as.integer(CLc$Length)
    if(!quiet) cat("OK\n")
    return(CLc)
  }

  getCLs <- function(sexes, lengths)
  {
    if(!quiet) cat("CLs       ")
    nsexes <- length(sexes)
    nlengths <- length(lengths)
    nobs <- readVector("Number_of_surveyC@L",same.line=FALSE)
    obs <- readMatrix("surveyyearsamplesizesex1l1sex1l2sex1l3", nrow=nobs)                     # "Observed_C@L"  !unique
    fit <- readMatrix("surveyyearsamplesizesex1l1sex1l2sex1l3", nrow=nobs, header=2*(nobs+1))  # "Predicted_C@L" !unique
    CLs <- data.frame(Series=rep(obs[,1],each=nsexes*nlengths), Year=rep(obs[,2],each=nsexes*nlengths),
                      SS=rep(obs[,3],each=nsexes*nlengths), Sex=rep(rep(sexes,each=nlengths),nobs),
                      Length=rep(lengths,nsexes*nobs), Obs=as.vector(t(obs[,-(1:3)])), Fit=as.vector(t(fit)), ...)
    CLs$Series <- as.integer(CLs$Series)
    CLs$Year <- as.integer(CLs$Year)
    CLs$Length <- as.integer(CLs$Length)
    if(!quiet) cat("OK\n")
    return(CLs)
  }

  getLA <- function(sexes, ages)
  {
    ## Sex | Age | Obs | Fit | CV
    if(!quiet) cat("LA        ")
    nsexes <- length(sexes)
    nages <- length(ages)
    nobs <- readVector("#femalesmales", same.line=FALSE, file=latage.file, vector=latage.vector)  # two elements
    obs <- readMatrix("VonBertalanfy--Lenght-at-agefit--Likelihood", nrow=sum(nobs), header=8)
    obs <- data.frame(Sex=rep(sexes,nobs), Age=obs[,1], Obs=obs[,2], ...)
    owarn <- options(warn=-1)                                               #\
    Linf <- readVector("VonBeratalanfy:Linf")[-(1:3)]                       # \
    K <- readVector("VonBeratalanfy:k")[-(1:3)]                             #  \  workaround for unusual input file
    t0 <- readVector("VonBeratalanfy:to")[-(1:3)]                           #   ) the five warnings are suppressed
    CV1 <- readVector("cvoftheFitbysex")[-(1:5)]                            #  /
    CVratio <- readVector("ratioofcv(L_an)/cv(L_a1)oftheFitbysex")[-(1:7)]  # /
    options(owarn)                                                          #/
    sigmaLA <- readVector("#LinearrelationshipofsigmaL@A(1=age;2=length---ignoreifW@Aissupplied)", same.line=FALSE,
                          file=txt.file, vector=txt.vector)[1]
    max.age <- c(max(obs$Age[obs$Sex==sexes[1]]), max(obs$Age[obs$Sex==sexes[2]]))  # max ages in LA data
    fit <- data.frame(Sex=rep(sexes,max.age), Age=c(1:max.age[1],1:max.age[2]), ...)
    fit$Fit[fit$Sex==sexes[1]] <- Linf[1]*(1-exp(-K[1]*(fit$Age[fit$Sex==sexes[1]]-t0[1])))
    fit$Fit[fit$Sex==sexes[2]] <- Linf[2]*(1-exp(-K[2]*(fit$Age[fit$Sex==sexes[2]]-t0[2])))
    if(sigmaLA == 1)  # CV ~ age
    {
      A <- rep(max(ages), 2)
      a <- cbind(fit$Age[fit$Sex==sexes[1]], fit$Age[fit$Sex==sexes[2]])
      fit$CV[fit$Sex==sexes[1]] <- CV1[1] + CV1[1]*(CVratio[1]-1)/(A[1]-1)*(a[,1]-1)
      fit$CV[fit$Sex==sexes[2]] <- CV1[2] + CV1[2]*(CVratio[2]-1)/(A[2]-1)*(a[,2]-1)
    }
    if(sigmaLA == 2)  # CV ~ length
    {
      L1 <- Linf*(1-exp(-K*(1-t0)))
      Ln <- Linf*(1-exp(-K*(max(ages)-t0)))
      fit$CV[fit$Sex==sexes[1]] <- CV1[1] + CV1[1]*(CVratio[1]-1)/(Ln[1]-L1[1])*(fit$Fit[fit$Sex==sexes[1]]-L1[1])
      fit$CV[fit$Sex==sexes[2]] <- CV1[2] + CV1[2]*(CVratio[2]-1)/(Ln[2]-L1[2])*(fit$Fit[fit$Sex==sexes[2]]-L1[2])
    }
    LA <- merge(obs, fit, by=c("Sex","Age"), all=TRUE)
    LA$Age <- as.integer(LA$Age)
    LA$Fit <- LA$Fit
    LA$CV <- LA$CV
    if(!quiet) cat("OK\n")
    return(LA)
  }

  ## 2  Parse args
  if(!file.exists(res.file)) stop("File ", res.file, " not found; use / or \\\\ separators")

  ## 3  Read dimensions
  res.vector <- readLines(res.file)                                    # string vector, one element being one line
  res.vector <- gsub("\"","", gsub("\t","",gsub(" ","",res.vector)))   # remove white space and quotes
  if(!quiet) cat("\nParsing text file ", res.file, ":\n\nPreamble  ", sep="")
  sexes <- if(readVector("Nsexes")==1) "Unisex" else c("Female","Male")
  gears <- seq(1, length.out=readVector("Nmethods"))        # possibly a vector of length zero
  surveys <- seq(1, length.out=readVector("Nsurveyindex"))  # "
  years <- seq(from=readVector("StartYear"), to=readVector("EndYear")+1)
  ages <- seq(from=1, to=readVector("Nages"))
  lengths <- seq(from=readVector("First_length"), by=readVector("Length_class_increment"),
                 length.out=readVector("Number_of_length_classes"))
  if(!quiet) cat("OK\n")

  ## 4  Read N, B, R, Sel
  model <- list()
  model$N   <- getN(sexes, years, ages)
  model$B   <- getB(years, gears)            # Recruits:
  rec       <- model$N[model$N$Age==1,]      #  by year and sex
  rec       <- tapply(rec$N, rec$Year, sum)  #  combine sexes
  model$B$R <- c(rec[-1], NA)                #  align with spawning stock
  model$Sel <- getSel(gears, surveys, years, sexes, ages)

  ## 5  Read Dev, CPUE, Survey, CAc, CAs, CLc, CLs, LA
  if(Dev)    model$Dev    <- getDev(ages, years)
  if(CPUE)   model$CPUE   <- getCPUE(gears, years)
  if(Survey) model$Survey <- getSurvey(years)
  if(CAc)    model$CAc    <- getCAc(sexes, ages)
  if(CAs)    model$CAs    <- getCAs(sexes, ages)
  if(CLc)    model$CLc    <- getCLc(sexes, lengths)
  if(CLs)    model$CLs    <- getCLs(sexes, lengths)
  if(LA)
  {
    latage.file <- paste(dirname(res.file),"l_at_age.dat",sep="/")
    if(!file.exists(latage.file)) stop("File ", latage.file, " not found; use / or \\\\ separators")
    latage.vector <- readLines(latage.file)  # LAnobs for each sex
    latage.vector <- gsub("\"","", gsub("\t","",gsub(" ","",latage.vector)))
    txt.file <- gsub("\\.res", "\\.txt", res.file)
    if(!file.exists(txt.file)) stop("File ", txt.file, " not found; use / or \\\\ separators")
    txt.vector <- readLines(txt.file)  # sigmaLA~x
    txt.vector <- gsub("\"","", gsub("\t","",gsub(" ","",txt.vector)))
    model$LA <- getLA(sexes, ages)
  }
  if(!quiet) cat("\n")

  ## 6  Create attributes
  attr(model,"call") <- match.call()
  attr(model,"scape.version") <- as.character(packageVersion("scape"))
  class(model) <- "scape"

  return(model)
}
