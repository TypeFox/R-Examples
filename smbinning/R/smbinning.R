#' Optimal Binning for Scoring Modeling
#'
#' \strong{Optimal Binning} categorizes a numeric characteristic into bins for ulterior usage in scoring modeling.
#' This process, also known as \emph{supervised discretization}, 
#' utilizes \href{http://cran.r-project.org/package=partykit}{Recursive Partitioning} to categorize 
#' the numeric characteristic.\cr
#' The especific algorithm is Conditional Inference Trees 
#' which initially excludes missing values (\code{NA}) to compute the cutpoints, adding them back later in the 
#' process for the calculation of the \emph{Information Value}.
#' @param df A data frame.
#' @param y Binary response variable (0,1). Integer (\code{int}) is required.
#' Name of \code{y} must not have a dot.
#' @param x Continuous characteristic. At least 10 different values. Value \code{Inf} is not allowed.
#' Name of \code{x} must not have a dot.
#' @param p Percentage of records per bin. Default 5\% (0.05). 
#' This parameter only accepts values greater that 0.00 (0\%) and lower than 0.5 (50\%).
#' @return The command \code{smbinning} generates and object containing the necessary info and utilities for binning.
#' The user should save the output result so it can be used 
#' with \code{smbinning.plot}, \code{smbinning.sql}, and \code{smbinning.gen}.
#' @examples
#' # Package loading and data exploration
#' library(smbinning) # Load package and its data
#' data(chileancredit) # Load smbinning sample dataset (Chilean Credit)
#' str(chileancredit) # Quick description of the data
#' table(chileancredit$FlagGB) # Tabulate target variable
#'  
#' # Training and testing samples (Just some basic formality for Modeling) 
#' chileancredit.train=subset(chileancredit,FlagSample==1)
#' chileancredit.test=subset(chileancredit,FlagSample==0)
#'  
#' # Package application
#' result=smbinning(df=chileancredit.train,y="FlagGB",x="TOB",p=0.05) # Run and save result
#' result$ivtable # Tabulation and Information Value
#' result$iv # Information value
#' result$bands # Bins or bands
#' result$ctree # Decision tree from partykit
smbinning = function(df,y,x,p=0.05){
  # Check data frame and formats
  if (!is.data.frame(df)){ # Check if data.frame
    return("Data is not a data.frame")
  } else if (is.numeric(y) | is.numeric(x)){ # Check if target variable is numeric
    return("Characteristic name not string")
  } else if (grepl("[.]",y) | grepl("[.]",x)){ # Check if there is a dot
    return("Name of a characteristic must not have a dot [.]")
    } else 
    i=which(names(df)==y) # Find Column for dependant
  j=which(names(df)==x) # Find Column for independant
  if (!is.numeric(df[,i])){ 
    return("Target (y) not found or it is not numeric")
  } else if (max(df[,i],na.rm=T)!=1){
    return("Maximum not 1")
  } else if (fn$sqldf("select count(*) from df where cast($x as text)='Inf' or cast($x as text)='-Inf'")>0){
    return("Characteristic (x) with an 'Inf' value (Divided by Zero). Replace by NA")  
  } else if (min(df[,i],na.rm=T)!=0){
    return("Minimum not 0")
  } else if (p<=0 | p>0.5){
    return("p must be greater than 0 and lower than 0.5 (50%)")
  } else if (!is.numeric(df[,j])){
    return("Characteristic (x) not found or it is not a number")
  } else if (length(unique(df[,j]))<10){
    return("Characteristic (x) has less than 10 uniques values")  
  } else { 
    ctree=ctree(formula(paste(y,"~",x)),
                data=df, 
                na.action=na.exclude,
                control=ctree_control(minbucket=ceiling(round(p*nrow(df)))))
    bins=width(ctree)
    if (bins<2){return("No Bins")}
    # Append cutpoints in a table (Automated)
    cutvct=data.frame(matrix(ncol=0,nrow=0)) # Shell
    n=length(ctree) # Number of nodes
    for (i in 1:n) {
      cutvct=rbind(cutvct,ctree[i]$node$split$breaks)
    }
    cutvct=cutvct[order(cutvct[,1]),] # Sort / converts to a ordered vector (asc)
    # Build Information Value Table #############################################
    # Counts per not missing cutpoint
    ivt=data.frame(matrix(ncol=0,nrow=0)) # Shell
    n=length(cutvct) # Number of cutpoits
    for (i in 1:n) {
      cutpoint=cutvct[i]
      ivt=rbind(ivt,
                fn$sqldf(
                  "select '<= $cutpoint' as Cutpoint,
                  NULL as CntRec,
                  NULL as CntGood,
                  NULL as CntBad,
                  sum(case when $x <= $cutpoint and $y in (1,0) then 1 else 0 end) as CntCumRec,
                  sum(case when $x <= $cutpoint and $y=1 then 1 else 0 end) as CntCumGood,
                  sum(case when $x <= $cutpoint and $y=0 then 1 else 0 end) as CntCumBad,
                  NULL as PctRec,
                  NULL as GoodRate,
                  NULL as BadRate,
                  NULL as Odds,
                  NULL as LnOdds,
                  NULL as WoE,
                  NULL as IV
                  from df where $x is not NULL and $y is not NULL")
                )
    }
    cutpoint=max(df[,j],na.rm=T) # Calculte Max without Missing
    mincutpoint=min(df[,j],na.rm=T) # Calculte Min without Missing for later usage
    ivt=rbind(ivt,
              fn$sqldf(
                "select '<= $cutpoint' as Cutpoint,
                NULL as CntRec,
                NULL as CntGood,
                NULL as CntBad,
                sum(case when $x <= $cutpoint and $y in (1,0) then 1 else 0 end) as CntCumRec,
                sum(case when $x <= $cutpoint and $y=1 then 1 else 0 end) as CntCumGood,
                sum(case when $x <= $cutpoint and $y=0 then 1 else 0 end) as CntCumBad,
                NULL as PctRec,
               NULL as GoodRate,
                NULL as BadRate,
                NULL as Odds,
                NULL as LnOdds,
                NULL as WoE,
                NULL as IV
                from df where $x is not NULL and $y is not NULL")
              )
    # Missing Data
    x.na=fn$sqldf("select count(*) from df where $x is null")  
    y.na=fn$sqldf("select count(*) from df where $y is null")
    if(x.na>0){
      ivt=rbind(ivt,
                fn$sqldf(
                  "select 'Missing' as Cutpoint,
                  sum(case when $x is NULL and $y in (1,0) then 1 else 0 end) as CntRec,
                  sum(case when $x is NULL and $y=1 then 1 else 0 end) as CntGood,
                  sum(case when $x is NULL and $y=0 then 1 else 0 end) as CntBad,
                  NULL as CntCumRec,
                  NULL as CntCumGood,
                  NULL as CntCumBad,
                  NULL as PctRec,
                  NULL as GoodRate,
                  NULL as BadRate,
                  NULL as Odds,
                  NULL as LnOdds,
                  NULL as WoE,
                  NULL as IV
                  from df where $y is not NULL")
                )
    } 
    
     else {
      ivt=rbind(ivt,
                c("Missing",0,0,0,NA,NA,NA,NA,NA,NA,NA,NA,NA))}
    
    # Total
    ivt=rbind(ivt,
              fn$sqldf(
                "select 'Total' as Cutpoint,
                count(*) as CntRec,
                sum(case when $y=1 then 1 else 0 end) as CntGood,
                sum(case when $y=0 then 1 else 0 end) as CntBad,
                NULL as CntCumRec,
                NULL as CntCumGood,
                NULL as CntCumBad,
                NULL as PctRec,
                NULL as GoodRate,
                NULL as BadRate,
                NULL as Odds,
                NULL as LnOdds,
                NULL as WoE,
                NULL as IV
                from df where $y is not NULL")
              )
    
    # Covert to table numeric
    ncol=ncol(ivt)
    for (i in 2:ncol){
      ivt[,i]=as.numeric(ivt[,i])
    }
    
    # Complete Table 
    ivt[1,2]=ivt[1,5] # Nbr Records
    ivt[1,3]=ivt[1,6] # Nbr Goods
    ivt[1,4]=ivt[1,7] # Nbr Bads
    
    # From 2nd row
    n=nrow(ivt)-2
    for (i in 2:n){ivt[i,2]=ivt[i,5]-ivt[i-1,5]
                   ivt[i,3]=ivt[i,6]-ivt[i-1,6]
                   ivt[i,4]=ivt[i,7]-ivt[i-1,7]}
    
    ivt[2,2]=ivt[2,5]-ivt[1,5]
    ivt[2,3]=ivt[2,6]-ivt[1,6]
    ivt[2,4]=ivt[2,7]-ivt[1,7]
    
    # Missing row.  Update: Added "if" statement
        ivt[i+1,5]=ivt[i,5]+ivt[i+1,2]
    ivt[i+1,6]=ivt[i,6]+ivt[i+1,3]
    ivt[i+1,7]=ivt[i,7]+ivt[i+1,4]

    # Calculating metrics
    options(scipen=999) # Remove Scientific Notation
    ivt[,8]=round(ivt[,2]/ivt[i+2,2],4) # PctRec
    ivt[,9]=round(ivt[,3]/ivt[,2],4) # GoodRate
    ivt[,10]=round(ivt[,4]/ivt[,2],4) # BadRate
    ivt[,11]=round(ivt[,3]/ivt[,4],4) # Odds
    ivt[,12]=round(log(ivt[,3]/ivt[,4]),4) # LnOdds
    G=ivt[i+2,3]
    B=ivt[i+2,4]
    LnGB=log(G/B) # IV Part 1
    ivt[,13]=round(log(ivt[,3]/ivt[,4])-LnGB,4) # WoE
    ivt[,14]=round(ivt[,13]*(ivt[,3]/G-ivt[,4]/B),4) # Mg IV
    # ivt[i+2,14]=round(sum(ivt[,13]*(ivt[,3]/G-ivt[,4]/B),na.rm=T),4) -- Old Calculation
    # Calculates Information Value even with undefined numbers
    ivt[i+2,14]=0.0000
    for (k in 1:(nrow(ivt)-1))
    {
      if(is.finite(ivt[k,14])) {mgiv=ivt[k,14]} else {mgiv=0.0000}
      ivt[i+2,14]=ivt[i+2,14]+mgiv
    }
    iv=ivt[i+2,14]
    # End Inf. Value Table ###################################################### 
    }
  bands=append(mincutpoint,cutvct)
  bands=append(bands,cutpoint)
  list(ivtable=ivt,iv=iv,ctree=ctree,bands=bands,x=x,col_id=j,cuts=cutvct)
  }

# Begin Custom Cutpoints 20150307 #############################################
#' Customized Binning
#'
#' It gives the user the ability to create customized cutpoints. In Scoring Modeling, the analysis
#' of a characteristic usually begins with intervals with the same length to understand its distribution,
#' and then intervals with the same proportion of cases to explore bins with a reasonable sample size.
#' @param df A data frame.
#' @param y Binary response variable (0,1). Integer (\code{int}) is required.
#' Name of \code{y} must not have a dot.
#' @param x Continuous characteristic. At least 10 different values. Value \code{Inf} is not allowed. 
#' Name of \code{x} must not have a dot.
#' @param cuts Vector with the cutpoints selected by the user. It does not have a default so user must define it. 
#' @return The command \code{smbinning.custom} generates and object containing the necessary info and utilities for binning.
#' The user should save the output result so it can be used 
#' with \code{smbinning.plot}, \code{smbinning.sql}, and \code{smbinning.gen}.
#' @examples
#' # Package loading and data exploration
#' library(smbinning) # Load package and its data
#' data(chileancredit) # Load smbinning sample dataset (Chilean Credit)
#' str(chileancredit) # Quick description of the data
#' table(chileancredit$FlagGB) # Tabulate target variable
#'  
#' # Training and testing samples (Just some basic formality for Modeling) 
#' chileancredit.train=subset(chileancredit,FlagSample==1)
#' chileancredit.test=subset(chileancredit,FlagSample==0)
#'
#' # Remove exclusions from chileancredit dataset
#' TOB.train=
#'   subset(chileancredit,(FlagSample==1 & (FlagGB==1 | FlagGB==0)), select=TOB)
#' TOB.test=
#'   subset(chileancredit,(FlagSample==0 & (FlagGB==1 | FlagGB==0)), select=TOB)
#'   
#' # Custom cutpoints using percentiles (20% each)
#' TOB.Pct20=quantile(TOB.train, probs=seq(0,1,0.2), na.rm=TRUE)
#' TOB.Pct20.Breaks=as.vector(quantile(TOB.train, probs=seq(0,1,0.2), na.rm=TRUE))
#' Cuts.TOB.Pct20=TOB.Pct20.Breaks[2:(length(TOB.Pct20.Breaks)-1)]
#' 
#' # Package application and results
#' result=
#'   smbinning.custom(df=chileancredit.train,
#'                    y="FlagGB",x="TOB",cuts=Cuts.TOB.Pct20) # Run and save
#' result$ivtable # Tabulation and Information Value
smbinning.custom = function(df,y,x,cuts){
  # Check data frame and formats
  if (!is.data.frame(df)){ # Check if data.frame
    return("Data is not a data.frame")
  } else if (is.numeric(y) | is.numeric(x)){ # Check if target vable is numeric
    return("Characteristic name not string")
  } else if (grepl("[.]",y) | grepl("[.]",x)){ # Check if there is a dot
    return("Name of a characteristic must not have a dot [.]")
  } else 
    i=which(names(df)==y) # Find Column for dependant
  j=which(names(df)==x) # Find Column for independant
  if (!is.numeric(df[,i])){ 
    return("Target (y) not found or it is not numeric")
  } else if (max(df[,i],na.rm=T)!=1){
    return("Maximum not 1")
  } else if (fn$sqldf("select count(*) from df where cast($x as text)='Inf' or cast($x as text)='-Inf'")>0){
    return("Characteristic (x) with an 'Inf' value (Divided by Zero). Replace by NA")  
  } else if (min(df[,i],na.rm=T)!=0){
    return("Minimum not 0")
  } else if (!is.numeric(df[,j])){
    return("Characteristic (x) not found or it is not a number")
  } else if (length(unique(df[,j]))<10){
    return("Characteristic (x) has less than 10 uniques values")  
  } else { 
    # Append cutpoints in a table (Automated)
    cutvct=data.frame(matrix(ncol=0,nrow=0)) # Shell
    n=length(cuts) # Number of cutpoints
    if (n<1){return("No Bins")} # At least 1 cutpoint
    for (i in 1:n) {
      cutvct=rbind(cutvct,cuts[i])
    }
    cutvct=cutvct[order(cutvct[,1]),] # Sort / converts to a ordered vector (asc)
    # Build Information Value Table #############################################
    # Counts per not missing cutpoint
    ivt=data.frame(matrix(ncol=0,nrow=0)) # Shell
    n=length(cutvct) # Number of cutpoits
    for (i in 1:n) {
      cutpoint=cutvct[i]
      ivt=rbind(ivt,
                fn$sqldf(
                  "select '<= $cutpoint' as Cutpoint,
                  NULL as CntRec,
                  NULL as CntGood,
                  NULL as CntBad,
                  sum(case when $x <= $cutpoint and $y in (1,0) then 1 else 0 end) as CntCumRec,
                  sum(case when $x <= $cutpoint and $y=1 then 1 else 0 end) as CntCumGood,
                  sum(case when $x <= $cutpoint and $y=0 then 1 else 0 end) as CntCumBad,
                  NULL as PctRec,
                  NULL as GoodRate,
                  NULL as BadRate,
                  NULL as Odds,
                  NULL as LnOdds,
                  NULL as WoE,
                  NULL as IV
                  from df where $x is not NULL and $y is not NULL")
                )
    }
    cutpoint=max(df[,j],na.rm=T) # Calculte Max without Missing
    mincutpoint=min(df[,j],na.rm=T) # Calculte Min without Missing for later usage
    ivt=rbind(ivt,
              fn$sqldf(
                "select '= $cutpoint' as Cutpoint,
                NULL as CntRec,
                NULL as CntGood,
                NULL as CntBad,
                sum(case when $x <= $cutpoint and $y in (1,0) then 1 else 0 end) as CntCumRec,
                sum(case when $x <= $cutpoint and $y=1 then 1 else 0 end) as CntCumGood,
                sum(case when $x <= $cutpoint and $y=0 then 1 else 0 end) as CntCumBad,
                NULL as PctRec,
                NULL as GoodRate,
                NULL as BadRate,
                NULL as Odds,
                NULL as LnOdds,
                NULL as WoE,
                NULL as IV
                from df where $x is not NULL and $y is not NULL")
              )
    # Missing Data
    x.na=fn$sqldf("select count(*) from df where $x is null")  
    y.na=fn$sqldf("select count(*) from df where $y is null")
    if(x.na>0){
      ivt=rbind(ivt,
                fn$sqldf(
                  "select 'Missing' as Cutpoint,
                  sum(case when $x is NULL and $y in (1,0) then 1 else 0 end) as CntRec,
                  sum(case when $x is NULL and $y=1 then 1 else 0 end) as CntGood,
                  sum(case when $x is NULL and $y=0 then 1 else 0 end) as CntBad,
                  NULL as CntCumRec,
                  NULL as CntCumGood,
                  NULL as CntCumBad,
                  NULL as PctRec,
                  NULL as GoodRate,
                  NULL as BadRate,
                  NULL as Odds,
                  NULL as LnOdds,
                  NULL as WoE,
                  NULL as IV
                  from df where $y is not NULL")
                )
    } else {
      ivt=rbind(ivt,
                c("Missing",0,0,0,NA,NA,NA,NA,NA,NA,NA,NA,NA))
    }
    # Total
    ivt=rbind(ivt,
              fn$sqldf(
                "select 'Total' as Cutpoint,
                count(*) as CntRec,
                sum(case when $y=1 then 1 else 0 end) as CntGood,
                sum(case when $y=0 then 1 else 0 end) as CntBad,
                NULL as CntCumRec,
                NULL as CntCumGood,
                NULL as CntCumBad,
                NULL as PctRec,
                NULL as GoodRate,
                NULL as BadRate,
                NULL as Odds,
                NULL as LnOdds,
                NULL as WoE,
                NULL as IV
                from df where $y is not NULL")
              )
    
    # Covert to table numeric
    ncol=ncol(ivt)
    for (i in 2:ncol){
      ivt[,i]=as.numeric(ivt[,i])
    }
    
    # Complete Table 
    ivt[1,2]=ivt[1,5] # Nbr Records
    ivt[1,3]=ivt[1,6] # Nbr Goods
    ivt[1,4]=ivt[1,7] # Nbr Bads
    
    # From 2nd row
    n=nrow(ivt)-2
    for (i in 2:n){ivt[i,2]=ivt[i,5]-ivt[i-1,5]
                   ivt[i,3]=ivt[i,6]-ivt[i-1,6]
                   ivt[i,4]=ivt[i,7]-ivt[i-1,7]}
    
    ivt[2,2]=ivt[2,5]-ivt[1,5]
    ivt[2,3]=ivt[2,6]-ivt[1,6]
    ivt[2,4]=ivt[2,7]-ivt[1,7]
    
    # Missing row
    ivt[i+1,5]=ivt[i,5]+ivt[i+1,2]
    ivt[i+1,6]=ivt[i,6]+ivt[i+1,3]
    ivt[i+1,7]=ivt[i,7]+ivt[i+1,4]
    
    # Calculating metrics
    options(scipen=999) # Remove Scientific Notation
    ivt[,8]=round(ivt[,2]/ivt[i+2,2],4) # PctRec
    ivt[,9]=round(ivt[,3]/ivt[,2],4) # GoodRate
    ivt[,10]=round(ivt[,4]/ivt[,2],4) # BadRate
    ivt[,11]=round(ivt[,3]/ivt[,4],4) # Odds
    ivt[,12]=round(log(ivt[,3]/ivt[,4]),4) # LnOdds
    G=ivt[i+2,3]
    B=ivt[i+2,4]
    LnGB=log(G/B) # IV Part 1
    ivt[,13]=round(log(ivt[,3]/ivt[,4])-LnGB,4) # WoE
    ivt[,14]=round(ivt[,13]*(ivt[,3]/G-ivt[,4]/B),4) # Mg IV
    # ivt[i+2,14]=round(sum(ivt[,13]*(ivt[,3]/G-ivt[,4]/B),na.rm=T),4) -- Old Calculation
    # Calculates Information Value even with undefined numbers
    ivt[i+2,14]=0.0000
    for (k in 1:(nrow(ivt)-1))
    {
      if(is.finite(ivt[k,14])) {mgiv=ivt[k,14]} else {mgiv=0.0000}
      ivt[i+2,14]=ivt[i+2,14]+mgiv
    }
    iv=ivt[i+2,14]
    # End Inf. Value Table ###################################################### 
    }
  bands=append(mincutpoint,cutvct)
  bands=append(bands,cutpoint)
  list(ivtable=ivt,iv=iv,bands=bands,x=x,col_id=j,cuts=cutvct)
  }
# End Custom Cutpoints 20150307 ###############################################


# Begin Binning Factors 20150407 #############################################
#' Binning on Factor Variables
#'
#' It generates the output table for the uniques values of a given factor variable.
#' @param df A data frame.
#' @param y Binary response variable (0,1). Integer (\code{int}) is required.
#' Name of \code{y} must not have a dot.
#' @param x A factor variable with at least 2 different values. Value \code{Inf} is not allowed. 
#' Name of \code{x} must not have a dot.
#' @return The command \code{smbinning.factor} generates and object containing the necessary info and utilities for binning.
#' The user should save the output result so it can be used 
#' with \code{smbinning.plot}, \code{smbinning.sql}, and \code{smbinning.gen}.
#' @examples
#' # Package loading and data exploration
#' library(smbinning) # Load package and its data
#' data(chileancredit) # Load smbinning sample dataset (Chilean Credit)
#' str(chileancredit) # Quick description of the data
#' table(chileancredit$FlagGB) # Tabulate target variable
#' 
#' # Data transformation. Data type must be factor.
#' chileancredit$IncomeLevel= factor(chileancredit$IncomeLevel, 
#'                                   levels=c(0,1,2,3,4,5),
#'                                   labels=c("00","01","02","03","04","05"))
#'  
#' # Training and testing samples (Just some basic formality for Modeling) 
#' chileancredit.train=subset(chileancredit,FlagSample==1)
#' chileancredit.test=subset(chileancredit,FlagSample==0)
#' 
#' # Package application and results
#' result.train=smbinning.factor(df=chileancredit.train,
#'                                y="FlagGB",x="IncomeLevel")
#' result.train$ivtable
#' result.test=smbinning.factor(df=chileancredit.test,
#'                                y="FlagGB",x="IncomeLevel")
#' result.test$ivtable
#' 
#' # Plots
#' par(mfrow=c(2,2))
#' smbinning.plot(result.train,option="dist",sub="Income Level (Tranining Sample)")
#' smbinning.plot(result.train,option="badrate",sub="Income Level (Tranining Sample)")
#' smbinning.plot(result.test,option="dist",sub="Income Level (Test Sample)")
#' smbinning.plot(result.test,option="badrate",sub="Income Level (Test Sample)")

smbinning.factor = function(df,y,x){
  # Check data frame and formats
  if (!is.data.frame(df)){ # Check if data.frame
    return("Data is not a data.frame")
  } else if (is.numeric(y) | is.numeric(x)){ # Check if target vable is numeric
    return("Characteristic name not string")
  } else if (grepl("[.]",y) | grepl("[.]",x)){ # Check if there is a dot
    return("Name of a characteristic must not have a dot [.]")
  } else 
    i=which(names(df)==y) # Find Column for dependant
  j=which(names(df)==x) # Find Column for independant
  if (!is.numeric(df[,i])){ 
    return("Target (y) not found or it is not numeric")
  } else if (max(df[,i],na.rm=T)!=1){
    return("Maximum not 1")
  } else if (fn$sqldf("select count(*) from df where cast($x as text)='Inf' or cast($x as text)='-Inf'")>0){
    return("Characteristic (x) with an 'Inf' value (Divided by Zero). Replace by NA")  
  } else if (min(df[,i],na.rm=T)!=0){
    return("Minimum not 0")
  } else if (!is.factor(df[,j])){
    return("Characteristic (x) not found or it is not a factor")
  } else if (length(unique(df[,j]))<=1){
    return("Characteristic (x) requires at leats 2 uniques values")  
  } else { 
    # Append cutpoints in a table (Automated)
    # cutvct=data.frame(matrix(ncol=0,nrow=0)) # Shell
    cutvct=c()
    cuts=fn$sqldf("select distinct $x from df where $x is not NULL")
    cuts=as.vector(as.matrix(cuts))
    n=length(cuts) # Number of cutpoints
    if (n<1){return("No Bins")} # At least 1 cutpoint
    for (i in 1:n) {
      cutvct=rbind(cutvct,cuts[i])
    }
    cutvct=cutvct[order(cutvct[,1]),] # Sort / converts to a ordered vector (asc)
    # Build Information Value Table #############################################
    # Counts per not missing cutpoint
    ivt=data.frame(matrix(ncol=0,nrow=0)) # Shell
    n=length(cutvct) # Number of cutpoits
    for (i in 1:n) {
      cutpoint=cutvct[i]
      ivt=rbind(ivt,
                fn$sqldf(
                  "select '= $cutpoint' as Cutpoint,
                  sum(case when $x = '$cutpoint' and $y in (1,0) then 1 else 0 end) as CntRec,
                  sum(case when $x = '$cutpoint' and $y=1 then 1 else 0 end) as CntGood,
                  sum(case when $x = '$cutpoint' and $y=0 then 1 else 0 end) as CntBad,
                  NULL as CntCumRec,
                  NULL as CntCumGood,
                  NULL as CntCumBad,
                  NULL as PctRec,
                  NULL as GoodRate,
                  NULL as BadRate,
                  NULL as Odds,
                  NULL as LnOdds,
                  NULL as WoE,
                  NULL as IV
                  from df where $x is not NULL and $y is not NULL")
                )
    }
    # Missing Data
    x.na=fn$sqldf("select count(*) from df where $x is null")  
    y.na=fn$sqldf("select count(*) from df where $y is null")
    if(x.na>0){
      ivt=rbind(ivt,
                fn$sqldf(
                  "select 'Missing' as Cutpoint,
                  sum(case when $x is NULL and $y in (1,0) then 1 else 0 end) as CntRec,
                  sum(case when $x is NULL and $y=1 then 1 else 0 end) as CntGood,
                  sum(case when $x is NULL and $y=0 then 1 else 0 end) as CntBad,
                  NULL as CntCumRec,
                  NULL as CntCumGood,
                  NULL as CntCumBad,
                  NULL as PctRec,
                  NULL as GoodRate,
                  NULL as BadRate,
                  NULL as Odds,
                  NULL as LnOdds,
                  NULL as WoE,
                  NULL as IV
                  from df where $y is not NULL")
                )
    } else {
      ivt=rbind(ivt,
                c("Missing",0,0,0,NA,NA,NA,NA,NA,NA,NA,NA,NA))
    }
    # Total
    ivt=rbind(ivt,
              fn$sqldf(
                "select 'Total' as Cutpoint,
                count(*) as CntRec,
                sum(case when $y=1 then 1 else 0 end) as CntGood,
                sum(case when $y=0 then 1 else 0 end) as CntBad,
                NULL as CntCumRec,
                NULL as CntCumGood,
                NULL as CntCumBad,
                NULL as PctRec,
                NULL as GoodRate,
                NULL as BadRate,
                NULL as Odds,
                NULL as LnOdds,
                NULL as WoE,
                NULL as IV
                from df where $y is not NULL")
              )
    
    # Covert table to numeric
    ncol=ncol(ivt)
    for (i in 2:ncol){
      ivt[,i]=as.numeric(ivt[,i])
    }
    
    # Complete Table: 1st row 
    ivt[1,5]=ivt[1,2] # Nbr Cum. Records
    ivt[1,6]=ivt[1,3] # Nbr Cum. Goods
    ivt[1,7]=ivt[1,4] # Nbr Cum. Bads
    
    # From 2nd row
    n=nrow(ivt)-2
    for (i in 2:n){ivt[i,5]=ivt[i,2]+ivt[i-1,5]
                   ivt[i,6]=ivt[i,3]+ivt[i-1,6]
                   ivt[i,7]=ivt[i,4]+ivt[i-1,7]}
    
    ivt[2,5]=ivt[2,2]+ivt[1,5]
    ivt[2,6]=ivt[2,3]+ivt[1,6]
    ivt[2,7]=ivt[2,4]+ivt[1,7]
    
    # Missing row
    ivt[i+1,5]=ivt[i,5]+ivt[i+1,2]
    ivt[i+1,6]=ivt[i,6]+ivt[i+1,3]
    ivt[i+1,7]=ivt[i,7]+ivt[i+1,4]
    
    # Calculating metrics
    options(scipen=999) # Remove Scientific Notation
    ivt[,8]=round(ivt[,2]/ivt[i+2,2],4) # PctRec
    ivt[,9]=round(ivt[,3]/ivt[,2],4) # GoodRate
    ivt[,10]=round(ivt[,4]/ivt[,2],4) # BadRate
    ivt[,11]=round(ivt[,3]/ivt[,4],4) # Odds
    ivt[,12]=round(log(ivt[,3]/ivt[,4]),4) # LnOdds
    G=ivt[i+2,3]
    B=ivt[i+2,4]
    LnGB=log(G/B) # IV Part 1
    ivt[,13]=round(log(ivt[,3]/ivt[,4])-LnGB,4) # WoE
    ivt[,14]=round(ivt[,13]*(ivt[,3]/G-ivt[,4]/B),4) # Mg IV
    # ivt[i+2,14]=round(sum(ivt[,13]*(ivt[,3]/G-ivt[,4]/B),na.rm=T),4) -- Old Calculation
    # Calculates Information Value even with undefined numbers
    ivt[i+2,14]=0.0000
    for (k in 1:(nrow(ivt)-1))
    {
      if(is.finite(ivt[k,14])) {mgiv=ivt[k,14]} else {mgiv=0.0000}
      ivt[i+2,14]=ivt[i+2,14]+mgiv
    }
    iv=ivt[i+2,14]
    # End Inf. Value Table ###################################################### 
    }
  list(ivtable=ivt,iv=iv,x=x,col_id=j,cuts=cutvct)
  }

# End Binning Factors 20150407 #################################################


# Begin Plotting ##############################################################
#' Plots after binning
#'
#' It generates plots for distribution, bad rate, and weight of evidence after running \code{smbinning} 
#' and saving its output.
#' @param ivout An object generated by binning.
#' @param option Distribution ("dist"), Good Rate ("goodrate"), Bad Rate ("badrate"), and Weight of Evidence ("WoE").
#' @param sub Subtitle for the chart (optional).
#' @examples
#' # Package loading and data exploration
#' library(smbinning) # Load package and its data
#' data(chileancredit) # Load smbinning sample dataset (Chilean Credit)
#' chileancredit.train=subset(chileancredit,FlagSample==1)
#' chileancredit.test=subset(chileancredit,FlagSample==0)
#' result=smbinning(df=chileancredit.train,y="FlagGB",x="TOB",p=0.05) # Run and save result
#' 
#' # Plots
#' par(mfrow=c(2,2))
#' boxplot(chileancredit.train$TOB~chileancredit.train$FlagGB,
#'         horizontal=TRUE, frame=FALSE, col="lightgray",main="Distribution")
#' mtext("Time on Books (Months)",3)
#' smbinning.plot(result,option="dist",sub="Time on Books (Months)")
#' smbinning.plot(result,option="badrate",sub="Time on Books (Months)")
#' smbinning.plot(result,option="WoE",sub="Time on Books (Months)")
smbinning.plot=function(ivout,option="dist",sub=""){
  r=ifelse(ivout$ivtable[nrow(ivout$ivtable)-1,2]==0,2,1)
  if (option=="dist"){
    # Distribution
    x_upper=nrow(ivout$ivtable)-r
    y_upper=max(ivout$ivtable[1:x_upper,8])*1.25
    ch_dist=barplot(ivout$ivtable[1:x_upper,8],
                    names.arg=ivout$ivtable[1:x_upper,1],
                    axes=F, 
                    main="Percentage of Cases", 
                    ylim=c(0,y_upper),
                    col=gray.colors(length(unique(ivout$ivtable[1:x_upper,1]))))
    text(x=ch_dist,y=ivout$ivtable[1:x_upper,8], label=round(ivout$ivtable[1:x_upper,8]*100,1), pos=3,cex=1)
    abline(h=0)
    mtext(sub,3)
  } else if (option=="goodrate"){
    # Good Rate
    x_upper=nrow(ivout$ivtable)-r
    y_upper=max(ivout$ivtable[1:x_upper,9],na.rm=T)*1.25
    ch_goodrate=barplot(ivout$ivtable[1:x_upper,9],
                        names.arg=ivout$ivtable[1:x_upper,1],
                        axes=F, 
                        main="Good Rate (%)",
                        ylim=c(0,y_upper),
                        col=gray.colors(length(unique(ivout$ivtable[1:x_upper,1]))))
    text(x=ch_goodrate,y=ivout$ivtable[1:x_upper,9], label=round(ivout$ivtable[1:x_upper,9]*100,1), pos=3,cex=1)
    abline(h=0)
    mtext(sub,3)
  } else if (option=="badrate"){
    # Bad Rate
    x_upper=nrow(ivout$ivtable)-r
    y_upper=max(ivout$ivtable[1:x_upper,10],na.rm=T)*1.25
    ch_badrate=barplot(ivout$ivtable[1:x_upper,10],
                       names.arg=ivout$ivtable[1:x_upper,1],
                       axes=F, 
                       main="Bad Rate (%)",
                       ylim=c(0,y_upper),
                       col=gray.colors(length(unique(ivout$ivtable[1:x_upper,1]))))
    text(x=ch_badrate,y=ivout$ivtable[1:x_upper,10], label=round(ivout$ivtable[1:x_upper,10]*100,1), pos=3,cex=1)
    abline(h=0)
    mtext(sub,3)
  } else if (option=="WoE") {
    # WoE
    x_upper=nrow(ivout$ivtable)-r
    y_upper=max(ivout$ivtable[1:x_upper,13],na.rm=T)*1.25
    y_lower=min(ivout$ivtable[1:x_upper,13],na.rm=T)*1.25
    ch_woe=barplot(ivout$ivtable[1:x_upper,13],
                   names.arg=ivout$ivtable[1:x_upper,1],
                   axes=F, 
                   main="Weight of Evidence", 
                   ylim=c(y_lower,y_upper),
                   col=gray.colors(length(unique(ivout$ivtable[1:x_upper,1]))))
    text(x=ch_woe,y=ivout$ivtable[1:x_upper,13], label=round(ivout$ivtable[1:x_upper,13],2), pos=3,cex=1)
    abline(h=0)
    mtext(sub,3)
  } else {
    return("Options are dist, goodrate, badrate, or WoE")
  }
}
# End Plotting ################################################################

# Begin Gen Characteristic #####################################################
#' Utility to generate a new characteristic
#'
#' It generates a data frame with a new predictive characteristic after the binning process.
#' @param df Dataset to be updated with the new characteristic.
#' @param ivout An object generated after \code{smbinning}.
#' @param chrname Name of the new characteristic.
#' @return A data frame with the binned version of the characteristic analyzed with \code{smbinning}.
#' @examples
#' # Package loading and data exploration
#' library(smbinning) # Load package and its data
#' data(chileancredit) # Load smbinning sample dataset (Chilean Credit)
#' chileancredit.train=subset(chileancredit,FlagSample==1)
#' chileancredit.test=subset(chileancredit,FlagSample==0)
#' result=smbinning(df=chileancredit.train,y="FlagGB",x="TOB",p=0.05) # Run and save result
#' 
#' # Generate new binned characteristic into a existing data frame
#' chileancredit.train=
#' smbinning.gen(chileancredit.train,result,"gTOB") # Update training sample
#' chileancredit=
#'   smbinning.gen(chileancredit,result,"gTOB") # Update population
#' sqldf("select gTOB,count(*) as Recs 
#'       from chileancredit group by gTOB") # Check new field counts 
smbinning.gen=function(df,ivout,chrname="NewChar"){
  df=cbind(df,tmpname=NA)
  ncol=ncol(df)
  col_id=ivout$col_id
  df[,ncol][df[,col_id]>=ivout$bands[1] & df[,col_id]<=ivout$bands[2]]=paste(sprintf("%02d",0),ivout$x,"<=",ivout$bands[2])
  dim=length(ivout$bands)-1
  for (i in 2:dim){
    df[,ncol][df[,col_id]>ivout$bands[i] & df[,col_id]<=ivout$bands[i+1]]=paste(sprintf("%02d",i-1),ivout$x,"<=",ivout$bands[i+1])
  }
  df[,ncol][is.na(df[,col_id])]=paste("99",ivout$x,"Is Null")
  names(df)[names(df)=="tmpname"]=chrname
  return(df)
}
# End Gen Characteristic #######################################################

# Begin: SQL Code #############################################################
#' SQL Code
#'
#' It outputs a SQL code to facilitate the generation of new binned characetristic 
#' in a SQL environment.
#' @param ivout An object generated by \code{smbinning}.
#' @return A text with the SQL code for binning.
#' @examples
#' # Package loading and data exploration
#' library(smbinning) # Load package and its data
#' data(chileancredit) # Load smbinning sample dataset (Chilean Credit)
#' chileancredit.train=subset(chileancredit,FlagSample==1)
#' chileancredit.test=subset(chileancredit,FlagSample==0)
#' result=smbinning(df=chileancredit.train,y="FlagGB",x="TOB",p=0.05) # Run and save result
#' 
#' # Generate SQL code
#' smbinning.sql(result)
smbinning.sql=function(ivout){
  lines=nrow(ivout$ivtable)-2
  sqlcodetable=as.list(matrix(ncol=0,nrow=0))
  sqlcodetable=rbind("case")
  for (k in 1:lines){
    sqlcodetable=rbind(sqlcodetable,paste("when",ivout$x,ivout$ivtable[k,1],"then","'",sprintf("%02d",k),":",ivout$x,ivout$ivtable[k,1],"'"))
  }
  sqlcodetable=rbind(sqlcodetable,paste("when",ivout$x,"Is Null then","'",ivout$x,"Is Null'"))
  sqlcodetable=rbind(sqlcodetable,paste("else '99: Error' end"))
  sqlcodetable
  # Some Make Up
  sqlcodetable=gsub(" '","'",sqlcodetable)
  sqlcodetable=gsub("' ","'",sqlcodetable)
  sqlcodetable=gsub("then","then ",sqlcodetable)
  sqlcodetable=gsub("else","else ",sqlcodetable)
  sqlcodetable=gsub("'end","' end ",sqlcodetable)
  sqlcodetable=gsub(" :",":",sqlcodetable)
  return(gsub(",","",toString(sqlcodetable)))
}
# End: SQL Code ###############################################################


# Begin: Chilean Credit Data ##################################################
#' Chilean Credit Data 
#'
#' A simulated dataset based on six months of information collected by a Chilean Bank
#' whose objective was to develop a credit scoring model to determine the probability
#' of default within the next 12 months. The target 
#' variable is FlagGB, which represents the binary status of default (0) and not default(1).
#'
#' \itemize{
#'   \item CustomerId. Customer Identifier.
#'   \item TOB. Time on books in months since first account was open.
#'   \item IncomeLevel. Income level from 00 (Low) to 05 (High).
#'   \item RevAccts. Number of open revolving accounts.
#'   \item InsAccs. Number of open installment accounts.
#'   \item RevBal. Outstanding balance in all open revolving accounts
#'   \item InsBal. Outstanding balance in all open installment accounts
#'   \item RevLim. Limit of all open revolving accounts
#'   \item SavingsAmtL3M. Amount saved in the last 3 months.
#'   \item Bal. Outstanding balance
#'   \item MaxDqBin. Max. delinquency bin. 0:No Dq., 1:1-29 ... 6:150-179.
#'   \item BureauBalDq1st. Outstanding balance 30-89 in Credit Bureau.
#'   \item BureauBalDq2nd.Outstanding balance 90-179 in Credit Bureau.
#'   \item BureauBalDq3rd.Outstanding balance 180+ in Credit Bureau.
#'   \item CntOtherLenders. Numer of other lenders (Credit Bureau).
#'   \item MtgBal. Mortgage outstanding balance at the Credit Bureau.
#'   \item NonBankTradesDq. Number of non-bank delinquent trades.
#'   \item Performance. Scoring model performance definition.
#'   \item CatGBI. Exclusion (NULL), Bad (00), Indet. (01), Good (02).
#'   \item FlagGB. 1: Good, 0: Bad.
#'   \item Random. Uniformly distributed random value for sampling purposes.
#'   \item FlagSample. Training and test sample indicator (1:75\%,0:25\%).
#'   }
#'
#'
#' @format Data frame with 29,519 rows and 57 columns.
#' @name chileancredit
NULL
# End: Chilean Credit Data ####################################################
