#' Calculate income after tax and benefits
#'
#' Based on current UK taxation rules this function 
#' calculates components that subtract from gross income
#' and provides net income.
#' 
#' Current, in the context of default values, is Tax Year 2014
#' 
#' @param persons                    Provide the information required for calculating income,  
#'                                   values should be provided as annual incomes
#'  
#' @param incomeGrain                Define the time period in which the income return  
#'                                   should be expressed i.e. "Annual", "Month", "Week"
#'                                   
#' @param financialYear              What financial year the calculation should be performed for.
#'                                   Can't go back further than 2014, if you need to go back please
#'                                   submit a pull request on the CSVs in inst/extdata with them filled in.
#'                                   
#' @param modelArgs                  Indicate whether a forward prediction with some 
#'                                   changing values should be performed, and what scenario values should be
#'                                   used
#' 
#' @param thresholdsTable            The values needed for calculating various components
#' 
#' @param taxRateTable               The values needed for calculating Income Tax and NI (Class 1 and 4).
#'                                   Rate tables contain lower bound (LB), upper bound (UB) 
#'                                   and the prevailing tax rates (Rate) at which portions of 
#'                                   income are taxed at. LB >= Income < UB
#' 
#' @return income                    Income components for each person at the relevant grain
#' 
#' @keywords financial tax income
#' @family tax
#' 
#' @export
calcNetIncome <- function(
  persons = data.table(
    personID                  = 1:2, 
    householdID               = 1, 
    employedIncome            = c(15000, 40000), 
    investmentIncome          = c(0, 5000), 
    nonTaxableIncome          = 0, 
    selfEmployedProfits       = 0, 
    taxCode                   = "1000L", 
    numberOfChildren          = 1, 
    salarySacrificePercentage = c(0, 0.05), 
    studentLoan               = 0:1  
  ),
  incomeGrain           = "Month" ,# c("Annual", "Month", "Week")
  financialYear         = taxYear(Sys.Date()),
  modelArgs             = list(model                  = FALSE,
                              inflation               = 1.00, 
                              years                   = 3, 
                              childBenefitChange      = 1.00, 
                              personalAllowanceChange = 500),
  
  thresholdsTable        = fread(system.file("extdata",
                                                 "annualthresholds.csv", package = "optiRum")),
  taxRateTable           = fread(system.file("extdata",
                                                 "annualtaxthresholds.csv", package = "optiRum"))
  ){
  
  # input checks
  persons_cols <- c("householdID", "personID", "employedIncome", "investmentIncome", 
                   "nonTaxableIncome", "selfEmployedProfits", "taxCode", 
                   "numberOfChildren", "salarySacrificePercentage", "studentLoan")
  stopifnot(
    is.data.table(persons),
    nrow(persons) > 0L,
    all(persons_cols %in% names(persons)),
    nrow(persons[,.N,householdID][N>2]) == 0L # only two people allowed per household - significant reduction in complexity
  )
  
  #CRAN check fudge for data.table columns
  householdID <- N <- employedIncome <- investmentIncome <- nonTaxableIncome <- selfEmployedProfits <- 
    salarySacrifice <- salarySacrificePercentage <- personalAllowance <- taxCode <- totalTaxableIncome <- 
    incomeTax <- incomeTaxable <- class1NI <- class1NITaxable <- class4NI <- class4NITaxable <- studentLoan <- 
    generalTaxable <- studentLoanRepayment <- numberOfChildren <- childBenefits <- childBenefitTax <- generalTaxableRank <- 
    householdChildBenefits <- householdChildBenefitTax <- totalIncome <- netIncome <- householdNetIncome <- personID <- 
    TaxYear <- Type <- NULL
  
  # persons table could have more info, reduce it for the purposes of ongoing calcs 
  income <- persons[, .SD, .SDcols = persons_cols]
  
  # get tax year tables
  thresholds       <- thresholdsTable[TaxYear==financialYear]
  taxBrackets      <- taxRateTable[TaxYear==financialYear&Type=="Income Tax"]
  class1NIBrackets <- taxRateTable[TaxYear==financialYear&Type=="NI Class1"]
  class4NIBrackets <- taxRateTable[TaxYear==financialYear&Type=="NI Class4"]
  
  # apply model
  if(modelArgs$model){
    inflation_over_years <- modelArgs$inflation^modelArgs$years
    income[,`:=`(employedIncome      = employedIncome      * inflation_over_years,
                 investmentIncome    = investmentIncome    * inflation_over_years,
                 nonTaxableIncome    = nonTaxableIncome    * inflation_over_years,
                 selfEmployedProfits = selfEmployedProfits * inflation_over_years,
                 taxCode             = "")]
    personalAllowanceValue <- thresholds$personalAllowanceValue + modelArgs$personalAllowanceChange*modelArgs$years
  }
  
  # sum incomes
  income[,`:=`(totalTaxableIncome       = employedIncome + investmentIncome
               ,totalIncome             = employedIncome + investmentIncome + nonTaxableIncome + selfEmployedProfits
               ,class1NationalInsurance = 0
               ,class4NationalInsurance = 0
               ,studentLoanRepayment    = 0
               ,childBenefits           = 0
               ,childBenefitTax         = 0)]
  
  # calc salarySacrifice
  income[,salarySacrifice := employedIncome * salarySacrificePercentage]
  
  # calc personalAllowance from taxCode or default
  
  income[,personalAllowance := taxCode_to_personalAllowance(taxCode)
         ][is.na(personalAllowance),
           personalAllowance := pmin(
             thresholds$personalAllowanceValue,
             pmax(0, thresholds$personalAllowanceValue - (totalTaxableIncome - salarySacrifice - thresholds$personalAllowanceThreshold)/2)
           )]
  
  # calc taxable amount
  income[,`:=`(generalTaxable  = employedIncome + investmentIncome - salarySacrifice,
               incomeTaxable   = pmax(employedIncome + investmentIncome - personalAllowance - salarySacrifice, 0),
               class1NITaxable = employedIncome - salarySacrifice,
               class4NITaxable = selfEmployedProfits)]
  
  # calc tax owed
  income[, incomeTax := TaxOwed(incomeTaxable,   taxBrackets)]
  income[, class1NI  := TaxOwed(class1NITaxable, class1NIBrackets)]
  income[, class4NI  := TaxOwed(class4NITaxable, class4NIBrackets)]
  
  # calc student loan repayment
  income[studentLoan == TRUE & generalTaxable >= thresholds$studentLoanThreshold,
         studentLoanRepayment := (generalTaxable - thresholds$studentLoanThreshold) * thresholds$studentLoanPercentage]
  
  # calc child benefits
  income[numberOfChildren > 0,
         childBenefits := (thresholds$childBenefitChild1 + (numberOfChildren - 1) * thresholds$childBenefitChildS)]
  # 1% less for every 100 over 50,000/year
  # Divide by 100 for every 100 over, then divide by 100 to get as a percent
  income[numberOfChildren > 0 & generalTaxable >= thresholds$childBenefitThreshold,
         childBenefitTax := pmin(
           childBenefits,
           pmax(0, childBenefits * (floor((generalTaxable - thresholds$childBenefitThreshold)/100)/100))
         )]
  
  
  # household child benefits, take max
  income[,`:=`(householdChildBenefits   = max(childBenefits),
               householdChildBenefitTax = max(childBenefitTax),
               childBenefits = 0,
               childBenefitTax = 0),
         by = householdID]
  
  # recalc household benefits - add only to one of the two persons
  income[, generalTaxableRank := frankv(generalTaxable, order=-1L, ties.method="first"), householdID
         ][generalTaxableRank==1L,
           `:=`(childBenefits = householdChildBenefits,
                childBenefitTax = householdChildBenefitTax)
           ][, generalTaxableRank := NULL]
  
  # apply model
  if(modelArgs$model){
    benefit_over_years <- modelArgs$childBenefitChange^modelArgs$years
    income[,`:=`(childBenefits   = childBenefits   * benefit_over_years,
                 childBenefitTax = childBenefitTax * benefit_over_years)]
  }
  
  # take child benefits into account
  income[, totalIncome := totalIncome + childBenefits]
  
  # calculate the final net income
  income[, netIncome := employedIncome + investmentIncome + nonTaxableIncome + 
           (selfEmployedProfits - class4NI) + (childBenefits - childBenefitTax) - 
           (incomeTax + class1NI) - studentLoanRepayment]
  
  # and household income
  income[, householdNetIncome := sum(netIncome), by = householdID]
  
  # adjust for grain
  inputcols<-c("householdID","personID","personalAllowance")
  outputcols<-c("totalIncome","netIncome","householdNetIncome","incomeTax",
                "class1NI","class4NI","childBenefits","childBenefitTax",
                "studentLoanRepayment")
  allcols<-c(inputcols,outputcols)
  
  income[,(outputcols):=lapply(.SD,grainAdjustment,grain=incomeGrain)
         ,.SDcols=outputcols]
  
  # setting the keys for Jan :-)
  income<-setkey(income[,allcols,with=FALSE], householdID, personID)
  return(income)
  
}

taxCode_to_personalAllowance <- function(x) trunc(as.numeric(gsub("[^\\d]+", "", x, perl=TRUE)))*10L

grainAdjustment<-function(x, grain = "Month"){ switch(grain,
                                                      Annual  = x/1
                                                      ,Month  = x/12
                                                      ,Week   = x/52
)  }

TaxOwed <- function(incomeCol, taxTable) {
  # Purpose : Based on someone's annual income, calculate how much they would
  # contribute to a form of taxation
  #
  # Must : Allow for different reference tax tables, handle vectorised income
  #
  # Known : an income column (NI or Tax determined), which rate table to
  # consider, the standard methodology for calculating tax
  stopifnot(
    is.data.table(taxTable),
    all(c("LB","UB","Rate") %in% names(taxTable))
  )
  LB <- UB <- Rate <- NULL
  toPay <- 0
  for (i in seq_len(nrow(taxTable))) {
    toPay <- taxTable[i, toPay + (incomeCol>=LB) * pmin(incomeCol - LB, UB - LB) * Rate]
  }
  toPay
}


