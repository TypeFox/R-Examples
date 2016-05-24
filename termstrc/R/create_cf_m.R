
## Cashflows matrix creation function 

create_cashflows_matrix <- function(group,include_price=FALSE) {
  
 
  n_of_cf <- summary(factor(group$CASHFLOWS$ISIN,levels=group$ISIN),maxsum=1000)
  n_of_bonds <- length(n_of_cf)
  max_cf <- max(n_of_cf)
  pos_cf <- c(0,cumsum(n_of_cf))
  
  # the number of rows of the matrix is the number of 
  # cashflows of the bond with the maximum of cashflows
  # all missing elements of the matrix are filled up with zeros 
       
  CASHFLOWMATRIX <-
    mapply(function(i) c(group$CASHFLOWS$CF[(pos_cf[i]+1):pos_cf[i+1]],
                         rep(0,max_cf-n_of_cf[i])),
           1:n_of_bonds)
  
  if (include_price == TRUE) {CASHFLOWMATRIX <- rbind(-(group[["PRICE"]] +
        group[["ACCRUED"]] ),CASHFLOWMATRIX)}              
  colnames(CASHFLOWMATRIX) <- group$ISIN
  CASHFLOWMATRIX
}

## Maturities matrix creation function 

create_maturities_matrix <- function(group,include_price=FALSE) {

  n_of_cf <- summary(factor(group$CASHFLOWS$ISIN,levels=group$ISIN),maxsum=1000)
  n_of_bonds <- length(n_of_cf)
  max_cf <- max(n_of_cf)
  pos_cf <- c(0,cumsum(n_of_cf))

  year_diff <- as.numeric(difftime(as.Date(group$CASHFLOW$DATE), 
        as.Date(group$TODAY), units = "days"))/365
  
  # RQuantLib version
  # DayCounter: 2 ActualActual
  # year_diff <- mapply(function(i) yearFraction(as.Date(group$TODAY),
  #                          as.Date(group$CASHFLOW$DATE)[[i]], 2), 1:sum(n_of_cf))
  
  # the number of rows of the matrix is the number of 
  # maturity dates of the bond with the longest maturity
  # all missing elements of the matrix are filled up with zeros 
  
  MATURITYMATRIX <-
     mapply(function(i) c(year_diff[(pos_cf[i]+1):pos_cf[i+1]],
                     rep(0,max_cf-n_of_cf[i])),
            1:n_of_bonds)
  
  if (include_price == TRUE) {MATURITYMATRIX <- rbind(rep(0,n_of_bonds),
                           MATURITYMATRIX)}  
  colnames(MATURITYMATRIX) <- group$ISIN
  MATURITYMATRIX
}
