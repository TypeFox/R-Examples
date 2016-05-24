#' Clean up the input data set
#' 
#' This function cleans up the input data set in the backtestGraphics function.
#' This function checks if all the columns required by the backtestGraphics
#' function exist in the input data set. If any columns are missing, this
#' function will try to calculate or fill in the missing columns. If the input
#' data set misses too much data/columns, this function will return an error 
#' indicating the missing columns.
#' 
#' @param x is a data frame that contains the necessary information.
#' @param name.var is the column name of the instrument name column in the input 
#' data set, character type. The default value of this variable is "name". The user has to
#' specify name.var if she passes in a data frame with a different column name for the 
#' instrument name. The input data set must contain either a name column or an ID column.
#' @param id.var is the column name of the instrument ID column, character type. 
#' The default value of this variable is "id". If the input data set labels the
#' column of instrument ID's with some other column name, the user has to pass
#' in the column name here. The input data set has to contain either a name
#' column or an ID column.
#' @param date.var is the column name of the date column, character type. The
#' default value of this variable is "date". This column has to exist, and the
#' column name has to be correct in order for the function to process the data
#' set properly.
#' 
#' @param nmv.var is the column name of the "net market value" column, character 
#' type. The default value of this variable is "start.nmv". The input data set
#' has to contain either a "net market value" column or a "number of contract"
#' column.
#' @param gmv.var The column name of the "gross market value" column if exists, 
#' character type. The default value of this variable is "gmv". Such column will 
#' be automatically calculated from the "net market value" column if it does not 
#' exist.
#' @param pnl.var The column name of the profit-and-loss column, character type.
#' The default value of this variable is "pnl.adj". The data set has to contain
#' such column so that the function can function properly. If such column is
#' missing, the function will return an error indicating the problem.
#' @param contract.var The column name of the contract number column, character 
#' type. The default value of this variable is "num.contract.start". If such 
#' column is missing, the function will instead use net market value as the
#' contract number for each day.
#' 
#' @param sector.var is the column name of the sector column, character type. The
#' default value of this variable is "sector". The sector column helps the user
#' to group instruments into big groups according to industries. The function
#' will still perform properly if the column name for sector is wrong.
#' @param strategy.var The column name of the strategy column, if any. Character 
#' type. The default value of this variable is "strategy". This column can be
#' missing from the input data set.
#' @param substrategy.var The column name of the substrategy column, if any. 
#' Character type. The default value of this variable is "substrategy". This
#' column can be missing from the input data set.
#' @param portfolio.var The column name of the portfolio number column, if any. 
#' Character type. The default value of this variable is "portfolio". This
#' column can be missing from the input data set.
#' 
#' @return x A data frame that only contains the needed data. The column names
#' are cleaned up here.

cleanup_column <- function(x,
                           name.var,   
                           id.var,
                           date.var,
                           nmv.var,
                           gmv.var,
                           pnl.var,
                           contract.var,
                           sector.var,
                           strategy.var,
                           substrategy.var,
                           portfolio.var){
  
  ## Set up variables to fool R CMD check
  
  nmv <- contract <- nmv.temp <- name <- id <- sector <- portfolio <- NULL 
  strategy <- substrategy <- gmv.temp <- pnl.temp <- contract.temp <- NULL
  
  ## Change the column names of date to "date" and use a temporary name "pnl.temp"
  ## for the P&L column. Return error if these two columns are missing. 
  
  if(pnl.var %in% colnames(x)){
    colnames(x)[which(colnames(x) == pnl.var)]  <- "pnl.temp"
  } 
  else {
    stop("Did not find P&L, or the column name for P&L is wrong")
  }
  
  if(date.var %in% colnames(x)){
    colnames(x)[which(colnames(x) == date.var)] <- "date"
  } else{
    stop("Did not find dates, or the column name for dates is wrong")
  }
  
  ## Change the column names of sectors, portfolios and strategies to the ones
  ## that are standard to this package, if these columns exist. If not, these
  ## columns will be created to fit into the later functions, and the elements
  ## of all the created entries will be "".
  
  if(sector.var %in% colnames(x)){
    colnames(x)[which(colnames(x) == sector.var)] <- "sector"
  } else {
    x <- mutate(x, sector = "")
  }
  
  if(portfolio.var %in% colnames(x)){
    colnames(x)[which(colnames(x) == portfolio.var)] <- "portfolio"
  } else {
    x <- mutate(x, portfolio = "")
  }
  
  if(strategy.var %in% colnames(x)){
    colnames(x)[which(colnames(x) == strategy.var)] <- "strategy"
  } else {
    x <- mutate(x, strategy = "")
  }
  
  if(substrategy.var %in% colnames(x)){
    colnames(x)[which(colnames(x) == substrategy.var)] <- "substrategy"
  } else {
    x <- mutate(x, substrategy = "")
  }
  
  ## If the input data frame does not contain instrument ID's, take instrument 
  ## names as the instruments' ID's. Otherwise, change the column name of ID's
  ## into the standard name
  
  if(id.var %in% colnames(x)){
    colnames(x)[which(colnames(x) == id.var)] <- "id"
  } else if(name.var %in% colnames(x)){
    colnames(x)[which(colnames(x) == name.var)] <- "id"
  } else {
    stop("Did not find names nor ID's")
  }
  
  ## Fill in the name column with instrument ID's if there does not exist a name
  ## column. Otherwise, change the column name of names into the standard name
  
  if(name.var %in% colnames(x)){
    colnames(x)[which(colnames(x) == name.var)] <- "name"
  } else {
    x <- mutate(x, name = "")
  }
  
  ## Fill in contract/net market value with one or the other's values, if one of 
  ## them is missing. If everything is good, just change the column names to
  ## the standard ones. If both number of contracts and net market values are
  ## missing, the function returns error.
  
  if(nmv.var %in% colnames(x)){
    colnames(x)[which(colnames(x) == nmv.var)] <- "nmv.temp"
    if(contract.var %in% colnames(x)){
      colnames(x)[which(colnames(x) == contract.var)] <- "contract.temp"
    }
    else {
      x <- mutate(x, contract.temp = abs(nmv.temp))
    }
  } 
  else if(contract.var %in% colnames(x)){
    colnames(x)[which(colnames(x) == contract.var)] <- "contract.temp"
    x <- mutate(x, nmv.temp = contract.temp)
  } 
  else {
    stop("Did not find net market values nor number of contracts")
  }
  
  ## Fill in the gross market value column if that does not exist. Otherwise just
  ## rename the column.
  
  if(gmv.var %in% colnames(x)){
    colnames(x)[which(colnames(x) == gmv.var)] <- "gmv.temp"
  } 
  else {
    x <- mutate(x, gmv.temp = abs(nmv.temp))
  }
  
  ## Summarize the data for individual commodities by adding up the data
  ## of different sub-ID's under the same commodity and ID, and then only take
  ## out the ones that are used by the backtestGraphics package.
  
  x <- x %>% 
    group_by(name, id, date, sector, portfolio, strategy, substrategy) %>%
    summarise(gmv = sum(gmv.temp, na.rm = TRUE),
              nmv = sum(nmv.temp, na.rm = TRUE),
              pnl = sum(pnl.temp, na.rm = TRUE),
              contract = sum(abs(contract.temp), na.rm = TRUE)) %>% 
    ungroup()
  
  return(x)
}