#' Read Bob Schiller's Data
schiller <- function() {
    Date <- NULL

    #schiller_data <- gdata::read.xls('../data-raw/ie_data.xls',
    #                                 sheet='Data', as.is=TRUE)
    schiller_data <- gdata::read.xls('http://www.econ.yale.edu/~shiller/data/ie_data.xls',
                                     sheet='Data', as.is=TRUE)

    sp_data     <- SP500()
    sp_data     <- dplyr::arrange(sp_data, Date)
    sp_data$chg <- calc_chg(sp_data$Adj.Close)

    # extract the monthly dividend data from Schiller
    # ===============================================
    # starting on January 1950
    idx_1950.01 <- which(schiller_data[,1]=="1950.01")
    # ending on the last entry afterwards with dividend data
    idx_endofdiv <- idx_1950.01 +
        min(which(schiller_data[idx_1950.01:nrow(schiller_data),3]==''))-2

    # data.frame with numeric dividend data
    monthly        <- schiller_data[idx_1950.01:idx_endofdiv, 1:3]
    names(monthly) <- c("dte","comp","div")
    monthly$comp   <- as.numeric(monthly$comp)
    monthly$div    <- as.numeric(monthly$div)

    # add monthly div data to sp500 daily data
    # ========================================
    sp_data$Year  <- lubridate::year(sp_data$Date)
    sp_data$Month <- lubridate::month(sp_data$Date)

    year_month_idx <- c(1, 1+which(diff(sp_data$Month)!=0)) # 1st of each month
    sp_data$div <- 0
    sp_data$div[year_month_idx] <- c(monthly$div,
                                    rep(NA,nrow(sp_data[year_month_idx,])-length(monthly$div)))
    sp_data$schiller_dte[year_month_idx] <- c(monthly$dte,
                                    rep(NA,nrow(sp_data[year_month_idx,])-length(monthly$div)))

    # in case NA is introduced into div column
    sp_data$div[is.na(sp_data$div)] <- 0


    # from https://dqydj.com/sp-500-return-calculator/
    # To calculate the 'dividend reinvested' price index, I take the trailing
    # twelve month dividend yield reported in any month of Shiller's data and
    # divide by 12 to get an approximate count of dividends paid out in the month.
    # Using that number, I calculate how many 'shares' of the S&P 500 index
    # I can buy, and run a cumulative count from 1876 to 2012.


    # add dividend payout ratio (dpr)
    sp_data$dpr <- 0
    sp_data$dpr[year_month_idx] <- sapply(year_month_idx,
        function(idx) {
            comp <- sp_data[idx,'Adj.Close']
            div  <- sp_data[idx,'div']
            return(
                ifelse(sp_data[idx,'Month']==12,(div/comp)/12,0)
            )
    })

    sp_data$cumsum <- cumsum(sp_data$dpr)

    #first_of_month <- sp_data[year_month_idx,]

    sp_data$tr <- sp_data$Adj.Close[1]
    for (idx in 2:nrow(sp_data)) {
        sp_data$tr[idx] <- sp_data$tr[idx-1]*sp_data$chg[idx] #+
            #ifelse(sp_data$Month[idx]==12, sp_data$div[idx]/60, 0)
        #ifelse(sp_data$Month[idx]==12, sp_data$div[idx]/
         #          (12+((289.2027-249.76)/249.76*12)), 0)
    }
    foo=sp_data[sp_data$Year==1988,]
}
