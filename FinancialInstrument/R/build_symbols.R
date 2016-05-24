#' construct a series of symbols based on root symbol and suffix letters
#'
#' The columns needed by this version of the function are \code{primary_id}
#' and \code{month_cycle}. \code{primary_id} should match the \code{primary_id}
#' of the instrument describing the root contract.
#' \code{month_cycle} should contain a comma delimited string describing the
#' month sequence to use, e.g. \code{"F,G,H,J,K,M,N,Q,U,V,X,Z"} for all months
#' using the standard futures letters, or \code{"H,M,U,Z"} for quarters, or
#' \code{"Mar,Jun,Sep,Dec"} for quarters as three-letter month abbreviations, etc.
#' The correct values will vary based on your data source.
#'
#' TODO add more flexibility in input formats for \code{roots}
#' #' @param yearlist vector of year suffixes to be applied, see Details
#' @param roots data.frame containing at least columns \code{primary_id} and \code{month_cycle}, see Details
#' @author Brian G. Peterson
#' @seealso \code{\link{load.instruments}}
#' @export
build_series_symbols <- function(roots, yearlist=c(0,1)) {
        symbols<-''
        id_col<-grep('primary_id',colnames(roots)) #TODO: check length
	date_col<-grep('month_cycle',colnames(roots)) #TODO: check length
	for (year_code in yearlist){
		for(i in 1:nrow(roots)) {
			symbols <- c(symbols, paste(paste(roots[i,id_col], strsplit(as.character(roots[i,date_col]),",")[[1]],sep=''),year_code,sep=''))
		}
	}
	return(symbols[-1])
}

#' build symbols for exchange guaranteed (calendar) spreads
#'
#' The columns needed by this version of the function are \code{primary_id},
#'  \code{month_cycle}, and code \code{contracts_ahead}.
#'
#' \code{primary_id} should match the \code{primary_id}
#' of the instrument describing the root contract.
#'
#' \code{month_cycle} should contain a comma delimited string describing the
#' month sequence to use, e.g. \code{"F,G,H,J,K,M,N,Q,U,V,X,Z"} for all months
#' using the standard futures letters, or \code{"H,M,U,Z"} for quarters, or
#' \code{"Mar,Jun,Sep,Dec"} for quarters as three-letter month abbreviations, etc.
#' The correct values will vary based on your data source.
#'
#' \code{contracts_ahead} should contain a comma-delimited string describing
#' the cycle on which the guaranteed calendar spreads are to be consructed,
#' e.g. '1' for one-month spreads, '1,3' for one and three month spreads,
#' '1,6,12' for 1, 6, and 12 month spreads, etc.  
#' For quarterly symbols, the correct \code{contracts_ahead} may be 
#' something like '1,2,3' for quarterly, bi-annual, and annual spreads.  
#'
#' \code{active_months} is a numeric field indicating how many months including  
#' the month of the \code{start_date} the contract is available to trade.  
#' This number will be used as the upper limit for symbol generation.
#' 
#' If \code{type} is also specified, it should be a specific instrument type, 
#' e.g. 'future_series','option_series','guaranteed_spread' or 'calendar_spread' 
#'
#' One of \code{data} or \code{file} must be populated for input data.
#'
#' @param data data.frame containing at least columns \code{primary_id}, \code{month_cycle}, amd \code{contracts_ahead}, see Details
#' @param file if not NULL, will read input data from the file named by this argument, in the same folrmat as \code{data}, above
#' @param outputfile if not NULL, will write out put to this file as a CSV
#' @param start_date date to start building from, of type \code{Date} or an ISO-8601 date string, defaults to \code{\link{Sys.Date}}
#' @author Ilya Kipnis <Ilya.Kipnis<at>gmail.com>
#' @seealso
#' \code{\link{load.instruments}}
#' \code{\link{build_series_symbols}}
# @examples 
# build_spread_symbols(data=data.frame(primary_id='CL',
#                                      month_sequence="F,G,H,J,K,M,N,Q,U,V,X,Z",
#                                      contracts_ahead="1,2,3",
#                                      type='calendar_spread')
#' @export
build_spread_symbols <- function(data=NULL,file=NULL,outputfile=NULL,start_date=Sys.Date())
{
        if(!is.null(data)) {
                Data<-data
        } else if(!is.null(file)) {
                Data<-read.csv(file,header=TRUE,stringsAsFactors=FALSE)
        } else {
                stop("you must either pass a data.frame as the 'data' parameter or pass the 'file' parameter")
        }



	yearCheck<-function(monthNum,yearNum){
		if(monthNum>12){
			yearNum=yearNum+1
		} else if(monthNum<0){
			yearNum=yearNum-1
		} else {
			yearNum=yearNum
		}
		return(yearNum)
	}

	makeNameCal<-function(primary_id,MonthOne,YearOne,MonthTwo,YearTwo){
		contractName<-NULL
		contractName<-paste(primary_id,MonthOne,YearOne,"-",MonthTwo,YearTwo,sep="")
		return(contractName)
	}

	makeCalRows<-function(dataRow){
		calFrame<-NULL
		calContracts<-as.numeric(unlist(strsplit(dataRow$contracts_ahead[1],",",fixed=TRUE)))
		contractIndex<-c(1:length(monthsTraded))
		contractTable<-cbind(monthsTraded,contractIndex)
		for(k in 1:length(calContracts)){
            if(length(monthsTraded)==12 & calContracts[k]>1){
                yearsAhead = trunc((calContracts[k]-1)/length(monthsTraded))%%10
                monthContractsAhead=(calContracts[k]-1)%%length(monthsTraded)
            } else {
                yearsAhead=trunc(calContracts[k]/length(monthsTraded))%%10
                monthContractsAhead=calContracts[k]%%length(monthsTraded)
            }
			workingContractNum<-as.numeric(contractTable[which(contractTable[,1]==workingContractMonthLet),2])
			monthIndex=workingContractNum+monthContractsAhead
			if(monthIndex>length(monthsTraded)){
				yearsAhead<-(yearsAhead+1)%%10
				monthIndex=monthIndex-length(monthsTraded)
			}
                        newCalMonthLetter<-contractTable[monthIndex,1]
                        newCalYearNumber<-workingContractYearNum+yearsAhead%%10
                        contractName<-makeNameCal(Data$primary_id[i],workingContractMonthLet,workingContractYearNum,newCalMonthLetter,newCalYearNumber)
                        calRow<-cbind(contractName,dataRow$type[1])
                        #Other data would be added here.  Take in a row from root contracts, add other details (EG timezone, currency, etc...)
                        calFrame<-rbind(calFrame,calRow)
                }
                return(calFrame)
        }

	PPFNames<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
	PPFNums<-c(-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
        PPFLets<-c("F","G","H","J","K","M","N","Q","U","V","X","Z","F","G","H","J","K","M","N","Q","U","V","X","Z","F","G","H","J","K","M","N","Q","U","V","X","Z")
        PPFFrame<-as.data.frame(cbind(PPFNames,PPFNums,PPFLets),stringsAsFactors=FALSE)  #PPF=Past, Present, Future

        contractFrame<-NULL
        today<-start_date
        for(i in 1:nrow(Data)){
		currentMonthNum<-as.numeric(substr(today,6,7))
		monthsTraded<-unlist(strsplit(Data$month_cycle[i],",",fixed=TRUE))
		monthNumsTraded<-as.numeric(PPFFrame[which(PPFFrame[,3] %in% monthsTraded),2])
		currentContractMonthNum<-as.numeric(monthNumsTraded[(min(which(monthNumsTraded>=currentMonthNum)))])
		currentContractMonthLet<-PPFFrame[currentContractMonthNum,3]
		currentYearNum<-as.numeric(substr(today,4,4))
		currentContractYearNum<-yearCheck(currentContractMonthNum,currentYearNum)

		workingContractMonthNum<-currentContractMonthNum
		workingContractMonthLet<-currentContractMonthLet
		workingContractYearNum<-currentContractYearNum

		contractIndex<-c(1:length(monthsTraded))
		contractTable<-cbind(monthsTraded,contractIndex)

		if(Data$type[i]=="calendar_spread" || Data$type[i]=='guaranteed_spread'){
			calFrame<-makeCalRows(Data[i,])
			contractFrame<-rbind(contractFrame,calFrame)
		} else {
			#make name for non-calendar
			#add other details to the row
			#rbind it to the contractFrame
		}

        for(j in 1:(Data$active_months[i]-1)){
            yearsAhead=trunc(j/length(monthsTraded))%%10 #aka if it's 10+ years ahead, since the year is one digit, for contracts further out than 10 years, you'll get 0, 1, 2, 3 instead of 10, 11, 12, etc.
            monthContractsAhead=j%%length(monthsTraded)
            currentContractNum<-as.numeric(contractTable[which(contractTable[,1]==currentContractMonthLet),2])
            monthIndex=currentContractNum+monthContractsAhead
            if(monthIndex>length(monthsTraded)){
				yearsAhead<-(yearsAhead+1)%%10
				monthIndex=monthIndex-length(monthsTraded)
			}
			workingContractMonthLet<-contractTable[monthIndex,1]
			workingContractYearNum<-currentContractYearNum+yearsAhead
			if(Data$type[i]=="calendar_spread" || Data$type[i]=='guaranteed_spread'){
				calFrame<-makeCalRows(Data[i,])
				contractFrame<-rbind(contractFrame,calFrame)
			} else {
				#make name for non-calendar
				#add other details to the row
				#rbind it to the contractFrame
			}
		}
    }

    rownames(contractFrame)<-NULL
    colnames(contractFrame)<-c("symbol","type")

    if(!is.null(outputfile)){
            write.csv(contractFrame,outputfile) 
    } else {
            return(contractFrame)
    }
}

###############################################################################
# R (http://r-project.org/) Instrument Class Model
#
# Copyright (c) 2009-2012
# Peter Carl, Dirk Eddelbuettel, Jeffrey Ryan, 
# Joshua Ulrich, Brian G. Peterson, and Garrett See
#
# This code is distributed under the terms of the GNU Public License (GPL)
# for full details see the file COPYING
#
# $Id: build_symbols.R 1638 2014-10-08 03:43:16Z gsee $
#
###############################################################################

