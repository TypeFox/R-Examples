#' Summary Method
#' 
#' \code{summary} method displays only the essential info about the CDS
#' class object.
#' 
#' @name summary
#' 
#' @aliases summary,CDS-method
#' @docType methods
#' @rdname summary-methods
#' @param object the input \code{CDS} object
#' @param ... additional arguments to pass in
#' 
#' @export

setMethod("summary",
          signature(object = "CDS"),
          function(object,
                   ...){
            
            endDate <- add_dates(data.frame(date = as.Date(object@date),
                                            tenor = object@tenor,
                                            currency = object@currency))$endDate
            
            cat(sprintf(paste("Contract Type:", object@contract,
                              sep = paste0(rep(" ",
                                40-nchar(as.character("Contract Type:")) -
                                   nchar(as.character(object@contract))), collapse = ""))),
                sprintf(paste("   Date:", object@date,
                              sep = paste0(rep(" ", 
                                40-nchar(as.character("   Date:")) -
                                   nchar(as.character(object@date))), collapse = ""))), "\n",
                
                
                sprintf(paste("Entity Name:", object@name,
                              sep = paste0(rep(" ",
                                40-nchar(as.character("Entity Name:")) -
                                   nchar(as.character(object@name))), collapse = ""))),
                sprintf(paste("   RED:", object@RED,
                              sep = paste0(rep(" ", 
                                40-nchar(as.character("   RED:")) -
                                   nchar(as.character(object@RED))), collapse = ""))), "\n",
                
                sprintf(paste("Currency:", object@currency,
                              sep = paste0(rep(" ",
                                40-nchar(as.character("Currency:")) -
                                   nchar(as.character(object@currency))), collapse = ""))),
                sprintf(paste("   End Date:", endDate,
                              sep = paste0(rep(" ", 
                                40-nchar(as.character("   End Date:")) -
                                   nchar(as.character(endDate))), collapse = ""))), "\n",
                sprintf(paste("Spread:",
                              format(round(object@spread, 4),big.mark = ",",
                                     scientific=F),
                              sep = paste0(rep(" ",
                                40-nchar(as.character("Spread:")) -
                                   nchar(as.character(
                                   format(round(object@spread, 4), big.mark=",",
                                                          scientific=F)))),
                                           collapse = ""))),
                sprintf(paste("   Coupon:",
                              format(object@coupon, big.mark = ",",
                                     scientific=F),
                              sep = paste0(rep(" ", 
                                40-nchar(as.character("   Coupon:")) -
                                   nchar(as.character(format(object@coupon, big.mark=",",
                                                              scientific=F)))),
                                           collapse = ""))), "\n",
                
                sprintf(paste("Upfront:",
                              format(round(object@upfront, 0),big.mark=",",
                                     scientific =F ),
                              sep = paste0(rep(" ",
                                40-nchar(as.character("Upfront:")) -
                                   nchar(as.character(
                                      format(round(object@upfront, 0),
                                                          big.mark=",", scientific=F)))),
                                           collapse = ""))),
                sprintf(paste("   Spread DV01:",
                              format(round(object@spread.DV01, 0), big.mark = ",",
                                     scientific=F),
                              sep = paste0(rep(" ", 
                                40-nchar(as.character("   Spread DV01:")) -
                                   nchar(as.character(
                                      format(round(object@spread.DV01, 0), big.mark=",",
                                                scientific=F)))),
                                           collapse = ""))), "\n",
                
                sprintf(paste("IR DV01:",
                              format(round(object@IR.DV01, 2),big.mark=",",
                                     scientific=F),
                              sep = paste0(rep(" ",
                                40-nchar(as.character("IR DV01:")) -
                                   nchar(as.character(
                                    format(round(object@IR.DV01, 2), big.mark=",",
                                                          scientific=F)))),
                                           collapse = ""))),
                sprintf(paste("   Rec Risk (1 pct):",
                              format(round(object@rec.risk.01, 2),big.mark=",",
                                     scientific=F),
                              sep = paste0(rep(" ", 
                                40-nchar(as.character("   Rec Risk (1 pct):")) -
                                   nchar(as.character(
                                     format(round(object@rec.risk.01, 2),
                                                          big.mark=",",
                                                          scientific=F)))),
                                           collapse = ""))), "\n",
                sep = ""
            )
            
            cat("\n")
          }
)
