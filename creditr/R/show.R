#' Show Method
#' 
#' \code{show} shows a CDS class object.
#' 
#' @name show
#' @aliases show,CDS-method
#' @docType methods
#' @rdname show-methods
#' @param object the input \code{CDS} object
#' @param ... additional arguments to pass in
#' 
#' @export

setMethod("show",
          signature(object = "CDS"),
          function(object){
            
            cat("CDS Contract \n")
            
            cat(sprintf(paste("Contract Type:", object@contract,
                        sep = paste0(rep(" ",
                        40-nchar(as.character("Contract Type:")) -
                        nchar(as.character(object@contract))), collapse = ""))),
                sprintf(paste("   Currency:", object@currency,
                        sep = paste0(rep(" ",
                        40-nchar(as.character("   Currency:")) -
                        nchar(as.character(object@currency))),
                        collapse = ""))), "\n",
                
                sprintf(paste("Entity Name:", object@name,
                        sep = paste0(rep(" ",
                        40-nchar(as.character("Entity Name:")) -
                        nchar(as.character(object@name))), collapse = ""))),
                sprintf(paste("   RED:", object@RED,
                              sep = paste0(rep(" ",
                                40-nchar(as.character("   RED:")) -
                                  nchar(as.character(object@RED))),
                                           collapse = ""))), "\n",
                
                sprintf(paste("date:", object@date,
                         sep = paste0(rep(" ", 
                                40-nchar(as.character("date:")) -
                         nchar(as.character(object@date))), 
                         collapse = ""))), "\n", sep = "")
            
            cat("\n")
            cat("Calculation \n")
            
            cat(sprintf(paste("price:", round(object@price, 2),
                         sep = paste0(rep(" ", 
                         40-nchar(as.character("price:")) -
                          nchar(as.character(
                                round(object@price, 2)))), 
                         collapse = ""))),
                
                sprintf(paste("   Spread:",
                    format(round(object@spread, 4), big.mark = ",",
                        scientific = F),
                    sep = paste0(rep(" ",
                      40-nchar(as.character("   spread:")) -
                         nchar(as.character(
                         format(round(object@spread, 4),
                                big.mark = ",", scientific = F)))),
                                      collapse = ""))), "\n",
                
                sprintf(paste("Principal:",
                      format(round(object@principal, 0), big.mark=",",
                                     scientific=F),
                        sep = paste0(rep(" ",
                                40-nchar(as.character("Principal:")) -
                                   nchar(as.character(
                              format(round(object@principal, 0), 
                                    big.mark = "F", scientific = F)))),
                                           collapse = ""))),
                sprintf(paste("   Spread DV01:",
                      format(round(object@spread.DV01, 0), big.mark=",",
                              scientific=F),
                      sep = paste0(rep(" ", 
                            40-nchar(as.character("   Spread DV01:")) -
                               nchar(as.character(
                              format(round(object@spread.DV01, 0), 
                                      big.mark=",", scientific=F)))),
                                  collapse = ""))), "\n",
                
                sprintf(paste("Accrual:",
                      format(round(object@accrual, 0), big.mark=",",
                                     scientific=F),
                      sep = paste0(rep(" ",
                          40-nchar(as.character("Accrual:")) -
                            nchar(as.character(
                                format(round(object@accrual, 0), 
                                    big.mark = ",", scientific=F)))),
                                  collapse = ""))), 
                sprintf(paste("   IR DV01:", format(
                      round(object@IR.DV01, 2),big.mark=",",
                                                    scientific=F),
                          sep = paste0(rep(" ",
                              40-nchar(as.character("   IR DV01:")) -
                                 nchar(as.character(
                                  format(round(object@IR.DV01, 2), 
                                    big.mark=",", scientific=F)))),
                                    collapse = ""))), "\n",
                
                sprintf(paste("Upfront:",
                    format(round(object@upfront, 0), big.mark = ",",
                                     scientific=F),
                    sep = paste0(rep(" ",
                          40-nchar(as.character("Upfront:")) -
                             nchar(as.character(
                          format(round(object@upfront, 0), big.mark=",",
                                        scientific=F)))),
                                           collapse = ""))), 
                sprintf(paste("   Rec Risk (1 pct):",
                    format(round(object@rec.risk.01, 2),big.mark=",",
                                     scientific=F),
                    sep = paste0(rep(" ", 
                        40-nchar(as.character("   Rec Risk (1 pct):")) -
                           nchar(as.character(
                              format(round(object@rec.risk.01, 2),
                                  big.mark=",", scientific=F)))),
                            collapse = ""))), "\n",
                
                sprintf(paste("Default Prob:", round(object@pd, 4),
                    sep = paste0(rep(" ",
                        40-nchar(as.character("Default Prob:")) -
                           nchar(as.character(round(object@pd, 4)))),
                                collapse = ""))), "\n",
                
                sep = ""
            )
            cat("\n")
          }
)