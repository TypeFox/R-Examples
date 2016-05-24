modelgen <- function(){
  mmtfaModels <- list()
  mmtfaModels[["modold"]]  <- c("UUUU", "UUUC", "UCUU", "UCUC", "UUCU", "UUCC", "UCCU", "UCCC", 
                            "CUUU", "CUUC", "CCUU", "CCUC", "CUCU", "CUCC", "CCCU", "CCCC",
                            "Mt1U", "Mt1C", "Mt2U", "Mt2C", "Mt3U", "Mt3C", "Mt4U", "Mt4C")
  mmtfaModels[["allmodels"]] <- c("UUUU", "UUUC", "UCCU","UCCC", "UUIU", "UUIC", "UCIU", "UCIC", 
                           "CUUU", "CUUC", "CCCU", "CCCC", "CUIU", "CUIC", "CCIU", "CCIC",
                           "CUCU","CUCC","UUCU","UUCC","UCUU","UCUC","CCUU","CCUC")
  mmtfaModels[["dfconstrained"]] <- c("UUUC", "UCCC", "UUIC", "UCIC", 
                                      "CUUC",  "CCCC", "CUIC", "CCIC",
                                      "CUCC", "UUCC","UCUC", "CCUC")
  mmtfaModels[["dfunconstrained"]] <- c("UUUU", "UCCU", "UUIU", "UCIU", 
                                         "CUUU", "CCCU", "CUIU", "CCIU",
                                         "CUCU", "UUCU", "UCUU", "CCUU")
  mmtfaModels
}
