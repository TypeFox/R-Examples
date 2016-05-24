if(version$language == "R") {
## BEGIN R SPECIFIC IMPLEMENTATION

json.parse <- function(string) fromJSON(string)

## END R SPECIFIC IMPLEMENTATION
}

## TERR IMPLEMENTATION IN "terrUtils" PACKAGE

