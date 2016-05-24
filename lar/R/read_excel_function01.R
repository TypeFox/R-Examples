# File:   read_excel_function01.R
# Date:   March 25, 2014
# Author: richard.zijdeman@iisg.nl
# Last change: 
  # 2014-04-29, cleaned up file for CRAN


# first load data with labels internally
# data(lar_classification_treemap08, envir = environment())

# defining labour.rel.1 globally to avoid NOTE for CRAN submission
if(getRversion() >= "2.15.1") utils::globalVariables("labour.rel.1")

# read in lar excel-input-file and provide proper format for treemap
read.lar <- function(file) {
  df <- read.xlsx(file, 
                  sheetName  = "Data entry",
                  encoding   = "UTF-8",
                  colClasses = c(rep("integer", 7), rep("character", 4), 
                                 rep("integer", 3), rep("character", 3),
                                 rep("integer", 2), rep("character", 11),
                                 rep("integer", 2), "character",
                                 rep(c("integer", "numeric"), 3),
                                 rep("character",5)))
  dt  <- data.table(df)
  
  ## issue 1: delete row 1-2 (auxiliairy column headers)
  dt1 <- dt[3:nrow(dt)]
  
  ## issue 2: delete auxiliary columns
  dt2 <- dt1[, 1:40, with = FALSE]
  
  ## issue 3: remove completely empty rows
  # code borrowed from kmiddleton's github page
  numNA  <- apply(dt2, 1, function(x) sum(is.na(x)))
  numvar <- dim(dt2)[2]
  dt3    <- dt2[(numNA < numvar),]
  
  ## issue 4: change labels
  # ignore first 2 chars
  setnames(dt3, c(names(dt3)), c(substring(names(dt3), 3))) # read from 3rd char
  
  # replace ".." with "pct" ## (gsub sees "." as wildcard therefore use: \\)
  setnames(dt3, c(names(dt3)), c(gsub("\\..", "pct", names(dt3)))) #
  
  # remove "." ## (gsub sees "." as wildcard therefore use: \\)
  setnames(dt3, c(names(dt3)), c(gsub("\\.", "", names(dt3)))) #
  
  # replace "_" by "." (google R-style guide convention)
  setnames(dt3, c(names(dt3)), c(gsub("_", "\\.", names(dt3)))) #
  
  # replace "0x" by "x"
  setnames(dt3, c(names(dt3)), c(gsub("0", "", names(dt3)))) #
  
  # convert to lower cases
  setnames(dt3, c(names(dt3)), c(tolower(names(dt3))))
  
  ## issue 5: removing rows not holding info on the first labour relation
  ## this removes data on population and 'negative' total numbers
  dt4 <- dt3[!(is.na(labour.rel.1))]
  
  ## always convert last dt to 'lar', for follow up functions
  lar <- as.data.frame(dt4)
  


#### Add labels and vars for treemaps ####

#labels <- lar_classification_treemap08

# merge data with labour relation "lar", with labour relation labels "labels"
# this process is repeated for the multiple labour relations

# merge labels with 1st labour relation
lar1 <- merge(lar, labels,
              by.x  = "labour.rel.1", 
              by.y  = "lr.cat",
              all.x = TRUE)

names(lar1) <- gsub("lr.lab", "lab1", names(lar1))
names(lar1) <- gsub("lr.txt", "txt1", names(lar1))
names(lar1) <- gsub("lr.no",  "no1" , names(lar1))

# merge labels with 2nd labour relation
lar2 <- merge(lar1, labels,
              by.x  = "labour.rel.2", 
              by.y  = "lr.cat",
              all.x = TRUE)
names(lar2) <- gsub("lr.lab", "lab2", names(lar2))
names(lar2) <- gsub("lr.txt", "txt2", names(lar2))
names(lar2) <- gsub("lr.no",  "no2" , names(lar2))

# merge labels with 3rd labour relation
lar3 <- merge(lar2,labels,
              by.x  = "labour.rel.3", by.y = "lr.cat",
              all.x = TRUE)
names(lar3) <- gsub("lr.lab", "lab3", names(lar3))
names(lar3) <- gsub("lr.txt", "txt3", names(lar3))
names(lar3) <- gsub("lr.no",  "no3",  names(lar3))


## add numbers  2: and 3: to indicate 2nd and 3rd labour relation
lar3$txt2.1 <- ifelse(is.na(lar3$txt2.1), 
                      lar3$txt2.1,  paste0("+",lar3$txt2.1))
lar3$txt3.1 <- ifelse(is.na(lar3$txt3.1), 
                      lar3$txt3.1,  paste0("++",lar3$txt3.1))


## now create a sorting order variable (required for treemap)
lar3$sortID2 <- NA
lar3$sortID2 <- ifelse(lar3$txt1.1 == "Non working",
                            1, lar3$sortID2)
lar3$sortID2 <- ifelse(lar3$txt1.1 == "Reciprocal",
                            2, lar3$sortID2)
lar3$sortID2 <- ifelse(lar3$txt1.1 == "Tributary",
                            3, lar3$sortID2)
lar3$sortID2 <- ifelse(lar3$txt1.1 == "Commodified",
                            4, lar3$sortID2)
lar3$sortID2 <- ifelse(lar3$txt1.1 == "EitherOr",
                            5, lar3$sortID2)
lar3$sortID2 <- ifelse(lar3$txt1.1 == "Unknown",
                            6, lar3$sortID2)

## add a new variable to easily subset country and time periods
lar3$ctry.time <- with(lar3, paste(country, year, sep = "."))

lar3$bmyear <- lar3$year
lar3$bmyear <- ifelse(lar3$year >= 1500 & 
                             lar3$year < 1575, 1500, lar3$bmyear)
lar3$bmyear <- ifelse(lar3$year >= 1575 & 
                             lar3$year < 1725, 1650, lar3$bmyear)
lar3$bmyear <- ifelse(lar3$year >= 1725 & 
                             lar3$year < 1850, 1800, lar3$bmyear)
lar3$bmyear <- ifelse(lar3$year >= 1850 & 
                             lar3$year < 1950, 1900, lar3$bmyear)
lar3$bmyear <- ifelse(lar3$year >= 1950 & 
                             lar3$year < 2050, 2000, lar3$bmyear)

## output: 
return(lar3)
}

# EOF