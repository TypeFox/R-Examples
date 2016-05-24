
###################################################
### code chunk: Ch01init
###################################################
options(width=65, digits=5, show.signif.stars = FALSE )
date()
packageVersion("nlmeU")
packageVersion("nlme")
sessionInfo()


###################################################
### code chunk: R2.1
###################################################
## dataDir <- file.path("C:", "temp")                # Data directory
dataDir <- file.path(.Library, "nlmeU", "csvData")   # Directory in nlmeU package
fp   <- file.path(dataDir, "armd240.data.csv")       # .csv file path
armd240.data <- read.csv(fp, header = TRUE)          # Loading csv data
dim(armd240.data)                                    # No. of rows and cols
(nms <- names(armd240.data))                         # Variables' names
unique(sapply(armd240.data, class))                  # Variables' classes
str(armd240.data)                                    # Data structure
names(armd240.data) <- abbreviate(nms)               # Variables' names shortened
head(armd240.data, 3)                                # First 3 records
names(armd240.data) <- nms                           # Variables' names reinstated


###################################################
### code chunk: R2.2
###################################################
data(armd.wide, package = "nlmeU")                   # armd.wide data loaded
str(armd.wide)                                       # Structure of data
head(armd.wide)                                      # First few records
(facs <- sapply(armd.wide, is.factor))               # Factors indicated 
names(facs[facs == TRUE])                            # Factor names displayed                  


###################################################
### code chunk: R2.3a
###################################################
attach(armd240.data)                                 # Attach armd240.data
treat.f <- factor(treat,                             # Factor created 
  labels = c("Placebo", "Active")                    # 1 -> Placebo, 2 -> Active
) 
levels(treat.f)
str(treat.f)


###################################################
### code chunk: R2.3b
###################################################
miss.pat <-                                # missing patterns
    nlmeU:::missPat(visual4, visual12, visual24, visual52)  
length(miss.pat)                           # Vector length
mode(miss.pat)                             # Vector mode
miss.pat                                   # Vector contents 
detach(armd240.data)                       # Detach armd240.data


###################################################
### code chunk: R2.4
###################################################
data(armd0, package = "nlmeU")             # From nlmeU package
dim(armd0)                                 # No. of rows and cols
head(armd0)                                # First six records
names(armd0)                               # Variables' names
str(armd0)                                 # Data structure    


###################################################
### code chunk: R2.5
###################################################
auxDt <- subset(armd0, time > 0)          # Post-baseline measures
dim(auxDt)                                # No. of rows & cols
levels(auxDt$time.f)                      # Levels of treat.f
armd <- droplevels(auxDt)                 # Drop unused levels   
levels(armd$time.f)                       # Baseline level dropped
armd <- within(armd, {                    # Contrasts assigned
   contrasts(time.f) <- contr.poly(4, scores = c(4, 12, 24, 52))
})
head(armd)                                # First six records


###################################################
### code chunk: R2.6a
###################################################
fp <- file.path(dataDir, "prt.subjects.data.csv")
prt.subjects.data <- read.csv(fp, header = TRUE, as.is = TRUE) 
dim(prt.subjects.data)
names(prt.subjects.data)
str(prt.subjects.data)
head(prt.subjects.data, 4)


###################################################
### code chunk: R2.6b
###################################################
fp <- file.path(dataDir, "prt.fiber.data.csv")
prt.fiber.data <- read.csv(fp, header = TRUE) 
str(prt.fiber.data)
head(prt.fiber.data, 4)


###################################################
### code chunk: R2.7a
###################################################
prt.subjects <- within(prt.subjects.data,{
   id    <- factor(id)
   bmi   <- weight/(height^2)
   sex.f <- factor(gender, labels = c("Female", "Male"))
   age.f <- factor(ageGrp, labels = c("Young", "Old"))
   prt.f <- factor(trainGrp, levels = c("1", "0"), 
             labels = c("High", "Low"))
   gender<- ageGrp <- trainGrp <- height <- weight <- NULL
})
str(prt.subjects)


###################################################
### code chunk: R2.7b
###################################################
prt.fiber  <- within(prt.fiber.data, {
   id      <- factor(id)
   fiber.f <- factor(fiber.type, 
               labels = c("Type 1", "Type 2"))
   occ.f   <- factor(train.pre.pos, 
               labels = c("Pre", "Pos"))     
   fiber.type <- train.pre.pos <- NULL
})
str(prt.fiber)


###################################################
### code chunk: R2.8
###################################################
prt <- merge(prt.subjects, prt.fiber, sort=FALSE)
dim(prt)
names(prt)
head(prt)


###################################################
### code chunk: R2.9
###################################################
data(classroom, package = "WWGbook")
dim(classroom)                       # Number of rows & variables
names(classroom)                     # Variable names
# classroom                          # Raw data (long output)
str(classroom)


###################################################
### code chunk: R2.10
###################################################
SIIdata <- within(classroom, {
   sex <- factor(sex, 
        levels = c(0, 1),
        labels = c("M", "F"))             # 0 -> 1(M),  1 -> 2(F) 
   minority <- factor(minority, 
        labels = c("Mnrt:No", "Mnrt:Yes"))# 0 -> 1(No), 1 -> 2(Yes)
   schoolid <- factor(schoolid)         
   classid  <- factor(classid)          
   childid  <- factor(childid)
})
str(SIIdata)


###################################################
### code chunk: R2.11
###################################################
rdaDir <- file.path("C:", "temp")         # Dir path
fp <- file.path(rdaDir, "SIIdata.Rdata")  # External file path
save(SIIdata, file = fp)                  # Save data
file.exists(fp)                           
(load(fp))                                # Load data


###################################################
### code chunk: R2.12
###################################################
data(SIIdata, package = "nlmeU")         # Load data
dtId <- subset(SIIdata, select = c(schoolid, classid, childid))
names(dtId)                              # id names
any(duplicated(dtId))                    # Any duplicate ids?
require(nlme)
names(gsummary(dtId, form = ~childid, inv = TRUE))
names(gsummary(dtId, form = ~classid, inv = TRUE))
names(gsummary(dtId, form = ~schoolid, inv = TRUE))


###################################################
### code chunk: R2.13a
###################################################
(nms1 <- names(gsummary(SIIdata,
    form = ~schoolid,              # schoolid-specific
    inv = TRUE))) 


###################################################
### code chunk: R2.13b
###################################################
nms2a <- names(gsummary(SIIdata, 
    form = ~classid,               # classid- and schoolid-specific 
    inv = TRUE)) 
idx1  <- match(nms1, nms2a)
(nms2 <- nms2a[-idx1])             # classid-specific


###################################################
### code chunk: R2.13c
###################################################
nms3a <- names(gsummary(SIIdata, 
                 form = ~childid,  # All
                 inv = TRUE)
)
idx12 <- match(c(nms1, nms2), nms3a)
nms3a[-idx12]                      # childid-specific


###################################################
### code chunk: R2.14
###################################################
fp <- file.path(dataDir, "crossreg.data.csv")
crossreg.data <- read.csv(fp, header = TRUE) 
dim(crossreg.data)                      # No. of rows and columns
names(crossreg.data)                    # Variable names
head(crossreg.data)                     # First six records
str(crossreg.data)                      # Data structure


###################################################
### code chunk: R2.15
###################################################
unique(crossreg.data$target)          # Unique values for target
(unique(crossreg.data$id))            # Unique values for id 
unique(crossreg.data$scorec)          # Unique values for scorec
summary(crossreg.data$scorec)         # Summary statistics for scorec


###################################################
### code chunk: R2.16a
###################################################
nItms  <- c(4, 6, 8, 5, 9, 6, 8, 6, 5)        # See Table 2.1
(lbls <- paste("T", 1:9, "(", nItms, ")", sep = ""))
fcat <- within(crossreg.data, {
        id <- factor(id)
        target <- factor(target, labels = lbls)
})
str(fcat)


###################################################
### code chunk: R2.16b
###################################################
tab1 <- xtabs(~ id + target, data = fcat)     # id by target table
tab1[c(1,2,539),]
all(tab1 > 0)                                 # All counts > 0? 
range(tab1)                                   # Range of counts


###################################################
### code chunk: Cleanup
###################################################
rm(armd240.data, armd.wide, armd0, armd)           # Data not needed
rm(prt.fiber.data, prt.subjects.data, prt.fiber, prt.subjects, prt)
rm(classroom, SIIdata) 
rm(crossreg.data, fcat) 
rm(dataDir, rdaDir, auxDt)
rm(treat.f, miss.pat, facs, lbls)
rm(nms,nms1, nms2, nms2a, nms3a)
rm(idx1, idx12, nItms, tab1, dtId)

### sessionInfo
  
sessionInfo()            # with packages attached
detach(package:nlme)
