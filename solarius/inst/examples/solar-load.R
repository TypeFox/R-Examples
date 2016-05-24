# @ http://helix.nih.gov/Documentation/solar-6.6.2-doc/91.appendix_1_text.html#load
#solar> help file-pedigree                                                       
#
# The pedigree file consists of one record for each individual in the data
# set.  Each record must include the following fields:
#
#    ego ID, father ID, mother ID, sex
#
# In addition, a family ID is required when ego IDs are not unique across
# the entire data set.  If the data set contains genetically identical
# individuals, an MZ-twin ID must be present (as described below).  If an
# analysis of household effects is planned, a household ID can be included
# (also described below).
#
# The default field names are ID, FA, MO, SEX, FAMID, MZTWIN, and HHID.
#

### par
solar.dir <- "solar"

### dir create
dir.create(solar.dir, showWarnings = FALSE)

### data
dat <- loadMulticData()
dat30 <- subset(dat, famid < 30)

dat30 <- rename(dat30, c("famid" = "FAMID", "id" = "ID", "fa" = "FA", "mo" = "MO", "sex" = "SEX"))

# ped
ped.cols <- c("FAMID", "ID", "FA", "MO", "SEX")
ped30 <- subset(dat30, select = ped.cols)

write.table(ped30, file.path(solar.dir, "dat.ped"), 
  row.names = FALSE, sep = ",", quote = FALSE)
  
# phen
ind.ID <- which(ped.cols == "ID")
stopifnot(length(ind.ID) == 1)
ped.cols2 <- ped.cols[-ind.ID]

ind <- which(names(dat30) %in% ped.cols2)
phen30 <- dat30[, -ind]

write.table(phen30, file.path(solar.dir, "dat.phe"), 
  row.names = FALSE, sep = ",", quote = FALSE, na = "")
  
### run solar
wd <- getwd()
setwd(solar.dir)

ret <- try({
  system("solar", input = c("load pedigree dat.ped", 
    "load phenotypes dat.phe",
    "trait trait1",
    "covariate age age^2 SEX trait2",
    "polygenic",
    "quit"))
})

setwd(wd)  
