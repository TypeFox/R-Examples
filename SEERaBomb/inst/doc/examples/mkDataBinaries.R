rm(list=ls())
require(dplyr)
library(SEERaBomb)
# (df=getFields("~/data/SEER")) # this doesn't fly in windows since ~ maps to /users/radivot/documents
(df=getFields("/Users/radivot/data/SEER")) # so use absolute path. Here df holds the full set of data field (df) options
(rdf=pickFields(df)) # this is a reduced data field (rdf) dataframe (note fewer rows AND one extra integer/string column)
# for (i in c("73","92","00")) mkSEERold(rdf,dataset=i,mkDFs=T) #(old way) populates binaries into these folders
mkSEER(rdf,seerHome="/Users/radivot/data/SEER") #(new way) makes merged (all cancer) binaries in SEER data folder mrgd

# these are the default picks for pickFields (see its help page)
defPicks=c("casenum","reg","race","sex","agedx",
        "yrbrth","seqnum","modx","yrdx","histo3","radiatn","recno",
        "agerec","ICD9","COD","surv")
# the required core of this default list is
c("reg","race","sex","agedx","histo3","radiatn","agerec","ICD9")

# the user is responsible for assuring that additional fields are inserted in
# their correct positions corresponding to the SAS file order.
pickFields(df,picks=c(defPicks,"ICD10","siterwho")) # e.g. this won't fly due to negative field widths near the bottom

# instead, the user has to just bite the bullet and do this
picks=c("casenum","reg","race","sex","agedx",
        "yrbrth","seqnum","modx","yrdx","histo3","radiatn","recno",
        "agerec","siterwho","ICD9","ICD10","COD","surv") # add ICD10 and who2008 to list, in correct positions
(rdf=pickFields(df,picks))


## STAGES: if you want to explore different stage fields, you may want something like this
picks=c("casenum","reg","race","sex","agedx","yrbrth",
        "seqnum","yrdx","histo2","histo3","eod10sz","eod10nd","cstumsiz","dajcct","dajccn","dajccm", "radiatn","agerec",
        "siterwho","ICD9","ICD10","histrec","hststga","ajccstg","aj3seer",
        "COD","odthclass","surv",
        "dajcc7t","dajcc7n","dajcc7m","dajcc7stg")
(rdf=pickFields(df,picks))
mkSEER(rdf,outFile="cancStgs",writePops=F) #80 secs + 120 secs


## ALL COLUMNS ######## make one with all columns just to see how big things 
#get. The following defeats one of the main points of this package, which is to 
#gain speed by focusing only on fields of interest. The resulting binaries are 
#thus not likely to be useful on a regular basis, but may be useful from a
#computer programming perspective for checking out speed in the limit of using all fields.
rdf=pickFields(df,picks=df$names)
mkSEER(rdf,outFile="cancAll",writePops=F) #210 secs, i.e. 3.5 minutes. 3-fold more columns=> ~4.5-fold more time
# making the SQL db was an additional 248 secs, i.e. ~4 minutes. 

# If you want to check to see what fields you have in a binary right now, you can do this
system.time(load("~/data/SEER/mrgd/cancDef.RData")) # 2 secs to load 9M cases is fast relative to having all fields
head(canc,2)
system.time(load("~/data/SEER/mrgd/cancALL.RData")) # which almost takes 30 secs
head(canc,2)

