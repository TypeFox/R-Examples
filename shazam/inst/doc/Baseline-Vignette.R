## ---- eval=TRUE, warning=FALSE, message=FALSE----------------------------
# Set example data
library(shazam)
db <- InfluenzaDb

## ---- eval=TRUE, warning=FALSE, results="hide"---------------------------
library(shazam)

# Collapse clonal groups into single sequences
db_clone <- collapseByClone(db, regionDefinition=IMGT_V_NO_CDR3, 
                            expandedDb=TRUE, nproc=1)

## ---- eval=TRUE, warning=FALSE, results="hide"---------------------------
# Count observed mutations and add OBSERVED columns to db_clone
db_obs <- calcDBObservedMutations(db_clone, 
                                  sequenceColumn="CLONAL_SEQUENCE",
                                  regionDefinition=IMGT_V_NO_CDR3, nproc=1)
# Calculate expected mutations and add EXPECTED columns to db_obs
db_exp <- calcDBExpectedMutations(db_obs, 
                                  sequenceColumn="CLONAL_SEQUENCE",
                                  targetingModel=HS5FModel,
                                  regionDefinition=IMGT_V_NO_CDR3, nproc=1)

## ---- eval=TRUE, warning=FALSE, results="hide"---------------------------
# Calculate selection scores using the output from calcDBExpectedMutations
baseline <- calcBaseline(db_exp, testStatistic="focused", 
                         regionDefinition=IMGT_V_NO_CDR3,
                         nproc=1)

# Subset the original data to two time-points and switched isotypes
db_sub <- subset(db, CPRIMER %in% c("IGHA", "IGHG") & 
                     BARCODE %in% c("RL013", "RL018"))
# Calculate selection scores from scratch on subset
baseline_sub <- calcBaseline(db_sub, testStatistic="focused", 
                             regionDefinition=IMGT_V_NO_CDR3, nproc=1)

## ---- eval=TRUE, warning=FALSE, results="hide"---------------------------
# Combine selection scores by time-point
baseline_one <- groupBaseline(baseline, groupBy=c("BARCODE"))

# Combine selection scores by time-point and isotype
baseline_two <- groupBaseline(baseline, groupBy=c("BARCODE", "CPRIMER"))
baseline_sub <- groupBaseline(baseline_sub, groupBy=c("BARCODE", "CPRIMER"))

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Plot mean and confidence interval by time-point
plotBaselineSummary(baseline_one, "BARCODE")

# Plot selection scores by time-point and isotype for only CDR
plotBaselineSummary(baseline_two, "BARCODE", "CPRIMER", subsetRegions="CDR")

# Plot only two time-points and recolor isotypes
group_colors <- c("IGHM"="darkorchid", "IGHD"="firebrick", 
                  "IGHG"="seagreen", "IGHA"="steelblue")
plotBaselineSummary(baseline_sub, "BARCODE", "CPRIMER", groupColors=group_colors)

# Group by CDR/FWR and facet by isotype
plotBaselineSummary(baseline_two, "BARCODE", "CPRIMER", facetBy="group")

## ---- eval=TRUE, warning=FALSE-------------------------------------------
# Plot selection PDFs for a subset of the data
group_colors <- c("RL013"="steelblue", "RL018"="firebrick")
plotBaselineDensity(baseline_sub, "CPRIMER", "BARCODE",
                    groupColors=group_colors, sigmaLimits=c(-3, 3))

