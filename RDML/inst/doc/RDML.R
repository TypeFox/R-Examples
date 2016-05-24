## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
library(RDML)

## ---- results = "hide"---------------------------------------------------
PATH <- path.package("RDML")
filename <- paste(PATH, "/extdata/", "lc96_bACTXY.rdml", sep ="")
lc96 <- RDML$new(filename = filename)

## ------------------------------------------------------------------------
lc96

## ------------------------------------------------------------------------
fdata <- 
  lc96$
    experiment$`ca1eb225-ecea-4793-9804-87bfbb45f81d`$
    run$`65aeb1ec-b377-4ef6-b03f-92898d47488b`$
    react$`45`$
    data$bACT$
    adp$fpoints #'adp' means amplification data points (qPCR)
head(fdata)

## ---- results = "hide", fig.width = 6, fig.height = 4--------------------
lc96$AsDendrogram()

## ------------------------------------------------------------------------
ref <- lc96$
          experiment$`ca1eb225-ecea-4793-9804-87bfbb45f81d`$
          run$`65aeb1ec-b377-4ef6-b03f-92898d47488b`$
          react$`39`$
          sample$id
sample <- lc96$sample[[ref]]
sample$quantity$value

## ------------------------------------------------------------------------
id1 <- idType$new("id_1")
id2 <- id1
id3 <- id1$clone(deep = TRUE)
id2$id <- "id_2"
id3$id <- "id_3"
cat(sprintf("Original object\t: %s ('id_1' bacame 'id_2')\nSimple copy\t\t: %s\nClone\t\t\t: %s\n",
            id1$id, id2$id, id3$id))

## ------------------------------------------------------------------------
# Create 'real' copy of object
experiment <- lc96$experiment$`ca1eb225-ecea-4793-9804-87bfbb45f81d`$clone(deep = TRUE)
# Try to set 'id' with wrong input type.
# Correct type 'idType' can be seen at error message.
tryCatch(experiment$id <- "exp1",
         error = function(e) print(e))

# Set 'id' with correct input type - 'idType'
experiment$id <- idType$new("exp1")

# Similar operations for 'run'
run <- experiment$run$`65aeb1ec-b377-4ef6-b03f-92898d47488b`$clone(deep = TRUE)
run$id <- idType$new("run1")

# Replace original elements with modified
experiment$run <- list(run)
lc96$experiment <- list(experiment)

## ---- results = "hide", fig.width = 6, fig.height = 4--------------------
lc96$AsDendrogram()

## ------------------------------------------------------------------------
tab <- lc96$AsTable(
  # Custom name pattern 'position~sample~sample.type~target~dye'
  name.pattern = paste(
             react$Position(run$pcrFormat),
             react$sample$id,
             private$.sample[[react$sample$id]]$type$value,
             data$tar$id,
             target[[data$tar$id]]$dyeId$id,
             sep = "~"),
  # Custom column 'quantity' - starting quantity of added sample 
  quantity = sample[[react$sample$id]]$quantity$value
)
# Remove row names for compact printing
rownames(tab) <- NULL
head(tab)

## ---- results = "hide", messages = FALSE, fig.width = 6, fig.height = 4----
library(dplyr)
library(ggplot2)

# Prepare request to get only 'std' type samples
filtered.tab <- filter(tab,
                       sample.type == "std")

fdata <- lc96$GetFData(filtered.tab,
                       # long table format for usage with ggplot2
                       long.table = TRUE)
ggplot(fdata, aes(cyc, fluor)) +
    geom_line(aes(group = fdata.name,
                  color = target))

## ------------------------------------------------------------------------
library(chipPCR)
tab <- lc96$AsTable(
  # Custom name pattern 'position~sample~sample.type~target~run.id'
  name.pattern = paste(
             react$Position(run$pcrFormat),
             react$sample$id,
             private$.sample[[react$sample$id]]$type$value,
             data$tar$id,
             run$id$id, # run id added to names
             sep = "~"))
# Get all fluorescence data
fdata <- lc96$GetFData(tab,
                       # We don't need long table format for CPP()
                       long.table = FALSE)

fdata.cpp <- cbind(cyc = fdata[, 1],
                   apply(fdata[, -1], 2,
                         function(x) CPP(fdata[, 1],
                                         x)$y))

## ---- fig.width = 6, fig.height = 4--------------------------------------
tab$run.id <- "run.cpp"
# Set fluorescence data from previous section
lc96$SetFData(fdata.cpp,
              tab)

# View setted data
fdata <- lc96$GetFData(tab,
                       long.table = TRUE)
ggplot(fdata, aes(cyc, fluor)) +
    geom_line(aes(group = fdata.name,
                  color = target))

## ---- results = "hide", fig.width = 7, fig.height = 6--------------------
# Load another built in RDML file
stepone <- RDML$new(paste0(path.package("RDML"),
                           "/extdata/", "stepone_std.rdml"))
# Merge it with our 'lc96' object
merged <- MergeRDMLs(list(lc96, stepone))
# View structure of new object
merged$AsDendrogram()

## ---- results = "hide"---------------------------------------------------
RDML$set("public", "CalcCq",
         function() {
           library(chipPCR)
           fdata <- self$GetFData(
             self$AsTable())
           fdata <- cbind(cyc = fdata[, 1],
                          apply(fdata[, -1],
                                2,
                                function(x)
                                  # Data preprocessing
                                  CPP(fdata[, 1],
                                      x)$y)
                          )
           
           apply(fdata[, -1], 2,
                 function(x) {
                   tryCatch(
                     # Calculate Cq
                     th.cyc(fdata[, 1], x,
                            auto = TRUE)@.Data[1],
                     error = function(e) NA)
                 })
         }
)

# Create new object with our advanced class
stepone <- RDML$new(paste0(path.package("RDML"),
                           "/extdata/", "stepone_std.rdml"))

## ------------------------------------------------------------------------
stepone$CalcCq()

## ---- fig.width = 6, fig.height = 3, results = "hide"--------------------
### Create simulated data with AmpSim() from chipPCR package
# Cq for data to be generated
Cqs <- c(15, 17, 19, 21)
# PCR si,ulation will be 35 cycles
fdata <- data.frame(cyc = 1:35)
for(Cq in Cqs) {
  fdata <- cbind(fdata,
                 AmpSim(cyc = 1:35, Cq = Cq)[, 2])
}
# Set names for fluorescence curves
colnames(fdata)[2:5] <- c("c1", "c2", "c3", "c4")

# Create minimal description
descr <- data.frame(
  fdata.name = c("c1", "c2", "c3", "c4"),
  exp.id = c("exp1", "exp1", "exp1", "exp1"),
  run.id = c("run1", "run1", "run1", "run1"),
  react.id = c(1, 1, 2, 2),
  sample = c("s1", "s1", "s2", "s2"),
  target = c("gene1", "gene2", "gene1", "gene2"),
  target.dyeId = c("FAM", "ROX", "FAM", "ROX"),
  stringsAsFactors = FALSE
)

# Create empty RDML object
sim <- RDML$new()
# Add fluorescence data
sim$SetFData(fdata, descr)

# Observe object
sim$AsDendrogram()
fdata <- sim$GetFData(sim$AsTable(),
                      long.table = TRUE)
ggplot(fdata, aes(cyc, fluor)) +
  geom_line(aes(group = fdata.name,
                color = target,
                linetype = sample))

