### This is a basic configuration for workflow "without phi".

### Specify the case name.
workflow.name <- "wophi"

### Specify model.
model <- "roc"
file.data.fasta <- "genome.fasta"
file.data.tsv <- ""

### Specify output file name.
fn.out.prefix <- "roc_"

### For default path.
prefix <- list()
prefix$root <- "./"
prefix$param <- paste(prefix$root, "param/", sep = "")

prefix$all.out <- "./all.out/"
prefix$data <- paste(prefix$all.out, "data/", sep = "")
prefix$subset <- paste(prefix$all.out, "subset/", sep = "")
prefix$output <- paste(prefix$all.out, "output/", sep = "")

### For table.
prefix$table <- paste(prefix$all.out, "table/", sep = "")
prefix$table.ps <- paste(prefix$table, "table_ps/", sep = "")

### For ploting.
prefix$plot <- paste(prefix$all.out, "plot/", sep = "")
prefix$plot.diag <- paste(prefix$all.out, "plot/diag/", sep = "")
prefix$plot.match <- paste(prefix$all.out, "plot/match/", sep = "")
prefix$plot.single <- paste(prefix$all.out, "plot/single/", sep = "")
prefix$plot.trace <- paste(prefix$all.out, "plot/trace/", sep = "")
prefix$plot.AA <- paste(prefix$all.out, "plot/AA/", sep = "")

### For ploting without scaling.
prefix$plot.ps <- paste(prefix$plot, "plot_ps/", sep = "")
prefix$plot.ps.diag <- paste(prefix$plot, "plot_ps/diag/", sep = "")
prefix$plot.ps.match <- paste(prefix$plot, "plot_ps/match/", sep = "")
prefix$plot.ps.single <- paste(prefix$plot, "plot_ps/single/", sep = "")
prefix$plot.ps.trace <- paste(prefix$plot, "plot_ps/trace/", sep = "")
prefix$plot.ps.AA <- paste(prefix$plot, "plot_ps/AA/", sep = "")

### For code.
prefix$code <- paste(cubfits::get.workflow(model = model), "/",
                     "code/", sep = "")
prefix$code.plot <- paste(cubfits::get.workflow(model = model), "/",
                          "code_plot/", sep = "")
prefix$code.plot.ps <- paste(cubfits::get.workflow(model = model), "/",
                             "code_plot_ps/", sep = "")

### Specify data files.
file.data <- list()
file.data$fasta <- paste(prefix$param, file.data.fasta, sep = "")
if(exists("file.data.tsv") && file.data.tsv != ""){
  file.data$tsv <- paste(prefix$param, file.data.tsv, sep = "")
} else{
  file.data$tsv <- ""
}

### All case names.
case.names <- c("wophi_scuo")
case.names <- paste(model, "_", case.names, sep = "")

### Basic information.
run.info <- list()
run.info$nIter <- 5000

### For configuration.
run.info$dump <- FALSE
run.info$prefix.dump <- paste(prefix$output, "tmp/dump_", sep = "")
# run.info$parallel <- "lapply"
run.info$parallel <- "task.pull"

### For binning plots
run.info$bin.class <- NULL
# run.info$bin.class <- c(0, seq(0.05, 0.95, length = 20), 1)
# run.info$bin.class <- c(seq(0, 0.1, by = 0.02),
#                         seq(0.1, 0.9, by = 0.1),
#                         seq(0.9, 1, by = 0.02))

### For MCMC.
range <- list()
range$subset <- 3001:5000
range$thinning <- 10

### For simulation only.
simulation <- list()
simulation$EPhi <- FALSE
simulation$Eb <- FALSE
simulation$seed <- 1234
simulation$sdlog <- 1.5

### For plotting.
ci.prob <- c(0.025, 0.975)

### For modeling of logmixture.
p.nclass <- 2


### CAUTION: for extra changes globally.
# suppressMessages(library(cubfits, quietly = TRUE))
# .CF.CT$init.fit <- "current"
# .CF.CT$init.fit <- "RW_Norm"
# .CF.CT$type.p <- "lognormal_fix"
# .CF.CT$type.p <- "lognormal_RW"
# .CF.CT$type.p <- "lognormal_bias"
# .CF.CT$model.Phi <- "logmixture"
# .CF.CONF$scale.phi.Obs <- FALSE
# .CF.CONF$estimate.bias.Phi <- TRUE
