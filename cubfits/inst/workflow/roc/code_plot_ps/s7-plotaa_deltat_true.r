### This is for simulation only, ploting delta.t * phi aganst true values.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

source("00-set_env.r")
source(paste(prefix$code.plot.ps, "u0-get_case_main.r", sep = ""))
source(paste(prefix$code, "u1-get_negsel.r", sep = ""))
source(paste(prefix$code.plot.ps, "u4-plot_aa_allinone.r", sep = ""))

### Load true.
fn.in <- paste(prefix$data, "simu_true_", model, ".rda", sep = "")
if(file.exists(fn.in)){
  load(fn.in)
} else{
  stop(paste(fn.in, " is not found.", sep = ""))
}

b.Init <- convert.b.to.bVec(Eb)

### Load initial.
fn.in <- paste(prefix$data, "/pre_process.rda", sep = "")
load(fn.in)

### Get AA and synonymous codons.
aa.names <- names(reu13.df.obs)
coef.names <- cubfits:::get.my.coefnames(model)
label <- NULL
b.names <- NULL
for(i.codon in aa.names){
  tmp <- sort(unique(reu13.df.obs[[i.codon]]$Codon))
  tmp <- tmp[-length(tmp)]
  label <- c(label, paste(i.codon, tmp, sep = "."))
  b.names <- c(b.names, rep(coef.names, each = length(tmp)))
}

### Get id.slop.
all.names <- b.names
id.slop <- grep("Delta.t", all.names)

### Convert true to negsel and delta.t only.
tmp <- get.negsel(b.Init, id.slop, aa.names, label)
b.Init <- tmp$b.negsel.PM
label.negsel.true <- tmp$b.negsel.label


### Plot by case.
for(i.case in case.names){
  ### Load subset mcmc run.
  fn.in <- paste(prefix$subset, i.case, ".rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  ### Convert unscaled result to negsel and delta.t only.
  tmp <- lapply(1:ncol(b.mcmc),
           function(i.iter){
             tmp <- get.negsel(b.mcmc[, i.iter], id.slop, aa.names, label)
             tmp$b.negsel.PM
           })
  b.mcmc <- do.call("cbind", tmp)

  ### Load scaled posterior mean (negsel).
  fn.in <- paste(prefix$subset, i.case, "_PM_scaling.rda", sep = "")
  load(fn.in)

  ### Plot.
  t.phi.mcmc <- t(phi.mcmc)
  fn.out <- paste(prefix$plot.ps.AA, i.case, "_deltat.phi.pdf", sep = "")
  pdf(fn.out, width = 12, height = 11)
    nf <- layout(matrix(c(rep(1, 5), 2:21), nrow = 5, ncol = 5, byrow = TRUE),
                 rep(1, 5), c(2, 8, 8, 8, 8), respect = FALSE)
    ### Plot title.
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
    text(0.5, 0.6,
         paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""))
    text(0.5, 0.4, date(), cex = 0.6)
    par(mar = c(5.1, 4.1, 4.1, 2.1))

    for(i.aa in aa.names){
      id.label <- grep(paste("^", i.aa, "\\.", sep = ""), label)
      tl.codon <- length(id.label)

      plot.aa.allinone(i.aa, id.label, tl.codon,      ### For AA.
                       label.negsel.true,             ### For label.
                       b.Init, EPhi,                   ### For true.
                       b.mcmc, t.phi.mcmc,            ### For unscaled results.
                       b.negsel.PM, phi.PM,           ### For scaled results.
                       workflow.name, i.case, model)
    }

    plot(NULL, NULL, type = "n", axes = FALSE,
         xlim = c(0, 1), ylim = c(0, 1),
         xlab = "", ylab = "", main = "")
    legend(0, 0.9, "1-to-1", col = 4, lty = 2)
  dev.off()
}

