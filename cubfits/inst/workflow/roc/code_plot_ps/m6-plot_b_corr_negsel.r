### This script plot correlation of Delta t with adjusted reference codon
### such that all selection coefficients are negative.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

source("00-set_env.r")
source(paste(prefix$code.plot.ps, "u2-plot_b_corr.r", sep = ""))

if(length(case.names) < 4){
  stop("Need 4 cases to match with.")
}

### Load data.
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

### Ordered by "wophi_pm", "wophi_scuo", "wphi_pm", and "wphi_scuo".
b.ci.org <- list(NULL, NULL, NULL, NULL)
b.mean.org <- list(NULL, NULL, NULL, NULL)
label.org <- list(NULL, NULL, NULL, NULL)
b.ci <- list(NULL, NULL, NULL, NULL)
b.mean <- list(NULL, NULL, NULL, NULL)
label <- list(NULL, NULL, NULL, NULL)
for(i.case in 1:4){
  ### Subset of mcmc output.
  fn.in <- paste(prefix$subset, case.names[i.case], "_PM.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  ### Subset of mcmc output with scaling.
  fn.in <- paste(prefix$subset, case.names[i.case], "_PM_scaling.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  b.ci[[i.case]] <- b.negsel.ci.PM
  b.mean[[i.case]] <- b.negsel.PM
  label[[i.case]] <- b.negsel.label
}


### Plot Delta.t.
if(!is.null(b.mean[[3]]) && !is.null(b.mean[[1]])){
  x.pm <- b.mean[[3]]
  y.pm <- b.mean[[1]]
  x.pm.label <- label[[3]]
  x.pm.ci <- b.ci[[3]]
  y.pm.ci <- b.ci[[1]]
  xlim <- my.range(c(x.pm))
  ylim <- my.range(c(y.pm))
}
if(!is.null(b.mean[[4]]) && !is.null(b.mean[[2]])){
  x.scuo <- b.mean[[4]]
  y.scuo <- b.mean[[2]]
  x.scuo.label <- label[[4]]
  x.scuo.ci <- b.ci[[4]]
  y.scuo.ci <- b.ci[[2]]
  xlim <- my.range(c(x.scuo))
  ylim <- my.range(c(y.scuo))
}
if(!is.null(b.mean[[3]]) && !is.null(b.mean[[1]]) &&
   !is.null(b.mean[[4]]) && !is.null(b.mean[[2]])){
  xlim <- my.range(c(x.pm, x.scuo))
  ylim <- my.range(c(y.pm, y.scuo))
}

### Convert selection to negative value by changing the relative base. 
### Note that this can fail in some inconsistent cases.
if(!is.null(b.mean[[3]]) && !is.null(b.mean[[1]])){
  if(any((x.pm < 0 & y.pm > 0) | (x.pm > 0 & y.pm < 0) |
          label[[1]] != label[[3]])){
    stop("Inconsistent cases (PM).")
  }
}
if(!is.null(b.mean[[4]]) && !is.null(b.mean[[2]])){
  if(any((x.scuo < 0 & y.scuo > 0) | (x.scuo > 0 & y.scuo < 0) |
          label[[2]] != label[[4]])){
    stop("Inconsistent cases (SCUO).")
  }
}

### Plot Delta.t.
if(!is.null(b.mean[[3]]) && !is.null(b.mean[[1]])){
  fn.out <- paste(prefix$plot.ps.match, "corr_negsel_deltat_pm.pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    plot.b.corr(x.pm, y.pm, x.pm.label,
                x.ci = x.pm.ci, y.ci = y.pm.ci,
                xlim = xlim, ylim = ylim,
                xlab = "Delta.t with phi", ylab = "Delta.t without phi",
                main = "roc_pm", add.lm = TRUE,
                workflow.name = workflow.name)
    mtext(date(), line = 2.5, cex = 0.4)
  dev.off()
}

if(!is.null(b.mean[[4]]) && !is.null(b.mean[[2]])){
  fn.out <- paste(prefix$plot.ps.match, "corr_negsel_deltat_scuo.pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    plot.b.corr(x.scuo, y.scuo, x.scuo.label,
                x.ci = x.scuo.ci, y.ci = y.scuo.ci,
                xlim = xlim, ylim = ylim,
                xlab = "Delta.t with phi", ylab = "Delta.t without phi",
                main = "roc_scuo", add.lm = TRUE,
                workflow.name = workflow.name)
    mtext(date(), line = 2.5, cex = 0.4)
  dev.off()
}
