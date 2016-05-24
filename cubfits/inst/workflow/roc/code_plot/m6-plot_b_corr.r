### This script plots correlation of b (log(mu), Delta t).

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

source("00-set_env.r")
source(paste(prefix$code.plot, "u2-plot_b_corr.r", sep = ""))

if(length(case.names) < 4){
  stop("Need 4 cases to match with.")
}

### Load data.
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

### Ordered by "wophi_pm", "wophi_scuo", "wphi_pm", and "wphi_scuo".
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
  # fn.in <- paste(prefix$subset, case.names[i.case], "_PM_scaling.rda", sep = "")
  # if(!file.exists(fn.in)){
  #   cat("File not found: ", fn.in, "\n", sep = "")
  #   next
  # }
  # load(fn.in)

  b.ci[[i.case]] <- b.ci.PM
  b.mean[[i.case]] <- b.PM
  label[[i.case]] <- b.label
}


### Plot logmu.
all.names <- names(b.PM)
id.intercept <- grep("log.mu", all.names)

if(!is.null(b.mean[[3]]) && !is.null(b.mean[[1]])){
  x.pm <- b.mean[[3]][id.intercept]
  y.pm <- b.mean[[1]][id.intercept]
  x.pm.label <- label[[3]]
  x.pm.ci <- b.ci[[3]][id.intercept,]
  y.pm.ci <- b.ci[[1]][id.intercept,]
  xlim <- my.range(c(x.pm))
  ylim <- my.range(c(y.pm))
}
if(!is.null(b.mean[[4]]) && !is.null(b.mean[[2]])){
  x.scuo <- b.mean[[4]][id.intercept]
  y.scuo <- b.mean[[2]][id.intercept]
  x.scuo.label <- label[[4]]
  x.scuo.ci <- b.ci[[4]][id.intercept,]
  y.scuo.ci <- b.ci[[2]][id.intercept,]
  xlim <- my.range(c(x.scuo))
  ylim <- my.range(c(y.scuo))
}
if(!is.null(b.mean[[3]]) && !is.null(b.mean[[1]]) &&
   !is.null(b.mean[[4]]) && !is.null(b.mean[[2]])){
  xlim <- my.range(c(x.pm, x.scuo))
  ylim <- my.range(c(y.pm, y.scuo))
}

if(!is.null(b.mean[[3]]) && !is.null(b.mean[[1]])){
  fn.out <- paste(prefix$plot.match, "corr_logmu_pm.pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    plot.b.corr(x.pm, y.pm, x.pm.label,
                x.ci = x.pm.ci, y.ci = y.pm.ci,
                xlim = xlim, ylim = ylim,
                xlab = "log(mu) with phi", ylab = "log(mu) without phi",
                main = "roc_pm", add.lm = TRUE,
                workflow.name = workflow.name)
    mtext(date(), line = 2.5, cex = 0.4)
  dev.off()
}

if(!is.null(b.mean[[4]]) && !is.null(b.mean[[2]])){
  fn.out <- paste(prefix$plot.match, "corr_logmu_scuo.pdf", sep = "")
  pdf(fn.out, width = 5, height = 5)
    plot.b.corr(x.scuo, y.scuo, x.scuo.label,
                x.ci = x.scuo.ci, y.ci = y.scuo.ci,
                xlim = xlim, ylim = ylim,
                xlab = "log(mu) with phi", ylab = "log(mu) without phi",
                main = "roc_scuo", add.lm = TRUE,
                workflow.name = workflow.name)
    mtext(date(), line = 2.5, cex = 0.4)
  dev.off()
}


### Plot Delta.t.
id.slop <- grep("Delta.t", all.names)

if(!is.null(b.mean[[3]]) && !is.null(b.mean[[1]])){
  x.pm <- b.mean[[3]][id.slop]
  y.pm <- b.mean[[1]][id.slop]
  x.pm.label <- label[[3]]
  x.pm.ci <- b.ci[[3]][id.slop,]
  y.pm.ci <- b.ci[[1]][id.slop,]
  xlim <- my.range(c(x.pm))
  ylim <- my.range(c(y.pm))
}
if(!is.null(b.mean[[4]]) && !is.null(b.mean[[2]])){
  x.scuo <- b.mean[[4]][id.slop]
  y.scuo <- b.mean[[2]][id.slop]
  x.scuo.label <- label[[4]]
  x.scuo.ci <- b.ci[[4]][id.slop,]
  y.scuo.ci <- b.ci[[2]][id.slop,]
  xlim <- my.range(c(x.scuo))
  ylim <- my.range(c(y.scuo))
}
if(!is.null(b.mean[[3]]) && !is.null(b.mean[[1]]) &&
   !is.null(b.mean[[4]]) && !is.null(b.mean[[2]])){
  xlim <- my.range(c(x.pm, x.scuo))
  ylim <- my.range(c(y.pm, y.scuo))
}

if(!is.null(b.mean[[3]]) && !is.null(b.mean[[1]])){
  fn.out <- paste(prefix$plot.match, "corr_deltat_pm.pdf", sep = "")
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
  fn.out <- paste(prefix$plot.match, "corr_deltat_scuo.pdf", sep = "")
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
