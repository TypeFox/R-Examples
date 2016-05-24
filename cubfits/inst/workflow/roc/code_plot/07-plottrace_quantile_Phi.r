rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

### Load environment and set data.
source("00-set_env.r")
source(paste(prefix$code.plot, "u0-get_case_main.r", sep = ""))
source(paste(prefix$code.plot, "u3-plot_trace.r", sep = ""))
fn.in <- paste(prefix$data, "init_", model, ".rda", sep = "")
load(fn.in)
fn.in <- paste(prefix$data, "simu_true_", model, ".rda", sep = "")
if(file.exists(fn.in)){
  load(fn.in)
}

q.probs <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)

### Trace each run.
for(i.case in case.names){
  ### All mcmc outputs.
  fn.in <- paste(prefix$output, i.case, "/output_mcmc.rda", sep = "")
  if(!file.exists(fn.in)){
    cat("File not found: ", fn.in, "\n", sep = "")
    next
  }
  load(fn.in)

  ### Since my.appr() doesn't have phi.Mat, but have phi.pred.Mat
  if(is.null(ret[["phi.Mat"]])){
    ret$phi.Mat <- ret$phi.pred.Mat
  }

  ### Find genes by quantile.
  ret.phi.Mat <- ret$phi.Mat
  # ret.phi.Mat.scaled <- lapply(ret.phi.Mat, function(x) x / mean(x))
  ret.phi.Mat.scaled <- ret.phi.Mat

  ret.phi.Mat.PM <- rowMeans(do.call("cbind", ret.phi.Mat.scaled))
  q.PM <- quantile(ret.phi.Mat.PM, probs = q.probs)

  id.gene <- NULL
  for(i.q in 1:length(q.PM)){
    id.gene <- c(id.gene, which.min(abs(ret.phi.Mat.PM - q.PM[i.q])))
  }

  ### Subset.
  phi.sub <- lapply(ret.phi.Mat, function(x){ x[id.gene] })
  phi.PM.sub <- lapply(ret.phi.Mat.scaled, function(x){ x[id.gene] })

  ### Get trace plot.
  trace.Mat <- log10(do.call("cbind", phi.sub))
  xlim.trace <- c(1, ncol(trace.Mat))
  ylim.trace <- range(trace.Mat)

  ### Get hist plot.
  hist.Mat <- log10(do.call("cbind", phi.PM.sub)[, range$subset])
  hist.list <- apply(hist.Mat, 1,
                     function(x){ hist(x, nclass = 50, plot = FALSE) })
  hist.mean <- rowMeans(hist.Mat)
  xlim.hist <- range(hist.Mat)
  ylim.hist <- range(unlist(lapply(hist.list,
                                   function(p){ range(p$counts) })))

  ### Get non-log hist plot.
  hist.Mat.nl <- do.call("cbind", phi.PM.sub)[, range$subset]
  hist.mean.nl <- rowMeans(hist.Mat.nl)
  hist.list.nl <- apply(hist.Mat.nl, 1,
                        function(x){ hist(x, nclass = 50, plot = FALSE) })
  xlim.hist.nl <- range(hist.Mat.nl)
  ylim.hist.nl <- range(unlist(lapply(hist.list.nl,
                                      function(p){ range(p$counts) })))

  ### Set layout.
  fn.out <- paste(prefix$plot.trace, "quantile_Phi_", i.case, ".pdf", sep = "")
  pdf(fn.out, width = 9, height = 10)

  ### Plot trace for each quantile.
  for(i.q in 1:length(q.PM)){
### New page.
    if(i.q %% 3 == 1){
      nf <- layout(matrix(c(1, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                          nrow = 4, ncol = 3, byrow = TRUE),
                   c(1, 1, 1), c(2, 8, 8, 8), respect = FALSE)
      ### Plot title.
      par(mar = c(0, 0, 0, 0))
      plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
      text(0.5, 0.6,
           paste(workflow.name, ", ", get.case.main(i.case, model), sep = ""))
      text(0.5, 0.4, date(), cex = 0.6)
      par(mar = c(5.1, 4.1, 4.1, 2.1))
    }

    ### Plot trace
    plot(1:ncol(trace.Mat), trace.Mat[i.q,],
         type = "l",
         xlim = xlim.trace, ylim = ylim.trace,
         xlab = "Iterations", ylab = "Production Rate (log10)",
         main = paste(names(id.gene)[i.q], ", q = ", q.probs[i.q], sep = ""))
    # abline(h = hist.mean[i.q], col = 2)
    # if(exists("EPhi")){
    #   abline(h = log10(EPhi[id.gene[i.q]]), col = 4, lty = 2)
    # }

    ### Plot hist
    plot(hist.list[[i.q]],
         # xlim = xlim.hist, ylim = ylim.hist,
         xlab = "Posterior Production Rate (log10)",
         main = paste(names(id.gene)[i.q], ", sdlog = ",
                      sprintf("%.4f", sd(hist.Mat[i.q,])), sep = ""))
    abline(v = hist.mean[i.q], col = 2)
    if(exists("EPhi")){
      abline(v = log10(EPhi[id.gene[i.q]]), col = 4, lty = 2)
    }

    ### Plot non-log hist
    plot(hist.list.nl[[i.q]],
         # xlim = xlim.hist.nl, ylim = ylim.hist.nl,
         xlab = "Posterior Production Rate",
         main = paste(names(id.gene)[i.q], ", sd = ",
                      sprintf("%.4f", sd(hist.Mat.nl[i.q,])), sep = ""))
    abline(v = hist.mean.nl[i.q], col = 2)
    if(exists("EPhi")){
      abline(v = EPhi[id.gene[i.q]], col = 4, lty = 2)
    }
  }

  dev.off()
}

