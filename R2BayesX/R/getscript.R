getscript <-
function(object, file = NULL, device = NULL, ...)
{
  if(!inherits(object, "bayesx"))
    stop("object must be of class \'bayesx\'")
  if(is.null(file))
    dir <- getwd()
  else 
    dir <- file
  dode <- !is.null(device)
  if(dode) {
    dcall <- as.character(match.call()$device)
    args <- list(...)
    na <- names(args)
    if(!is.null(na)) {
      na <- paste(na, "=", args, collapse = ", ")
      na <- paste(na, ", ", sep = "")
    }
    dcall = paste(dcall, "(", na, sep = "")
  }
  on <- deparse(substitute(object), backtick = TRUE, width.cutoff = 500L)
  mn <- names(object)
  object <- get.model(object, NULL)
  n <- length(object)
  script <- file("script", "w")
  cat("## The following R code may be applied on object ", on, ".\n", sep = "", file = script)
  cat("## Also see the help files of function plot.bayesx()\n## for all available plots!\n\n", 
    file = script)
  if(n > 1L)
    cat("## Note: object contains of", n, "models:\n\n", file = script)
  cat("## Save the model in working directory:\n", file = script)
  cat("save(", on, ", file = \"", dir, "/", on, ".rda\")\n\n", sep = "", file = script) 
  cat("## Note: model is loaded again with:\n", file = script)
  cat("## load(file = \"", dir,"/", on, ".rda\")\n\n", sep = "", file = script) 
  cat("## Summary statistics are provided with:\n", file = script)
  cat("summary(", on, ")\n", sep = "", file = script) 
  for(i in 1L:length(object)) {
    if(!is.null(object[[i]]$model.fit$method) && (object[[i]]$model.fit$method == "MCMC" ||
      object[[i]]$model.fit$method == "HMCMC")) {
      cat("\n", file = script)
      if(n > 1L)
        cat("## MCMC model diagnostic plots of model", mn[i], "\n", file = script)
      else
        cat("## MCMC model diagnostic plots:\n", file = script)
      if(dode) {
        cat(dcall, "file = \"", paste(dir, "/model-", i, "-max-acf\"", sep = ""), 
          ")\n", sep = "", file = script)
      }
      cat("plot(", on, ", model = ", i, ", which = \"max-acf\")\n", sep = "", file = script)
      if(dode)
        cat("graphics.off()\n", file = script)
    }
    cat("\n", file = script)
    if(!is.null(object[[i]]$effects)) {
      if(n > 1L) {
        cat("## Plots of estimated effects of model: ", mn[i], "\n", sep = "", file = script) 
      } else {
        cat("## Plots of estimated effects\n", sep = "", file = script)
      } 
      ccheck <- vcheck <- FALSE
      neff <- length(object[[i]]$effects)
      for(k in 1L:neff) {
        specs <- attr(object[[i]]$effects[[k]], "specs")
        xlab <- specs$term
        ylab <- specs$label
        if(dode) {
          cat(dcall, "width = 5, height = 4, file = \"", paste(dir, "/model-", i, "-term-", specs$label,"\"", sep = ""), 
            ")\n", sep = "", file = script)
          cat("par(mar = c(4.1, 4.1, 0.1, 1.1))\n", file = script)
        }
        if(inherits(object[[i]]$effects[[k]], c("sm.bayesx", "linear.bayesx")) && specs$dim < 2L) {
          cat("plot(", on, ", model = ", i, ", term = ", "\"", specs$label, 
            "\", xlab = \"", xlab, "\", ylab = \"", ylab,"\")\n", sep = "", file = script)
        } else {
          cat("plot(", on, ", model = ", i, ", term = ", "\"", specs$label, "\")\n", 
            sep = "", file = script)
        }
        if(dode)
          cat("graphics.off()\n", file = script)

        ## now for maps
        if(!is.null(map.name <- attr(object[[i]]$effects[[k]], "map.name"))) {
          if(dode) {
            cat(dcall, "width = 5, height = 4, file = \"", paste(dir, "/model-", i, "-term-", specs$label,"\"", sep = ""), 
              ")\n", sep = "", file = script)
            cat("par(mar = c(0, 0, 0, 0))\n", file = script)
          }
          cat("plot(", on, ", model = ", i, ", term = ", "\"", specs$label, "\", map = ",
            map.name,")\n", sep = "", file = script)
          if(dode)
            cat("graphics.off()\n", file = script)
        }
        if(!is.null(attr(object[[i]]$effects[[k]], "sample")))
          ccheck <- TRUE
        if(!is.null(attr(object[[i]]$effects[[k]], "variance.sample")))
          vcheck <- TRUE
      }
      if(ccheck) {
        cat("\n## Coefficient sampling diagnostics\n", file = script)
        if(!is.null(attr(object[[i]]$fixed.effects, "sample"))) {
          if(dode) {
            cat(dcall, "file = \"", paste(dir, "/model-", i, "-intcpt-samples\"", sep = ""), ")\n", 
              sep = "", file = script)
          }
          cat("plot(", on, ", model = ", i, ", which = \"intcpt-samples\")\n", sep = "", 
            file = script)
          if(dode)
            cat("graphics.off()\n", file = script)
          if(dode) {
            cat(dcall, "file = \"", paste(dir, "/model-", i, "-intcpt-acf\"", sep = ""), ")\n", 
              sep = "", file = script)
          }
          cat("plot(", on, ", model = ", i, ", which = \"intcpt-samples\", acf = TRUE)\n", sep = "", 
            file = script)
          if(dode)
            cat("graphics.off()\n", file = script)
        }
        for(k in 1L:neff) {
          specs <- attr(object[[i]]$effects[[k]], "specs")
          if(!is.null(attr(object[[i]]$effects[[k]], "sample"))) {
            if(dode) {
              cat(dcall, "file = \"", paste(dir, "/model-", i, "-term-", 
                specs$label,"-coef.samples\"", sep = ""), ")\n", sep = "", file = script)
            }
            cat("plot(", on, ", model = ", i, ", term = ", "\"", specs$label, 
              "\", which = \"coef-samples\")\n", sep = "", file = script)
            if(dode)
              cat("graphics.off()\n", file = script)
            if(dode) {
              cat(dcall, "file = \"", paste(dir, "/model-", i, "-term-", 
                specs$label,"-coef.acf\"", sep = ""), ")\n", sep = "", file = script)
            }
            cat("plot(", on, ", model = ", i, ", term = ", "\"", specs$label, 
              "\", which = \"coef-samples\", acf = TRUE)\n", sep = "", file = script)
            if(dode)
              cat("graphics.off()\n", file = script)
            if(dode) {
              cat(dcall, "file = \"", paste(dir, "/model-", i, "-term-", 
                specs$label,"-coef.max.acf\"", sep = ""), ")\n", sep = "", file = script)
            }
            cat("plot(", on, ", model = ", i, ", term = ", "\"", specs$label, 
              "\", which = \"coef-samples\", max.acf = TRUE)\n", sep = "", file = script)
            if(dode)
              cat("graphics.off()\n", file = script)
          }
        }
      }
      if(vcheck) {
        cat("\n## Variance sampling diagnostics\n", file = script)
        for(k in 1L:neff) {
          specs <- attr(object[[i]]$effects[[k]], "specs")
          if(!is.null(attr(object[[i]]$effects[[k]], "variance.sample"))) {
            if(dode) {
              cat(dcall, "file = \"", paste(dir, "/model-", i, "-term-", 
                specs$label,"-variance.samples\"", sep = ""), ")\n", sep = "", file = script)
            }
            cat("plot(", on, ", model = ", i, ", term = ", "\"", specs$label, 
              "\", which = \"var-samples\")\n", sep = "", file = script)
            if(dode)
              cat("graphics.off()\n", file = script)
            if(dode) {
              cat(dcall, "file = \"", paste(dir, "/model-", i, "-term-", 
                specs$label,"-variance.acf\"", sep = ""), ")\n", sep = "", file = script)
            }
            cat("plot(", on, ", model = ", i, ", term = ", "\"", specs$label, 
              "\", which = \"var-samples\", acf = TRUE)\n", sep = "", file = script)
            if(dode)
              cat("graphics.off()\n", file = script)
          }
        }
      }
    }
  }
  close(script)
  script <- rval <- readLines("script")
  unlink("script")
  if(!is.null(file)) {
    file <- path.expand(file)
    cat("", file = file)    
    writeLines(script, con = file)
  }
  class(rval) <- "bayesx.script"

  return(rval)
}


print.bayesx.script <- function(x, ...)
{
  writeLines(x)
  return(invisible(NULL))
}

