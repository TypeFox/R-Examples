#' Class \code{"count_test"} 
#' 
#' Objects of this class are created by \code{\link{test_counts}}.
#' 
#' @name count_test
#' @aliases count_test-class count_test
#' @rdname count_test-class
#' @docType class
#' @section Slots: 
#' \describe{ 
#' \item{group_coef}{\code{"data.frame"} containing experiments, groups to which they
#' belong and calculated values of rate (lambda).}
#' \item{test_res}{\code{"matrix"} containing result of multiple comparisions t-test.} 
#' \item{model}{\code{"character"} name of GLM used to compare experiments.} }
#' @author Michal Burdukiewicz.
#' @seealso \code{\link{test_counts}}.
#' @export
#' @keywords classes
setClass("count_test", representation(group_coef = "data.frame", 
                                      test_res = "matrix",
                                      model = "character"))

#' @describeIn count_test Summary statistics of assigned groups.
#' @export
setMethod("summary", signature(object = "count_test"), function(object) {
  aggregate(. ~ group, slot(object, "group_coef"), mean)
})

#' @describeIn count_test Extract coefficients of groups.
#' @param object of class \code{count_test}.
#' @export
setMethod("coef", signature(object = "count_test"), function(object) {
  slot(object, "group_coef")
})


#' @describeIn count_test Print both \code{group_coef} and \code{test_res}.
#' @export
setMethod("show", "count_test", 
          function(object) {
            cat("Groups:\n")
            print(slot(object, "group_coef"))
            
            cat("\nResults of multiple comparison:\n")
            
            signif_stars <- symnum(slot(object, "test_res")[, "p_value"], corr = FALSE, na = FALSE, 
                                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                    symbols = c("***", "**", "*", ".", " "))
            print(data.frame(slot(object, "test_res"), signif = as.vector(signif_stars)))

            cat("---\nSignif. codes:  ", attr(signif_stars, "legend"), sep = "", 
                fill = getOption("width") + 4 + max(nchar(attr(signif_stars, "legend"), "bytes") 
                                                    - nchar(attr(signif_stars, "legend"))))

            cat("\nTest: ", slot(object, "model"))
          })




#' @describeIn count_test plots mean number of molecules per partition and its 
#' confidence intervals.
#' @aliases plot.count_test plot,count_test-method plot,count_test,ANY-method
#' @param x object of class \code{count_test}.
#' @param aggregate logical, if \code{TRUE} experiments are aggregated according
#' to their group.
#' @param nice logical, if \code{TRUE} a more aesthetically pleasing (but harder to 
#' customize) version of the plot is created.
#' @details In case of the aggregated plot, mean confidence intervals for groups are 
#' presented as dashed lines.
#' @export
setMethod("plot", signature(x = "count_test"), function(x, aggregate = FALSE, 
                                                        nice = TRUE) {
  group_coef <- slot(x, "group_coef")
  if (aggregate) {
    summ <- aggregate(. ~ group, group_coef, mean)
    #possible groups
    pos_groups <- group_coef[["group"]]
    
    plot(c(0.55, nrow(summ) + 0.45), range(summ[, c("lambda.low", "lambda.up")]), 
         xlab = "Group", ylab = expression(lambda), xaxt = "n", cex = 0)
    
    #colors and axis setup
    colors <- if(nice) {
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
           col = adjustcolor("grey", alpha.f = 0.30))
      axis(1, tck = 1, col.ticks = "white", labels = FALSE, at = 1L:nrow(summ))
      axis(2, tck = 1, col.ticks = "white", labels = FALSE)
      box()
      rainbow(nlevels(pos_groups))
    } else {
      rep("black", nlevels(pos_groups))
    }
    
    sapply(1L:nrow(summ), function(i)
      axis(side = 1, labels = summ[i, "group"], at = i, col.axis = colors[i]))
    
    sapply(1L:length(levels(pos_groups)), function(i) {
      points(i + seq(-2, 2, length.out = sum(levels(pos_groups)[i] == pos_groups))/10,
             group_coef[group_coef[["group"]] == levels(pos_groups)[i], "lambda"],
             col = colors[i], pch = ifelse(nice, 16, 1))
    })
    
    abline(h = summ[1L:nrow(summ), "lambda.low"], lty = "dashed", col = colors)
    abline(h = summ[1L:nrow(summ), "lambda.up"], lty = "dashed", col = colors)
    
    
  } else {
    #positions of experiments
    exp_pos <- 1L:nrow(group_coef)
                       
    plot(exp_pos, group_coef[["lambda"]], 
         ylim = range(group_coef[, c("lambda.low", "lambda.up")]), xaxt = "n",
         xlab = "Experiment", ylab = expression(lambda), cex = 0)
    
    #colors and axis setup
    colors <- if(nice) {
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
           col = adjustcolor("grey", alpha.f = 0.30))
      axis(1, tck = 1, col.ticks = "white", labels = FALSE, at = exp_pos)
      axis(2, tck = 1, col.ticks = "white", labels = FALSE)
      box()
      rainbow(nlevels(group_coef[["group"]]))
    } else {
      rep("black", nlevels(group_coef[["group"]]))
    } 
    
    
    sapply(1L:nlevels(group_coef[["group"]]), function(i) {
      lab_pos <- exp_pos[group_coef[["group"]] == levels(group_coef[["group"]])[i]]
      points(lab_pos, group_coef[lab_pos, "lambda"], col = colors[i],
             pch = ifelse(nice, 16, 1))
      #errorbars
      sapply(lab_pos, function(j)
        lines(c(j, j), group_coef[j, c("lambda.low", "lambda.up")],
              col = colors[i], lwd = 1.5))
      #top axis with group names
      axis(side = 3, labels = rep(levels(group_coef[["group"]])[i], length(lab_pos)), 
           at = lab_pos, col.axis = colors[i])
      })
    
    axis(1, at = exp_pos, labels = rownames(group_coef))
  }
})


