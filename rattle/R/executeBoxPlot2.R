#' Perform the required operations for displaying boxplots using ggplot2.
#' 
#' Time-stamp: <2015-09-17 21:08:32 gjw>
#' 
executeBoxPlot2 <- function(dataset, vars, target, targets, stratify, sampling, pmax)
{
  startLog(Rtxt("Display box plots for the selected variables."))

  # We start a new plot since we could be drawing multiple types of
  # plots.
  
  newPlot()

  for (i in seq_along(vars))
  {
    title.txt <- genPlotTitleCmd(generateTitleText(vars[i],
                                                   target,
                                                   sampling,
                                                   stratify && length(targets)),
                                 vector=TRUE)

    plot.cmd <- stringr::str_c('# Generate a box plot.\n\n',
                               sprintf("p%02d", i), ' <- crs %>%\n',
                               '  with(', dataset, ') %>%\n',
                               '  ggplot2::ggplot(ggplot2::aes(y=', vars[i], ')) +\n',
                               '  ggplot2::geom_boxplot(ggplot2::aes(x="All"), ',
                               'notch=TRUE, fill="grey") +\n',
                               '  ggplot2::stat_summary(ggplot2::aes(x="All"), ',
                               'fun.y=mean, geom="point", shape=8) +\n',
                               if (length(target))
                                 stringr::str_c('  ggplot2::geom_boxplot(',
                                                'ggplot2::aes(x=', target, ', ',
                                                'fill=', target, '), notch=TRUE) +\n',
                                                '  ggplot2::stat_summary(',
                                                'ggplot2::aes(x=', target, '), ',
                                                'fun.y=mean, geom="point", ',
                                                'shape=8) +\n'),
                               '  ggplot2::xlab("',
                               if (length(target))
                                 stringr::str_c(target, '\\n\\n'),
                               title.txt[2], '") +\n',
                               '  ggplot2::ggtitle("', title.txt[1], '") +\n',
                               '  ggplot2::theme(legend.position="none")')
  
    comment <- paste(Rtxt("Use ggplot2 to generate box plot for"), vars[i])
    appendLibLog(comment, plot.cmd, include.libs=(i==1))
    eval(parse(text=plot.cmd))
  }

  display.cmd <-
    "gridExtra::grid.arrange(" %s+%
    paste(sprintf("p%02d", seq_len(i)), collapse=", ") %s+%
    ")"

  appendLibLog("Display the plots.", display.cmd)
  eval(parse(text=display.cmd))

}

