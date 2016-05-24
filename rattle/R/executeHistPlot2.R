#' Perform the required operations for displaying histograms using ggplot2.
#' 
#' Time-stamp: <2015-09-17 21:08:44 gjw>
#' 
executeHistPlot2 <- function(dataset, vars, target, targets, stratify, sampling, pmax)
{
  startLog(Rtxt("Display histogram plots for the selected variables."))

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

    plot.cmd <- stringr::str_c('# Generate the plot.\n\n',
                               sprintf("p%02d", i), ' <- crs %>%\n',
                               '  with(', dataset, ') %>%\n',
                               '  dplyr::select(', vars[i],
                               ifelse(length(target), stringr::str_c(", ", target), ""),
                               ') %>%\n',
                               '  ggplot2::ggplot(ggplot2::aes(x=', vars[i], ')) +\n',
                               '  ggplot2::geom_density(lty=3) +\n',
                               ifelse(length(target),
                                      stringr::str_c('  ggplot2::geom_density(ggplot2',
                                                     sprintf("::aes(fill=%s, colour=%s)",
                                                             target, target),
                                                     ', alpha=0.55) +\n'),
                                      ""),
                               '  ggplot2::xlab("', vars[i],
                               '\\n\\n', title.txt[2], '") +\n',
                               '  ggplot2::ggtitle("', title.txt[1], '") +\n',
                               '  ggplot2::labs(',
                               ifelse(length(target),
                                      stringr::str_c('fill="', target, '", '),
                                      ""),
                               'y="Density")')

    ## plot.cmd <- stringr::str_c('# Calculate the variable value range.\n\n',
    ##                            'vrange <- crs %>%\n',
    ##                            '  with(', dataset, ') %>%\n', # Need access to crs vars
    ##                            '  dplyr::select(', vars[i], ') %>%\n',
    ##                            '  range(na.rm=TRUE)\n\n',
    ##                            '# Then detemine a good bin width for the bars.\n\n',
    ##                            'bwidth <- crs %>%\n',
    ##                            '  with(', dataset, '$', vars[i], ') %>%\n',
    ##                            '  na.omit() %>%\n',
    ##                            '  nclass.FD() %>%\n',
    ##                            '  magrittr::divide_by(vrange[2]-vrange[1], .)\n\n',
    ##                            '# Generate the plot.\n\n',
    ##                            sprintf("p%02d", i), ' <- crs %>%\n',
    ##                            '  with(', dataset, ') %>%\n',
    ##                            '  dplyr::select(', vars[i],
    ##                            ifelse(length(target), stringr::str_c(", ", target), ""),
    ##                            ') %>%\n',
    ##                            '  ggplot2::ggplot(ggplot2::aes(x=', vars[i], ')) +\n',
    ##                            '  ggplot2::geom_histogram(ggplot2::aes(y=..density..), ',
    ##                            'binwidth=bwidth, fill="grey", colour="black") +\n',
    ##                            '  ggplot2::geom_density(', 
    ##                            ifelse(length(target),
    ##                                   sprintf("ggplot2::aes(colour=%s)", target), ""),
    ##                            ') +\n',
    ##                            '  ggplot2::xlab("', vars[i],
    ##                            '\\n\\n', title.txt[2], '") +\n',
    ##                            '  ggplot2::ggtitle("', title.txt[1], '") +\n',
    ##                            '  ggplot2::labs(colour="", y="Density")')

    comment <- paste(Rtxt("Use ggplot2 to generate histogram plot for"), vars[i])
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
