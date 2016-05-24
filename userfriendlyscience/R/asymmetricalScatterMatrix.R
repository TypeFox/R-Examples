asymmetricalScatterMatrix <- function(dat, cols, rows,
#                                      scaleLimits = NULL,
#                                      powerHist=TRUE,
                                      theme=dlvTheme(),
                                      autoSize = TRUE,
                                      txtHeight = 1,
                                      histHeight = 3,
                                      scatterWidth = 6,
                                      scatterHeight = 6,
                                      unit = 'cm',
                                      dpi=200,
                                      showCorrelations=c('top-left',
                                                         'top-right',
                                                         'bottom-left',
                                                         'bottom-right'),
                                      correlationSize = 15,
                                      correlationColor = "#cadded",
                                      pointSize = 1.5) {
  
  ### Generate object with 3 sub-objects to store input,
  ### intermediate results, and output
  res <- list(input = as.list(environment()),
              intermediate = list(),
              output = list());
  
  ### Extract dataframe and select only complete cases
  res$intermediate$dat <-
    dat <-
    na.omit(dat[, c(cols, rows)]);

  ### Convert all variables to numeric vectors, if they weren't already
  res$intermediate$dat <-
    dat <-
    massConvertToNumeric(res$intermediate$dat);
  
  res$intermediate$plots <- list();
  res$intermediate$plotList <- list();
  for (currentRowVar in -1:length(rows)) {
    res$intermediate$plots[[currentRowVar+2]] <- list();
    for (currentColVar in -1:length(cols)) {
      
      if ((currentRowVar < 1) && (currentColVar < 1)) {
        ### Top-left corner, display nothing
        res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]] <-
          grid.rect(gp=gpar(col="white"));
        
      } else if (currentRowVar == -1) {
        ### In the first row: show column variable name
        res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]] <-
          textGrob(cols[currentColVar]);
        
      } else if (currentColVar == -1) {
        ### In the first column: show row variable name
        res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]] <-
          textGrob(rows[currentRowVar], rot=90);
        
      } else if (currentRowVar == 0) {
        ### In the second row: show column variable histogram
        res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]] <-
          powerHist(dat[[cols[[currentColVar]]]], xLabel=FALSE,
                    yLabel=FALSE, distributionLineSize=.5,
                    normalLineSize = .5)$plot +
          theme(axis.title=element_blank(),
                plot.margin=unit(rep(1, 4), "mm"));
        
      } else if (currentColVar == 0) {
        
        ### In the second column: show row variable histogram
        res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]] <-
          powerHist(dat[[rows[[currentRowVar]]]], xLabel=FALSE,
                    yLabel=FALSE, distributionLineSize=.5,
                    normalLineSize = .5)$plot +
          scale_y_reverse() + coord_flip() +
          theme(axis.title=element_blank(),
                plot.margin=unit(rep(1, 4), "mm"));
        
      } else {
        ### We're beyond the second row or column; show scatterplot
        
        ggplotDf <- data.frame(x = dat[, cols[[currentColVar]] ],
                               y = dat[, rows[[currentRowVar]] ]);
#         
#         corLabel <- data.frame(x = min(dat[, cols[[currentColVar]] ]),
#                                y = min(dat[, y[[currentRowVar]] ]),
#                                label = noZero(round(cor(dat[[x[[currentColVar]]]],
#                                                         dat[[y[[currentRowVar]]]]), 3)));
        
        ### Create the plot; note that we start with adding the scatterplot
        ### to get the dimensions, then we add the text, and then the scatterplot
        ### again so it overlays the text.
        
        jitteredPointsLayer <- geom_point(aes(x=x, y=y), position='jitter',
                                          size=pointSize);

        res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]] <-
          #ggplot(dat, aes_string(x=x[[currentColVar]], y=y[[currentRowVar]])) +
          ggplot(data = ggplotDf, aes(x=x, y=y));
        # + jitteredPointsLayer +
        
        if ("bottom-left" %IN% showCorrelations) {
          res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]] <-
            res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]] +
          geom_text(x = min(dat[, cols[[currentColVar]] ]),
                    y = min(dat[, rows[[currentRowVar]] ]),
                    label = noZero(round(cor(dat[[cols[[currentColVar]]]],
                                             dat[[rows[[currentRowVar]]]]), 3)),
                    size=correlationSize, color=correlationColor, vjust=0, hjust=0);
        }
        if ("top-left" %IN% showCorrelations) {
          res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]] <-
            res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]] +
            geom_text(x = min(dat[, cols[[currentColVar]] ]),
                    y = max(dat[, rows[[currentRowVar]] ]),
                    label = noZero(round(cor(dat[[cols[[currentColVar]]]],
                                             dat[[rows[[currentRowVar]]]]), 3)),
                    size=correlationSize, color=correlationColor, vjust=1, hjust=0);
        }
        if ("bottom-right" %IN% showCorrelations) {
          res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]] <-
            res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]] +
            geom_text(x = max(dat[, cols[[currentColVar]] ]),
                    y = min(dat[, rows[[currentRowVar]] ]),
                    label = noZero(round(cor(dat[[cols[[currentColVar]]]],
                                             dat[[rows[[currentRowVar]]]]), 3)),
                    size=correlationSize, color=correlationColor, vjust=0, hjust=1);
        }
        if ("top-right" %IN% showCorrelations) {
          res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]] <-
            res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]] +
          geom_text(x = max(dat[, cols[[currentColVar]] ]),
                    y = max(dat[, rows[[currentRowVar]] ]),
                    label = noZero(round(cor(dat[[cols[[currentColVar]]]],
                                             dat[[rows[[currentRowVar]]]]), 3)),
                    size=correlationSize, color=correlationColor, vjust=1, hjust=1);
        }
        
        res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]] <-
          res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]] +
          jitteredPointsLayer +
          theme_bw() +
          theme(axis.title=element_blank(),
                plot.margin=unit(rep(1, 4), "mm"));
        
#         print(res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]]);

        ### Add text - has to be separate, otherwise it changes the axes
#         res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]] <-
#           res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]] +
#           geom_text(x = min(dat[, x[[currentColVar]] ]),
#                     y = min(dat[, y[[currentRowVar]] ]),
#                     label = noZero(round(cor(dat[[x[[currentColVar]]]],
#                                              dat[[y[[currentRowVar]]]]), 3)),
#                     size=20, color="#cadded", vjust=0, hjust=0);
#           geom_text(data = corLabel,
#                     mapping = aes(x=x, y=y, label=label),
#                     size=20, color="#cadded", vjust=0, hjust=0) +
          
      }
      ### Now reorganise into one long list, from top-left to bottom-right
      res$intermediate$plotList[[length(res$intermediate$plotList) + 1]] <-
        res$intermediate$plots[[currentRowVar+2]][[currentColVar+2]];
    }
  }
  
  ### Create the final plot
  if (autoSize) {
    res$output$scatterMatrix <-
      do.call(arrangeGrob, c(res$intermediate$plotList,
                             list(ncol=length(cols) + 2,
                                  widths=c(txtHeight, histHeight, rep.int(scatterWidth, length(cols))),
                                  heights=c(txtHeight, histHeight, rep.int(scatterHeight, length(rows))))));
  } else {
    res$output$scatterMatrix <-
      do.call(arrangeGrob, c(res$intermediate$plotList,
                             list(ncol=length(cols) + 2,
                                  widths=unit(c(txtHeight, histHeight, rep.int(scatterWidth, length(cols))), unit),
                                  heights=unit(c(txtHeight, histHeight, rep.int(scatterHeight, length(rows))), unit))));
  }
  
  ### Store the size of the plot
  res$output$plotSize <- list(width = txtHeight + histHeight + scatterWidth * length(cols),
                              height = txtHeight + histHeight + scatterHeight * length(rows),
                              unit=unit,
                              res=dpi);
  
  ### Set class and return result
  class(res) <- "asymmetricalScatterMatrix";
  return(res);
  
}

print.asymmetricalScatterMatrix <- function(x, ...) {
  grid.newpage();
  grid.draw(x$output$scatterMatrix, ...);
  invisible();
}
