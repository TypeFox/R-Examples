#
# This file belongs to mldrGUI, an EDA GUI for multilabel datasets developed on top of the mldr package.
# (c) 2015 - Francisco Charte Ojeda (fcharte@ugr.es), David Charte Luque (fdavidcl@outlook.com)
# See the package LICENSE file for license information
#

shinyServer(function(input, output, session) {
  selected <- NULL  # mldr selected by the user in the drop-down list

  observe({
    if(input$loadButton != 0) {
      isolate({
        arfffile <- input$arffname
        xmlfile <- input$xmlname

        if (!is.null(arfffile)) {
          # Load the dataset in the global environment
          if (is.null(xmlfile)) # MEKA
            .GlobalEnv[[arfffile$name]] <-
              mldr(arfffile$datapath, auto_extension = FALSE, use_xml = FALSE)
          else # MULAN
            .GlobalEnv[[arfffile$name]] <-
              mldr(arfffile$datapath, auto_extension = FALSE, xml_file = xmlfile$datapath)

          selected <- arfffile$name
        }
      })
    }

    # Make sample datasets in package mldr available in global environment
    tryCatch(
      for (obj in ls("package:mldr"))
        if (class(get(obj, "package:mldr")) == "mldr")
          assign(obj, get(obj, "package:mldr"), .GlobalEnv),
      error = function(e) {
        warning("Couldn't load sample datasets. Make sure package mldr is attached (library(mldr)).")
      }
    )

    # Get available mldr objects in the global environment
    availableMLDs <- as.list(
      ls(.GlobalEnv)[unlist(sapply(ls(.GlobalEnv),
                                   function(obj) class(get(obj)) == "mldr"))
                     ]
    )

    # Select first available mldr or an empty one
    if (is.null(selected))
      selected <- if (length(availableMLDs) > 0)
        availableMLDs[[1]]
      else
        mldr::mldr()

    updateSelectInput(session, "mldrs",
                      choices = availableMLDs,
                      selected = selected)
  })

  title <- reactive({ paste(input$mldrs, " - mldr", sep = "")})
  output$title <- renderText(title())

  # Table with summary information about the mldr
  summaryTable <- reactive({
    input$loadButton

    if(!is.null(input$mldrs) && input$mldrs != "") {
      mld <- get(input$mldrs)
      titles <- c("Number of attributes",
                  "Number of instances",
                  "Number of labels",
                  "Number of labelsets")
      values <- c(mld$measures$num.attributes,
                  mld$measures$num.instances,
                  mld$measures$num.labels,
                  mld$measures$num.labelsets)
      table <- data.frame(Description = titles, Value = values)
      table
    }
  })
  output$summaryGeneral <- renderTable(summaryTable(), include.rownames = FALSE, digits = 0)

  summaryLabelsTable <- reactive({
    input$loadButton

    if(!is.null(input$mldrs) && input$mldrs != "") {
      mld <- get(input$mldrs)
      titles <- c("Cardinality",
                  "Density",
                  "Most frequent %",
                  "Least frequent %")
      values <- c(mld$measures$cardinality,
                  mld$measures$density,
                  max(mld$labels$freq)*100,
                  min(mld$labels$freq)*100)
      table <- data.frame(Description = titles, Value = values)

      titles <- c("Max imbalance ratio",
                  "Mean imbalance ratio",
                  "Max SCUMBLE",
                  "Mean SCUMBLE")
      values <- c(max(mld$labels$IRLbl),
                  mld$measures$meanIR,
                  max(mld$dataset$.SCUMBLE),
                  mld$measures$scumble)
      table <- cbind(table, data.frame(Description = titles, Value = values))
      table
    }
  })
  output$summaryLabels <- renderTable(summaryLabelsTable(), include.rownames = FALSE, digits = 4)

  summaryLabelsetsTable <- reactive({
    input$loadButton

    if(!is.null(input$mldrs) && input$mldrs != "") {
      mld <- get(input$mldrs)
      titles <- c("Number of single labelsets",
                  "Most frequent %",
                  "Least frequent %")
      values <-
        c(mld$measures$num.single.labelsets,
                  max(mld$labelsets) / mld$measures$num.instances * 100,
                  min(mld$labelsets) / mld$measures$num.instances * 100)

      table <- data.frame(Description = titles, Values = values)
      table
    }
  })
  output$summaryLabelsets <- renderTable(summaryLabelsetsTable(), include.rownames = FALSE, digits = 4)

  # Table with data about the labels in the mldr
  labelsTable <- reactive({
    if(!is.null(input$mldrs) && input$mldrs != "") {
      mld <- get(input$mldrs)
      tbl <-  cbind(Label = rownames(mld$labels), mld$labels)
      tbl
    }
  })
  output$labels <- renderDataTable(labelsTable(), options = list(
    "dom" = 'T<"clear">lfrtip',
    "oTableTools" = list(
      "sSwfPath" = "//cdnjs.cloudflare.com/ajax/libs/datatables-tabletools/2.1.5/swf/copy_csv_xls.swf",
      "aButtons" = list(
        "copy",
        "print",
        list("sExtends" = "collection",
             "sButtonText" = "Save",
             "aButtons" = c("csv","xls")
        )
      )
    )
  ))

  labelsNum <- reactive({
    if(!is.null(input$mldrs) && input$mldrs != "") {
      mld <- get(input$mldrs)
      mld$measures$num.labels
    }
  })
  output$labelRange <- renderUI({
    if(!is.null(input$mldrs) && input$mldrs != "") {
      sliderInput("labelRange", label = h5("Choose range of labels to plot"),
                  min = 1, max = labelsNum(), step = 1,
                  value = c(1, if(labelsNum() < 25) labelsNum() else 25),
                  width = "100%")
    }
  })

  labelHC <- reactive({
    if(!is.null(input$mldrs) && input$mldrs != "" && !is.null(input$labelRange)) {
      mld <- get(input$mldrs)
      labelRange <- input$labelRange
      plot(mld, title = mld$name, type = "LB",
           labelIndices = (labelRange[1] + mld$labels$index[1] - 1):(labelRange[2] + mld$labels$index[1] -1))
    }
  })
  output$labelHC <- renderPlot(labelHC(), height = 800, width = 1024)

  output$saveLabels <- downloadHandler(
    filename = "labels.png",
    content = function(file) {
      mld <- get(input$mldrs)
      labelRange <- input$labelRange
      png(file, type = 'cairo', width = 1024, height = 800)
      plot(mld, title = mld$name, type = "LB",
           labelIndices = (labelRange[1] + mld$labels$index[1] - 1):(labelRange[2] + mld$labels$index[1] -1))
      dev.off()
    },
    contentType = 'image/png'
  )

  attributeByType <- reactive({
    if(!is.null(input$mldrs) && input$mldrs != "") {
      mld <- get(input$mldrs)
      plot(mld, title = mld$name, type = "AT")
    }
  })
  output$attributeByType <- renderPlot(attributeByType(), height = 600, width = 600)

  cardHistogram <- reactive({
    if(!is.null(input$mldrs) && input$mldrs != "") {
      mld <- get(input$mldrs)
      plot(mld, title = mld$name, type = "CH")
    }
  })
  output$cardHistogram <- renderPlot(cardHistogram(), height = 600, width = 600)

  labelHistogram <- reactive({
    if(!is.null(input$mldrs) && input$mldrs != "") {
      mld <- get(input$mldrs)
      plot(mld, title = mld$name, type = "LH")
    }
  })
  output$labelHistogram <- renderPlot(labelHistogram(), height = 600, width = 600)

  labelsetHistogram <- reactive({
    if(!is.null(input$mldrs) && input$mldrs != "") {
      mld <- get(input$mldrs)
      plot(mld, title = mld$name, type = "LSH")
    }
  })
  output$labelsetHistogram <- renderPlot(labelsetHistogram(), height = 600, width = 600)

  labelsetHC <- reactive({
    if(!is.null(input$mldrs) && input$mldrs != "") {
      mld <- get(input$mldrs)
      plot(mld, title = mld$name, type = "LSB")
    }
  })
  output$labelsetHC <- renderPlot(labelsetHC(), height = 600, width = 600)

  # Table with data about labelsets in the mldr
  labelsetsTable <- reactive({
    if(!is.null(input$mldrs) && input$mldrs != "") {
      mld <- get(input$mldrs)
      if(mld$measures$num.instances > 0)
        data.frame(LabelSet = names(mld$labelsets), Count = mld$labelsets)
      else
        data.frame(LabelSet = character(0), Count = numeric(0))
    }
  })
  output$labelsets <- renderDataTable(labelsetsTable(), options = list(
    "dom" = 'T<"clear">lfrtip',
    "oTableTools" = list(
      "sSwfPath" = "//cdnjs.cloudflare.com/ajax/libs/datatables-tabletools/2.1.5/swf/copy_csv_xls.swf",
      "aButtons" = list(
        "copy",
        "print",
        list("sExtends" = "collection",
             "sButtonText" = "Save",
             "aButtons" = c("csv","xls")
        )
      )
    )
  ))

  # Table with data about the attributes in the mldr
  attributesTable <- reactive({
    if(!is.null(input$mldrs) && input$mldrs != "") {
      mld <- get(input$mldrs)
      tbl <- mld$attributes[-mld$labels$index]
      sum <- lapply(names(mld$dataset[-c(mld$labels$index,(length(mld$dataset)-1):length(mld$dataset))]),
                    function(column.name) {
                      tmpsum <- if(mld$attributes[column.name] == 'numeric')
                        summary(mld$dataset[,column.name])
                      else
                        summary(as.factor(mld$dataset[,column.name]))

                      paste('<table class="table table-striped"><tr><td><b>',
                            paste(names(tmpsum), collapse = '</b></td><td><b>'),
                            '</td></tr><tr><td>',
                            paste(tmpsum, collapse = '</td><td>'),
                            '</td></tr></table>')
                    })

      tbl <- data.frame(Attribute = names(tbl), Type = tbl, Summary = unlist(sum))
      tbl
    }
  })
  output$attributes <- renderDataTable(attributesTable(), options = list(
    "dom" = 'T<"clear">lfrtip',
    "oTableTools" = list(
      "sSwfPath" = "//cdnjs.cloudflare.com/ajax/libs/datatables-tabletools/2.1.5/swf/copy_csv_xls.swf",
      "aButtons" = list(
        "copy",
        "print",
        list("sExtends" = "collection",
             "sButtonText" = "Save",
             "aButtons" = c("csv","xls")
        )
      )
    )
  ), escape = FALSE)

  # Table with data about the labels in the mldr
  concurrenceTable <- reactive({
    if(!is.null(input$mldrs) && input$mldrs != "") {
      mld <- get(input$mldrs)
      tbl <- data.frame(Index = mld$labels$index,
                        Label = rownames(mld$labels),
                        Count = mld$labels$count,
                        SCUMBLE = mld$labels$SCUMBLE)
      tbl <- tbl[order(tbl$Count, tbl$SCUMBLE),]

      tbl
    }
  })

  concurrenceAnalysisTable <- reactive({
    if(!is.null(input$mldrs) && input$mldrs != "") {
      mld <- get(input$mldrs)
      lblint <- labelInteractions(mld)
      titles <- c("Dataset",
                  "Mean SCUMBLE",
                  "SCUMBLE CV",
                  "Minority labels with high SCUMBLE",
                  names(lblint$interactions)
                )
      values <- c(mld$name,
                  mld$measures$scumble,
                  mld$measures$scumble.cv,
                  paste(rownames(mld$labels[mld$labels$index %in% lblint$indexes,]), collapse = ", "),
                  sapply(lblint$interactions, function(x) paste(rownames(mld$labels[mld$labels$index %in% row.names(x), ]), collapse = ", "))
                )
      table <- data.frame(Description = titles, Value = values)
      table
    }
  })
  output$ConcurrenceAnalysis <- renderTable(concurrenceAnalysisTable(), include.rownames = FALSE)

  highScumbleLabels <- reactive({
    if(!is.null(input$mldrs) && input$mldrs != "") {
      mld <- get(input$mldrs)
       tbl <- data.frame(Label = rownames(mld$labels),
                         Count = mld$labels$count,
                         SCUMBLE = mld$labels$SCUMBLE)
       tbl <- tbl[order(tbl$Count, tbl$SCUMBLE),]

       ScumbleList <- sort(tbl$SCUMBLE, decreasing = TRUE)
       ScumbleList <- ScumbleList[1:10]
       ScumbleList <- ScumbleList[!is.na(ScumbleList)]
       ScumbleList <- paste(which(tbl$SCUMBLE %in% ScumbleList) - 1, collapse = ",")
#       lblint <- labelInteractions(mld)
#       ScumbleList <- paste(as.numeric(c(lblint$indexes, unique(unlist(lapply(lblint$interactions, names))))) - mld$labels$index[1], collapse = ",")
      paste("function(settings, json) {
            $('.dataTable').DataTable().rows([", ScumbleList, "]).nodes().to$().addClass('selected');
            Shiny.onInputChange('labels', [", ScumbleList, "]);
            }")
    }
  })
  tblConcurrenceOptions <- reactive({
    list(paging = FALSE, searching = FALSE,
         ordering = FALSE, info = FALSE,
         initComplete = I(highScumbleLabels()))
  })

  output$tblConcurrence <- renderDataTable(
    concurrenceTable()[,-1],
    options = tblConcurrenceOptions,
    callback = "function(table) {
      table.on('click.dt', 'tr', function() {
        $(this).toggleClass('selected');
        Shiny.onInputChange('labels',
          table.rows('.selected').indexes().toArray());
      });
    }")

  output$saveConcurrence <- downloadHandler(
    filename = "concurrence.png",
    content = function(file) {
      mld <- get(input$mldrs)
      labels <- concurrenceTable()[input$labels+1,1]
      png(file, type = 'cairo', width = 1024, height = 1024)
      plot(mld, title = mld$name, labelIndices = labels)
      dev.off()
    },
    contentType = 'image/png'
  )

  output$saveAT <- downloadHandler(
    filename = "attributeByType.png",
    content = function(file) {
      mld <- get(input$mldrs)
      png(file, type = 'cairo', width = 1024, height = 1024)
      plot(mld, title = mld$name, type = "AT")
      dev.off()
    },
    contentType = 'image/png'
  )

  output$saveCH <- downloadHandler(
    filename = "cardHistogram.png",
    content = function(file) {
      mld <- get(input$mldrs)
      png(file, type = 'cairo', width = 1024, height = 1024)
      plot(mld, title = mld$name, type = "CH")
      dev.off()
    },
    contentType = 'image/png'
  )

  output$saveLH <- downloadHandler(
    filename = "labelHistogram.png",
    content = function(file) {
      mld <- get(input$mldrs)
      png(file, type = 'cairo', width = 1024, height = 1024)
      plot(mld, title = mld$name, type = "LH")
      dev.off()
    },
    contentType = 'image/png'
  )

  output$saveLSH <- downloadHandler(
    filename = "labelsetHistogram.png",
    content = function(file) {
      mld <- get(input$mldrs)
      png(file, type = 'cairo', width = 1024, height = 1024)
      plot(mld, title = mld$name, type = "LSH")
      dev.off()
    },
    contentType = 'image/png'
  )

  output$saveLabelsets <- downloadHandler(
    filename = "labelsetBarPlot.png",
    content = function(file) {
      mld <- get(input$mldrs)
      png(file, type = 'cairo', width = 1024, height = 1024)
      plot(mld, title = mld$name, type = "LSB")
      dev.off()
    },
    contentType = 'image/png'
  )

  observe({
    input$labels
    labelLC <- reactive({
      if(!is.null(input$mldrs) && input$mldrs != "" && !is.null(input$labels)) {
        mld <- get(input$mldrs)
        labels <- concurrenceTable()[input$labels+1,1]
        plot(mld, title = mld$name, labelIndices = labels)
      }
    })

    output$labelLC <- renderPlot(labelLC(), height = 800, width = 800)
  })


  observe({
    if (is.null(input$pages) || input$pages != "finish")
      return()

    stopApp(0)
  })

})
