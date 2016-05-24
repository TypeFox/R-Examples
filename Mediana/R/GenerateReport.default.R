######################################################################################################################

# Function: GenerateReport.
# Argument: ResultDes returned by the CSE function and presentation model and Word-document title and Word-template.
# Description: This function is used to create a summary table with all results
#' @export
#' @import ReporteRs

GenerateReport.default = function(presentation.model = NULL, cse.results, report.filename, report.template = NULL){
  # Add error checks
  if (!is.null(presentation.model) & class(presentation.model) != "PresentationModel") stop("GenerateReport: the presentation.model parameter must be a PresentationModel object.")
  if (class(cse.results) != "CSE") stop("GenerateReport: the cse.results parameter must be a CSE object.")
  if (!is.character(report.filename)) stop("GenerateReport: the report.filename parameter must be character.")

  # Create the structure of the report
  # If no presentation model, initialize a presentation model object
  if (is.null(presentation.model)) presentation.model = PresentationModel()
  report = CreateReportStructure(cse.results, presentation.model)
  report.results = report$result.structure
  report.structure = report$report.structure

  # Delete an older version of the report
  if (!is.null(report.filename)){
    if (file.exists(report.filename)) file.remove(report.filename)
  }

  # Create a DOCX object
  if (!missing(report.template)) {
    doc = docx(title = report.structure$title, template = report.template)
  } else {
    # Use standard template
    doc = docx(title = report.structure$title)
  }

  # Report's title
  doc = addParagraph(doc, value = report.structure$title, stylename = "TitleDoc")

  # Text formatting
  my.text.format = parProperties(text.align = "left")

  # Table formatting
  header.cellProperties = cellProperties(border.left.width = 0, border.right.width = 0, border.bottom.width = 2, border.top.width = 2, padding = 5, background.color = "#eeeeee")
  data.cellProperties = cellProperties(border.left.width = 0, border.right.width = 0, border.bottom.width = 1, border.top.width = 0, padding = 3)

  header.textProperties = textProperties(font.size = 10, font.weight = "bold")
  data.textProperties = textProperties(font.size = 10)

  leftPar = parProperties(text.align = "left")
  rightPar = parProperties(text.align = "right")
  centerPar = parProperties(text.align = "center")

  # Number of sections in the report (the report's title is not counted)
  n.sections = length(report.structure$section)

  # Loop over the sections in the report
  for(section.index in 1:n.sections) {

    # Section's title (if non-empty)
    if (!is.na(report.structure$section[[section.index]]$title)) doc = addTitle(doc, value = report.structure$section[[section.index]]$title, 1)

    # Number of subsections in the current section
    n.subsections = length(report.structure$section[[section.index]]$subsection)

    # Loop over the subsections in the current section
    for(subsection.index in 1:n.subsections) {

      # Subsection's title (if non-empty)
      if (!is.na(report.structure$section[[section.index]]$subsection[[subsection.index]]$title)) doc = addTitle(doc, value = report.structure$section[[section.index]]$subsection[[subsection.index]]$title, 2)

      # Number of subsubsections in the current section
      n.subsubsections = length(report.structure$section[[section.index]]$subsection[[subsection.index]]$subsubsection)

      if (n.subsubsections>0){
        # Loop over the subsubsection in the current section
        for(subsubsection.index in 1:n.subsubsections) {

          # Subsubsection's title (if non-empty)
          if (!is.na(report.structure$section[[section.index]]$subsection[[subsection.index]]$subsubsection[[subsubsection.index]]$title)) doc = addTitle(doc, value = report.structure$section[[section.index]]$subsection[[subsection.index]]$subsubsection[[subsubsection.index]]$title, 3)

          # Number of subsubsubsections in the current section
          n.subsubsubsection = length(report.structure$section[[section.index]]$subsection[[subsection.index]]$subsubsection[[subsubsection.index]]$subsubsubsection)

          if (n.subsubsubsection>0){
            # Loop over the subsubsubsection in the current section
            for(subsubsubsection.index in 1:n.subsubsubsection) {

              # Subsubsubsection's title (if non-empty)
              if (!is.na(report.structure$section[[section.index]]$subsection[[subsection.index]]$subsubsection[[subsubsection.index]]$subsubsubsection[[subsubsubsection.index]]$title)) doc = addTitle(doc, value = report.structure$section[[section.index]]$subsection[[subsection.index]]$subsubsection[[subsubsection.index]]$subsubsubsection[[subsubsubsection.index]]$title, 4)

              # Number of items in the current subsubsection
              n.items = length(report.structure$section[[section.index]]$subsection[[subsection.index]]$subsubsection[[subsubsection.index]]$subsubsubsection[[subsubsubsection.index]]$item)

              # Loop over the items in the current subsection
              for(item.index in 1:n.items) {

                # Create paragraphs for each item

                # Determine the item's type (text by default)
                type = report.structure$section[[section.index]]$subsection[[subsection.index]]$subsubsection[[subsubsection.index]]$subsubsubsection[[subsubsubsection.index]]$item[[item.index]]$type
                if (is.null(type)) type = "text"

                label = report.structure$section[[section.index]]$subsection[[subsection.index]]$subsubsection[[subsubsection.index]]$subsubsubsection[[subsubsubsection.index]]$item[[item.index]]$label
                value = report.structure$section[[section.index]]$subsection[[subsection.index]]$subsubsection[[subsubsection.index]]$subsubsubsection[[subsubsubsection.index]]$item[[item.index]]$value
                param = report.structure$section[[section.index]]$subsection[[subsection.index]]$subsubsection[[subsubsection.index]]$subsubsubsection[[subsubsubsection.index]]$item[[item.index]]$param

                if (type == "table" & is.null(param)) param = list(span.columns = NULL, groupedheader.row = NULL)

                switch( type,
                        text = {
                          if (label != "") doc = addParagraph(doc, value = paste(label, value), stylename = "Normal", par.properties = my.text.format)
                          else  doc = addParagraph(doc, value = value, stylename = "Normal", par.properties = my.text.format)
                        },
                        table = {
                          header.columns = (is.null(param$groupedheader.row))
                          summary_table = FlexTable(data = value, body.cell.props = data.cellProperties, header.cell.props = header.cellProperties, header.columns = header.columns )
                          if (!is.null(param$span.columns)) {
                            for (ind.span in 1:length(param$span.columns)){
                              summary_table = spanFlexTableRows(summary_table, j = param$span.columns[ind.span], runs = as.character(value[,ind.span]) )
                            }
                          }
                          summary_table = setFlexTableBorders(summary_table, inner.vertical = borderNone(),
                                                              outer.vertical = borderNone())
                          if (!is.null(param$groupedheader.row)) {
                            summary_table = addHeaderRow(summary_table, value = param$groupedheader.row$values, colspan = param$groupedheader.row$colspan)
                            summary_table = addHeaderRow(summary_table, value = colnames( value ))
                          }
                          doc = addParagraph(doc, value = label, stylename = "rTableLegend", par.properties = my.text.format)
                          doc = addFlexTable(doc, summary_table)
                        },
                        plot =  {
                          doc = addPlot(doc, fun = print, x = value, width = 6, height = 5, main = label)
                          doc = addParagraph(doc, value = label, stylename = "rPlotLegend", par.properties = my.text.format)
                        }
                )
              }
            }
          }
          else {
            # Number of items in the current subsubsection
            n.items = length(report.structure$section[[section.index]]$subsection[[subsection.index]]$subsubsection[[subsubsection.index]]$item)

            # Loop over the items in the current subsection
            for(item.index in 1:n.items) {

              # Create paragraphs for each item

              # Determine the item's type (text by default)
              type = report.structure$section[[section.index]]$subsection[[subsection.index]]$subsubsection[[subsubsection.index]]$item[[item.index]]$type
              if (is.null(type)) type = "text"

              label = report.structure$section[[section.index]]$subsection[[subsection.index]]$subsubsection[[subsubsection.index]]$item[[item.index]]$label
              value = report.structure$section[[section.index]]$subsection[[subsection.index]]$subsubsection[[subsubsection.index]]$item[[item.index]]$value
              param = report.structure$section[[section.index]]$subsection[[subsection.index]]$subsubsection[[subsubsection.index]]$item[[item.index]]$param

              if (type == "table" & is.null(param)) param = list(span.columns = NULL, groupedheader.row = NULL)

              switch( type,
                      text = {
                        if (label != "") doc = addParagraph(doc, value = paste(label, value), stylename = "Normal", par.properties = my.text.format)
                        else  doc = addParagraph(doc, value = value, stylename = "Normal", par.properties = my.text.format)
                      },
                      table = {
                        header.columns = (is.null(param$groupedheader.row))
                        summary_table = FlexTable(data = value, body.cell.props = data.cellProperties, header.cell.props = header.cellProperties, header.columns = header.columns )
                        if (!is.null(param$span.columns)) {
                          for (ind.span in 1:length(param$span.columns)){
                            summary_table = spanFlexTableRows(summary_table, j = param$span.columns[ind.span], runs = as.character(value[,ind.span]) )
                          }
                        }
                        summary_table = setFlexTableBorders(summary_table, inner.vertical = borderNone(),
                                                            outer.vertical = borderNone())
                        if (!is.null(param$groupedheader.row)) {
                          summary_table = addHeaderRow(summary_table, value = param$groupedheader.row$values, colspan = param$groupedheader.row$colspan)
                          summary_table = addHeaderRow(summary_table, value = colnames( value ))
                        }
                        doc = addParagraph(doc, value = label, stylename = "rTableLegend", par.properties = my.text.format)
                        doc = addFlexTable(doc, summary_table)
                      },
                      plot =  {
                        doc = addPlot(doc, fun = print, x = value, width = 6, height = 5, main = label)
                        doc = addParagraph(doc, value = label, stylename = "rPlotLegend", par.properties = my.text.format)
                      }
              )
            }
          }
        }
      }
      else {

        # Number of items in the current subsection
        n.items = length(report.structure$section[[section.index]]$subsection[[subsection.index]]$item)

        # Loop over the items in the current subsection
        for(item.index in 1:n.items) {

          # Create paragraphs for each item

          # Determine the item's type (text by default)
          type = report.structure$section[[section.index]]$subsection[[subsection.index]]$item[[item.index]]$type
          if (is.null(type)) type = "text"

          label = report.structure$section[[section.index]]$subsection[[subsection.index]]$item[[item.index]]$label
          value = report.structure$section[[section.index]]$subsection[[subsection.index]]$item[[item.index]]$value
          param = report.structure$section[[section.index]]$subsection[[subsection.index]]$item[[item.index]]$param

          if (type == "table" & is.null(param)) param = list(span.columns = NULL, groupedheader.row = NULL)

          switch( type,
                  text = {
                    if (label != "") doc = addParagraph(doc, value = paste(label, value), stylename = "Normal", par.properties = my.text.format)
                    else  doc = addParagraph(doc, value = value, stylename = "Normal", par.properties = my.text.format)
                  },
                  table = {
                    header.columns = (is.null(param$groupedheader.row))
                    summary_table = FlexTable(data = value, body.cell.props = data.cellProperties, header.cell.props = header.cellProperties, header.columns = header.columns )
                    if (!is.null(param$span.columns)) {
                      for (ind.span in 1:length(param$span.columns)){
                        summary_table = spanFlexTableRows(summary_table, j = param$span.columns[ind.span], runs = as.character(value[,ind.span]) )
                      }
                    }
                    summary_table = setFlexTableBorders(summary_table, inner.vertical = borderNone(),
                                                        outer.vertical = borderNone())
                    if (!is.null(param$groupedheader.row)) {
                      summary_table = addHeaderRow(summary_table, value = param$groupedheader.row$values, colspan = param$groupedheader.row$colspan)
                      summary_table = addHeaderRow(summary_table, value = colnames( value ))
                    }
                    doc = addParagraph(doc, value = label, stylename = "rTableLegend", par.properties = my.text.format)
                    doc = addFlexTable(doc, summary_table)
                  },
                  plot =  {
                    doc = addPlot(doc, fun = print, x = value, width = 6, height = 5, main = label)
                    doc = addParagraph(doc, value = label, stylename = "rPlotLegend", par.properties = my.text.format)
                  }
          )
        }

      }

    }
  }

  # Save the report
  writeDoc(doc, report.filename)

  # Return
  return(invisible(report.results))

}
# End of GenerateReport
