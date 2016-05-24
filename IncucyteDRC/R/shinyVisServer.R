#' Shiny Visualisation Server
#'
#' Creates a shiny server for the interactive web app for the running the IncucyteDRC workflow
#'
#' @param input Shiny input list
#' @param output Shiny output list
#' @return A shiny server
#' @export
#' @import shiny
shinyVisServer <- function(input, output) {

    user_pm <- reactive({

        out1 <- try(importPlatemapXML(input$platemap_file$datapath), silent=TRUE)
        out2 <- try(importPlatemap(input$platemap_file$datapath), silent=TRUE)

        if(class(out1) == 'data.frame') {
            return(out1)
        } else  if (class(out2) == 'data.frame') {
            return(out2)
        } else {
            return(NULL)
        }
    })

    user_data <- reactive({
        out <- try(importIncucyteData(input$data_file$datapath, metric='pc', plateid=input$data_file$name), silent=TRUE)
        if(class(out) == 'IncucyteDRCPlateData') {
            return(out)
        } else {
            return(NULL)
        }
    })

    res_list <- reactive({
        if(is.null(user_pm()) | is.null(user_data())){
            return(NULL)
        } else {
            splitIncucyteDRCPlateData(user_pm(), user_data(), group_columns=input$group_columns_select)
        }
    })

    list_element <- reactive({
        if(is.null(input$metadata_row_last_clicked)) {
            1
        } else {
            as.numeric(input$metadata_row_last_clicked)
        }
    })

    res <- reactive({
        idrc_set <- res_list()

        if (class(idrc_set) == 'IncucyteDRCSetList') {
            idrc_set <- idrc_set[[list_element()]]
        }

        idrc_set <- idrc_set %>%
            fitGrowthCurvesGrouped() %>%
            fitGrowthCurvesIndividual()

        if (input$cut_time_mode) {
                idrc_set <- idrc_set %>%
                    calculateCutTimeForIDRCSet(baseline_time=input$baseline_time_slider,
                                               no_doublings=input$no_doublings_slider,
                                               max_val=input$max_val_slider) %>%
                    calculateDRCData()
        } else {
            idrc_set <- calculateDRCData(idrc_set, input$cut_time_slider)
        }
        return(idrc_set)
        })

    res_ec50 <- reactive({
        res() %>% fitDoseResponseCurve(include_control = input$include_control_mode) %>%
            calculateEC50()
    })

    dotmatics_data <- reactive({
        exportDRCDataToDotmatics(res(), user_pm())
    })

    drc_data <- reactive({
        if(input$data_format_select == 'dataframe'){
            exportDRCDataToDataFrame(res(), add_metadata=TRUE, include_control = input$include_control_mode)
        } else if (input$data_format_select == 'prism') {
            exportDRCDataToPRISM(res(), add_metadata=TRUE, include_control = input$include_control_mode)
        } else {
            dotmatics_data()$data
        }
    })

    ec50_data <- reactive({
        exportEC50Data(res_ec50(), add_metadata=FALSE)
    })


    metadata_df <- reactive({
        if(class(res_list()) != 'IncucyteDRCSetList') {
            res()$metadata
        } else {
            lapply(res_list(), function(x) x$metadata) %>% dplyr::bind_rows()
        }
    })

    output$mainpage_ui <- renderUI({
        if(is.null(user_pm()) & is.null(user_data())) {
            mainPanel(
                helpText('Upload a valid platemap and data file to start')
            )
        } else if (!is.null(user_pm()) & is.null(user_data())) {
            mainPanel(
                helpText("Platemap uploaded successfully, now upload a data file"),
                plotOutput('platemap_plot'),
                downloadButton('download_platemap_data', 'Download Platemap')
            )
        } else {
            mainPanel(
                helpText("Click on the table to select a dataset:"),
                DT::dataTableOutput('metadata'),
                plotOutput('plot'))

        }

    })

    output$datapage_ui <- renderUI({
        if(is.null(res_list())) {
            mainPanel(
                helpText('Upload a valid platemap and data file to start')
            )
        } else if (input$data_format_select == 'dotmatics') {
            if(class(res_list()) == 'IncucyteDRCSetList') {
                dotmatics_warn <- span("WARNING: this doesn't work properly when there are is more than one cell growth condition in a plate", style = "color:red")
            } else {
                dotmatics_warn <- ''
            }
            mainPanel(
                helpText(dotmatics_warn),
                downloadButton('download_samplelist_data', 'Download Samplelist'),
                downloadButton('download_drc_data', 'Download Data'),
                tableOutput('samplelist_data'),
                tableOutput('drc_data')
            )
        } else {
            mainPanel(
                tableOutput('drc_data'),
                downloadButton('download_drc_data', 'Download Data')
            )
        }
    })

    output$ec50page_ui <- renderUI({
        if(is.null(res_list())) {
            mainPanel(
                helpText('Upload a valid platemap and data file to start')
            )
        } else {
            mainPanel(
                plotOutput('ec50_plot'),
                DT::dataTableOutput('ec50_data'),
                downloadButton('download_ec50_data', 'Download Data')
            )
        }
    })

    output$platemap_ui <- renderUI({
        if(is.null(user_pm())) {
            mainPanel(
                helpText('Upload a valid platemap file to start')
            )
        } else {
            mainPanel(
                helpText("Platemap uploaded successfully"),
                plotOutput('platemap_plot'),
                downloadButton('download_platemap_data', 'Download Platemap')
            )
        }
    })

    output$cut_time_res_ui <- renderUI({
        if(input$cut_time_mode) {
            mainPanel(
                helpText("Cut time algorithm details"),
                plotOutput('cut_time_plot'),
                tableOutput('cut_time_table')
            )
        } else {
            mainPanel(
                helpText('Cut time calculation not enabled')
                )

        }
    })



    output$plot <- renderPlot({
        plotIncucyteDRCSet(res(), grouped=TRUE)
    })

    output$platemap_plot <- renderPlot({
        plotPlatemap(user_pm())
    })

    output$download_platemap_data <- downloadHandler(
        filename = 'platemap_download.txt',
        content = function(file) {
            write.table(user_pm(), file=file, sep='\t', row.names = FALSE, col.names = TRUE, na='')
        }
    )

    output$metadata <- DT::renderDataTable({
        metadata_df()
    }, filter='none', selection=list(mode = 'single', selected = 1),
        options=list(searching=FALSE, sorting=FALSE, paging=FALSE, info=FALSE, ordering=FALSE))

    output$drc_data <- renderTable({
        drc_data()
    })

    output$samplelist_data <- renderTable({
        dotmatics_data()$samplelist
    })

    output$download_drc_data <- downloadHandler(
        filename = 'drc_data_download.txt',
        content = function(file) {
            write.table(drc_data(), file=file, sep='\t', row.names = FALSE, col.names = TRUE, na='')
        }
    )

    output$download_samplelist_data <- downloadHandler(
        filename = 'samplelist_data_download.txt',
        content = function(file) {
            write.table(dotmatics_data()$samplelist, file=file, sep='\t', row.names = FALSE, col.names = TRUE, na='')
        }
    )

    output$ec50_plot <- renderPlot({

        rowid <- input$ec50_data_row_last_clicked
        rowid <- ifelse(is.null(rowid),1,rowid)
        sid <- ec50_data()[rowid,'sampleid']
        message(sid)
        plotDoseResponseCurve(res_ec50(), sid)
    })

    output$ec50_data <- DT::renderDataTable({
        ec50_data()
    }, filter='none', selection=list(mode = 'single', selected = 1),
    options=list(searching=FALSE, sorting=FALSE, paging=FALSE, info=FALSE, ordering=FALSE))

    output$download_ec50_data <- downloadHandler(
        filename = 'ec50_data_download.txt',
        content = function(file) {
            write.table(exportEC50Data(res_ec50(), add_metadata=TRUE),
                        file=file, sep='\t', row.names = FALSE, col.names = TRUE, na='')
        }
    )

    output$cut_time_table <- renderTable({
        t(res()$calculated_cut)
    })

    output$cut_time_plot <- renderPlot({
        print(res()$cut_plot)
    })

    output$help_ui <- renderUI({
        if(input$help_select == 'none') {
            helpText('Click on a drop down to select help')
        } else if (input$help_select == 'overview') {
            div(includeHTML(system.file('doc/Overview.html', package = 'IncucyteDRC')))
        } else if (input$help_select == 'zoom') {
            div(includeHTML(system.file('doc/IncucyteZoomExport.html', package = 'IncucyteDRC')))
        } else if (input$help_select == 'video') {
            fluidRow(helpText('View the video below to find out more about this web app:'),
                     tags$iframe(width="800", height="600",
                                 src="http://www.youtube.com/embed/EQzeKefmZVw?autoplay=0"))
        }

    })


}
