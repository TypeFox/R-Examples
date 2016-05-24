

#' Launch Plotter
#' @description Launches the plotting app
#' @param list_of_datasets List of datasets which should be made available for selection when the app is launched
#' @return Launches App
#' @import dplyr shinydashboard ggplot2 shiny
#' @examples
#' \dontrun{
#' diamonds_sample <- ggplot2::diamonds[sample(1:nrow(diamonds),size=1000),]
#' launch_plotter(list(diamonds_sample = diamonds_sample, mtcars = mtcars, iris = iris))}
#' @export

launch_plotter <- function(list_of_datasets) {

#   devtools::use_data(aesthetics=aesthetics,basic_operations=basic_operations,dplyr_tutorial=dplyr_tutorial,loops_tutorial=loops_tutorial,
#                      model_tutorial=model_tutorial, readme_file=readme_file, readme_file_plotter=readme_file_plotter,themes=themes,all_plots=all_plots, col_palette=col_palette, overwrite=T,internal=T)
#

  #### Display reccomended plot lists #####

  Type <- Aes <- Widget <- Palette <- NULL
  g <- ggplot(aesthetics) + geom_bar(aes(Type)) + ggthemes::theme_hc()

  shinyApp(ui <- dashboardPage(
    dashboardHeader(title = "Plotter"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Read Me", tabName = "readme", icon = icon("dashboard")),
        menuItem("Build a Plot", tabName = "datasets", icon = icon("dashboard")),
        menuItem("Plot Zoom View", tabName = "plot_zoom", icon = icon("dashboard"))

      )
    ),
    dashboardBody(
      shinyjs::useShinyjs(),

      tabItems(
        # Dataset Tab
        tabItem(tabName = "datasets",

                fluidRow(
                  tabBox(title = "", width = 6,
                         tabPanel("Select data",
                                  uiOutput("choose_dataset")
                         ),


                         tabPanel("Variables & Geom type",
                                  sliderInput("add_layers","Add additional layer",min=1,max=3,step=1,value=1),

                                  uiOutput("select_var"),
                                  uiOutput("select_plot"),
                                  uiOutput("show_all_plots")


                         ),
                         tabPanel("Set Aesthetics",

                                  uiOutput("show_aesthetics_control"),
                                  uiOutput("map_aesthetic_to_var"),
                                  uiOutput("manual_aesthetic_control"),
                                  uiOutput("reset_buttons")

                         ),

                         tabPanel("Other Plot Controls",

                              uiOutput("collapsible_panel"),
                              shinyBS::bsButton("reset_collapse_panel","Reset All", style = "primary")
                         )

                         ),

                  box(title = "Plot Output", status = "primary", solidHeader = T, width = 6,

                      conditionalPanel(
                        condition = "input.switch_interactivity == false",
                        plotOutput("plot_output")
                      ),
                      conditionalPanel(
                        condition = "input.switch_interactivity == true",
                        plotly::plotlyOutput("interactive_plot_output")
                      ),
                      checkboxInput("switch_interactivity","Switch ON Interactivity", value = F),
                      h4(strong("Underlying Code"),style="color:brown"),
                      textOutput("plot_code")
                  )
                )

        ),
        tabItem(tabName = "readme",
                fluidRow(
                  box(title = "Read Me", status = "primary", solidHeader = T, width = 12,
                      htmlOutput("read_me")))),
        tabItem(tabName = "plot_zoom",
                fluidRow(
                  box(title = "Plot (Enlarged View)", status = "primary", solidHeader = T, width = 12, height = 700,
                      plotOutput("plot_out_zoom"))))
      ))),


    server <- function(input, output, session) {


      #### Read Me ####################
      output$read_me <- reactive({


        val <- as.character(readme_file_plotter[1,1])
        gen_markdown(val)

      })



      #### Tab 1 - Choose a dataset from the dropdown ####################

      output$choose_dataset <- renderUI({

        selectInput("get_dataset","Choose a dataset", choices = names(list_of_datasets))

      })


      #### Display variables from the selected dataset ################

      output$select_var <- renderUI({

        generateVarDropdown(input$get_dataset,input$add_layers)

      })

      ######### Generate Plot dropdown #################

      output$show_all_plots <- renderUI({

        layerCount <- as.integer(input$add_layers)

        lapply(1:layerCount,function(i) {
          checkboxInput(inputId = paste("show_all_plots",i,sep=""), "Show all plots", value = F)
        })
      })



      ####### Show recommended plots basis selected variables ##########

      output$select_plot <- renderUI({

        validate(
          need(try(input[[paste("xaxis_var",1)]] != ""), "")
        )

        selected_dataset <- eval(parse(text = input$get_dataset))
        layerCount <- as.integer(input$add_layers)

        lapply(1:layerCount,function(i) {
          generatePlotList(selected_dataset,input[[paste("xaxis_var",i)]],input[[paste("yaxis_var",i)]],
                           input[[paste("show_all_plots",i,sep="")]],i)
        })
      })


      #### Toggle button - Hide/Show aesthetic controls ############

      output$show_aesthetics_control <- renderUI({

        layerCount <- as.integer(input$add_layers)

        lapply(1:layerCount,function(i) {
          checkboxInput(inputId = paste("show_aesthetic_control",i,sep=""), paste("Show Aesthetic Control for"," Layer",i,sep=""), value = F)
        })
      })

      #### Dynamically generate aesthetic controls (Aesthetic mapped to variable) ######

      output$map_aesthetic_to_var <- renderUI({

        selected_dataset <- eval(parse(text = input$get_dataset))
        vars <- colnames(selected_dataset)

        layerCount <- as.integer(input$add_layers)


        lapply(1:layerCount,function(i) {

          mapAes2Var(input[[paste("plot_type",i)]],vars,i)
        })
      })

      #### Dynamically generate Manual aesthetic Controls ####

      output$manual_aesthetic_control <- renderUI({

        layerCount <- as.integer(input$add_layers)

        lapply(1:layerCount,function(i) {

          applicable_aes <- aesthetics %>% filter(Geom == input[[paste("plot_type",i)]]) %>% select(Aes)
          applicable_aes <- as.character(unlist(applicable_aes))

          aes_widget_form <- aesthetics %>% filter(Geom == input[[paste("plot_type",i)]]) %>% select(Widget, Type)
          length_aes <- nrow(aes_widget_form) - nrow(aes_widget_form[aes_widget_form$Type == "Controls",])
          aes_widget_form <- as.character(unlist(aes_widget_form %>% select(Widget)))


          # reset_value <- input[[paste("reset_button",i,sep="")]]
          id_val <- paste("manual_control_group",i,sep="")
          div(id = id_val,
          lapply(1:length(applicable_aes),function(j) {

            if (j == 1) {
              div(
                h4(strong(paste("LAYER ",i, " - MANUALLY SET AESTHETIC"), style="color:darkblue")),
                eval(parse(text = aes_widget_form[j]))
              )

            } else if (j == (length_aes+1)) {

              div(
              h4(strong(paste("LAYER ",i, " - PLOT SPECIFIC CONTROLS"), style="color:darkblue")),
              eval(parse(text = aes_widget_form[j]))
              )

            } else {
              eval(parse(text = aes_widget_form[j]))
            }


          })
          )

        })
      })

      ############ Generate reset buttons #################

      output$reset_buttons <- renderUI({

        layerCount <- as.integer(input$add_layers)

        lapply(1:layerCount,function(i) {
          shinyBS::bsButton(inputId = paste("reset_button",i,sep=""), "Reset All", style = "primary")
        })
      })

      ############## Reset all values when button is clicked ######

      observe({

        layerCount <- as.integer(input$add_layers)

        lapply(1:layerCount,function(i) {
          observeEvent(input[[paste("reset_button",i,sep="")]], {
            shinyjs::reset(paste("manual_control_group",i,sep=""))
            shinyjs::reset(paste("map_aes_controls",i,sep=""))

          })
        })
      })

      ############ Reset color picker ####################

            observe({

              layerCount <- as.integer(input$add_layers)

              lapply(1:layerCount,function(i) {
                observeEvent(input[[paste("reset_button",i,sep="")]], {
                  applicable_aes <- aesthetics %>% filter(Geom == input[[paste("plot_type",i)]]) %>% select(Aes)
                  applicable_aes <- as.character(unlist(applicable_aes))
                  lapply(1:length(applicable_aes), function(j) {


                    shinyjs::updateColourInput(session,inputId = paste("manual_adjustment",i,j,sep=""), value = "#FFFFFF")
                  })

                })

              })

            })

      ############ Show/Hide Aesthetic Controls ############

      observe({

          layerCount <- as.integer(input$add_layers)

          lapply(1:layerCount,function(i) {
            if (!is.null(input[[paste("show_aesthetic_control",i,sep="")]])) {
              id_val_map <- paste("map_aes_controls",i,sep="")
              id_val_manual <- paste("manual_control_group",i,sep="")
              id_val_reset <- paste("reset_button",i,sep="")
              if (input[[paste("show_aesthetic_control",i,sep="")]] == F) {
                shinyjs::hide(id = id_val_map)
                shinyjs::hide(id = id_val_manual)
                shinyjs::hide(id = id_val_reset)

              } else {
                shinyjs::show(id = id_val_map)
                shinyjs::show(id = id_val_manual)
                shinyjs::show(id = id_val_reset)

              }
            }

          })

      })

      ########### Generate Collapsable Panel ##############

      output$collapsible_panel <- renderUI({

        if (!is.null(input$reset_collapse_panel)) {
          div(id = paste("collapse_panel_controls", input$reset_collapse_panel,sep = ""),
              shinyBS::bsCollapse(id = "collapseExample", open = "Panel 2",
                         shinyBS::bsCollapsePanel("Themes & Plot Title",
                                         selectInput("select_theme","Select theme", choices = c("",as.character(unlist(themes)))),
                                         textInput("add_title","Add Plot Title", value = ""), style = "success"),
                         shinyBS::bsCollapsePanel("Axis Controls",
                                         checkboxInput("manual_adjust_range","Manually adjust axis range?", value = F),
                                         uiOutput("axis_controls")
                                         , style = "success"),
                         shinyBS::bsCollapsePanel("Facetting",
                                         uiOutput("facetting_var"), style = "success"),
                         shinyBS::bsCollapsePanel("Color Palettes (Discrete Variables)",
                                         radioButtons("color_scheme", "Select Color scheme", choices = c("Sequential","Qualitative","Diverging")),
                                         uiOutput("col_palette"),
                                         radioButtons("fill_color_select", "Fill or Color?", choices = c("Fill","Color"))
                                         , style = "success"),
                         shinyBS::bsCollapsePanel("Color Palettes (Contnuous Variables)",
                                         checkboxInput("revert_default","Revert to default?",value = F),
                                         shinyjs::colourInput("low_color", "Select low color", value="white"),
                                         shinyjs::colourInput("high_color", "Select high color", value="white")
                                         , style = "success"))
          )


        }


      })


      ########## Facetting ###############################

      output$facetting_var <- renderUI({


        selected_dataset <- eval(parse(text = input$get_dataset))
        vars <- colnames(selected_dataset)

        selectInput(inputId = "facet_var","Select a facetting variable",choices = c("",vars))

      })

      ########### Color palette #####################

      output$col_palette <- renderUI({

        palette_choices <- col_palette %>% filter(Type == input$color_scheme) %>% select(Palette)
        selectInput("color_palette", "Choose a color palette", choices = c("", as.character(unlist(palette_choices))))

      })

      ########### Axis controls #####################

      output$axis_controls <- renderUI({

        plot_template <- plot_with_themes()
        plot_out <- eval(parse(text = plot_template))
        selected_dataset <- eval(parse(text = input$get_dataset))

        genAxisSlider(selected_dataset,plot_out)

      })


      ######## Construct Plot ###############

      genPlot <- reactive({

        validate(
          need(try(input[[paste("plot_type",1)]] != ""), "Please select a variable")
        )

          layerCount <- as.integer(input$add_layers)
          layers <- list()

          layers <- lapply(1:layerCount,function(i) {

            # Extract Variable and Plot Type selection
            plot_type <- input[[paste("plot_type",i)]]
            x_var <- input[[paste("xaxis_var",i)]]
            y_var <- input[[paste("yaxis_var",i)]]

            if (!is.null(input[[paste("select_aes",i,sep="")]])) {
              if(input[[paste("aes_mapping_var",i,sep="")]] != "") {

                selected_aes <- input[[paste("select_aes",i,sep="")]]
                mapping_var <- input[[paste("aes_mapping_var",i,sep="")]]

              } else {
                selected_aes <- ""
              }


              # Extract Manual aesthetic control settings
              applicable_aes <- aesthetics %>% filter(Geom == input[[paste("plot_type",i)]]) %>% select(Aes)
              applicable_aes <- as.character(unlist(applicable_aes))
              manual_aes_control <- vector()
              manual_aes_value <- vector()

              manual_aes_control <- sapply(1:length(applicable_aes),function(j) {
                if (!is.null(input[[paste("manual_adjustment",i,j,sep="")]]) & input[[paste("manual_adjustment",i,j,sep="")]] != "#FFFFFF") {
                  if (input[[paste("manual_adjustment",i,j,sep="")]] != "") {
                    manual_aes_control[j] <- applicable_aes[j] } else {
                      manual_aes_control[j] <- ""
                    }

                } else {
                  manual_aes_control[j] <- ""
                }
              })

              manual_aes_value <- sapply(1:length(applicable_aes),function(j) {

                if (!is.null(input[[paste("manual_adjustment",i,j,sep="")]])) {
                  if (input[[paste("manual_adjustment",i,j,sep="")]] != "" & input[[paste("manual_adjustment",i,j,sep="")]] != "#FFFFFF") {
                    if (applicable_aes[j] %in% c('color','fill','linetype','position','method','stat')) {
                      manual_aes_value[j] <- paste("'", input[[paste("manual_adjustment",i,j,sep="")]] , "'",sep="")

                    } else {
                      manual_aes_value[j] <- input[[paste("manual_adjustment",i,j,sep="")]]
                    }

                  } else {
                    manual_aes_value <- ""
                  }

                } else {
                  manual_aes_value[j] <- ""
                }

              })

              # Remove aesthetics not manipulated
              manual_aes_control <- manual_aes_control[manual_aes_control != ""]
              manual_aes_value <- manual_aes_value[manual_aes_value != ""]

              manual_adjust_concat <- paste(manual_aes_control,manual_aes_value,sep="=")
              manual_adjust_concat <- paste(manual_adjust_concat,collapse = ",")


              if (selected_aes == "" & manual_adjust_concat == "") {
                layers[i] <- paste(plot_type,"(aes(",x_var,",",y_var,"))")
              } else if(selected_aes != "" & manual_adjust_concat == "") {
                layers[i] <- paste(plot_type,"(aes(",x_var,",",y_var,",",selected_aes,"=",mapping_var,"))")
              } else if(selected_aes == "" & manual_adjust_concat != "") {
                layers[i] <- paste(plot_type,"(aes(",x_var,",",y_var,")", ",",manual_adjust_concat,")")
              } else {
                layers[i] <- paste(plot_type,"(aes(",x_var,",",y_var,",",selected_aes,"=",mapping_var,")", ",",manual_adjust_concat,")")
              }


            } else {
              layers[i] <- paste(plot_type,"(aes(",x_var,",",y_var,"))")
            }

          })

          plot_template <- paste(layers,collapse = "+")
          plot_template <- paste("ggplot(dataset)",plot_template,sep="+")

          plot_template <- gsub("dataset",input$get_dataset,plot_template)
          plot_template

        })


      #### Add theme to plot ################

      plot_with_themes <- reactive({

        plot_template <- genPlot()

        if (!is.null(input$select_theme)) {
          if(input$select_theme != "") {
            plot_template <- paste(plot_template, "+", input$select_theme)
          } else {
            plot_template <- plot_template
          }
        }


        ### Add Title #############

        if(!is.null(input$add_title)) {
          if (input$add_title != "") {
            plot_template <- paste(plot_template, " + ggtitle(' ", input$add_title," ') ")

          }
        }

        ##### Add facet ###############

        if (!(is.null(input$facet_var))) {
          if(input$facet_var != "") {
            plot_template <- paste(plot_template, " + facet_wrap(~", input$facet_var, ")")
          }
          else { plot_template = plot_template}
        }

       ####### Color palette - Discreete ############

        if (!(is.null(input$color_palette))) {
          if(input$color_palette != "") {
            if(input$fill_color_select == "Fill") {
              plot_template <- paste(plot_template, " + scale_fill_brewer(palette='", input$color_palette, "') ",sep="")
            } else {
              plot_template <- paste(plot_template, " + scale_color_brewer(palette='", input$color_palette, "') ",sep="")
            }

          }
          else { plot_template = plot_template}
        }

      ########## Color palette - Continuous #############

        if (!is.null(input$low_color)) {
          if(input$low_color != '#FFFFFF' & input$revert_default != T) {
            plot_template <- paste(plot_template, " + scale_color_gradient(low='", input$low_color, "', high = '", input$high_color, "')", sep="")

          } else {
            plot_template <- plot_template
          }

        }

      ############### Axis Range ############################

        observeEvent(input$y_slider,{

          updateSliderInput(session,inputId = "y_slider", value = input$y_slider)

        })

        observeEvent(input$x_slider,{

          updateSliderInput(session,inputId = "x_slider", value = input$x_slider)

        })

      if (!is.null(input$y_slider)) {
        if (input$manual_adjust_range == T) {
          val_slider <- input$y_slider
          plot_template <- paste(plot_template, " + ylim(c(",val_slider[1],",",val_slider[2],"))")
        } else {
          plot_template <- plot_template
        }

      }

        if (!is.null(input$x_slider)) {
          if (input$manual_adjust_range == T) {
            val_slider <- input$x_slider
            plot_template <- paste(plot_template, " + xlim(c(",val_slider[1],",",val_slider[2],"))")
          } else {
            plot_template <- plot_template
          }

        }


       plot_template

      })

      ######## Generate Plot (static) #####

      output$plot_output <- renderPlot({

        validate(
          need(try(input$get_dataset != ""), "Please select a data set")
        )

        plot_template <- plot_with_themes()
        eval(parse(text = plot_template))

      })


      ########## Generate Plot (Interactive) ####

      output$interactive_plot_output <- plotly::renderPlotly({

        plot_template <- plot_with_themes()
        g <- eval(parse(text = plot_template))
        plotly::ggplotly(g)

      })


      ########## Display the underlying code ###########

      output$plot_code <- renderText({

        validate(
          need(try(input$get_dataset != ""), "Please select a data set")
        )

        if (input$switch_interactivity  == F) {
          plot_with_themes()
        } else {
          paste("plotly::ggplotly(",plot_with_themes(),")")

        }

      })

      ########### Display the enlarged plot view #########

      output$plot_out_zoom <- renderPlot({
        plot_template <- plot_with_themes()
        eval(parse(text = plot_template))

      }, height = 600)

      session$onSessionEnded(function() {
        stopApp()
  })


    }

)

}
