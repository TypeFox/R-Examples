#'Animate Stream Clustering.
#'
#'A function to plot data streams and clusterings. The visualisation is based on 
#'\link[shiny]{shiny} and \link[ggplot2]{ggplot}. Data is plotted as a
#'scatterplot matrix and individual scatterplots can be selected for a more
#'detailed view that includes tooltips. Please note that this function was
#'developed for the Streaming algorithms in the subspaceMOA package and may or may
#'not work for streams and clustering algorithms.
#'
#'@param dsc a DSC object representing the clustering of a data stream.
#'@param dsd a DSD object representing a data stream.
#'@param step the step size used in \link{animate_stream_interactive}. This 
#'  regulates how many points will be taken out of the stream, clustered and the
#'  plotted along with their clusters every time a step is performed.
#'@param delay time between two clustering steps
#'@param launch.browser will be passed on to \link[shiny]{runApp}, so that the 
#'  visualisation can be shown in e.g. RStudio's Viewer pane, if this is 
#'  desired.
#'@export
#'@import ggplot2
#'@import shiny
#'@import magrittr
#'@import stream
#'@import streamMOA
 animate_stream_interactive <- function(dsc,dsd,step=1500,delay=10000,launch.browser=getOption("shiny.launch.browser",interactive())) {
   #Create a shiny UI in which to display the streaming data
   ui <- makeUI(show_animate_buttons=T)
   server <- makeServer(dsc,dsd,step,delay=delay)
   
   onStart <- function(){}
   app <- shinyApp(ui=ui,server=server,onStart=onStart)
   runApp(app,launch.browser=launch.browser)
 }
#'Show Stream Clustering.
#'
#'A non-animated version of \link{animate_stream_interactive}.
#'
#'@param dsc a DSC object representing the clustering of a data stream.
#'@param points a \link{data.frame} of points that will be plotted along with
#'  the clustering.
#'@param launch.browser will be passed on to \link[shiny]{runApp}, so that the 
#'  visualisation can be shown in e.g. RStudio's Viewer pane, if this is 
#'  desired.
#'@export
plot_stream_interactive <- function(dsc,points,launch.browser=getOption("shiny.launch.browser",interactive())) {
  ui <- makeUI(show_animate_buttons=F)
  server <- makeServer(dsc,points)
  onStart <- function(){}
  app <- shinyApp(ui=ui,server=server,onStart=onStart)
  runApp(app,launch.browser=launch.browser)
}

makeUI <- function(show_animate_buttons) {
  ui <- fluidPage(
    #This dummy input exists because a conditional panel can only depend on
    #values in the input or the output object, so we encode part of the 
    #application's state in an invisible selectInput. This is, of course, a very
    #horrible way of doing it, but it works.
    conditionalPanel("false",
                     selectInput("dummyInput",
                                 label="You should not be seeing this",
                                 choices=c("matrix","detail"))),
    conditionalPanel("input.dummyInput == 'matrix'",
                     plotOutput("plot_matrix",click="plot_matrix_click",
                                width="95%",
                                height="600px"),
                     conditionalPanel(r_logical_to_js_boolean_string(show_animate_buttons),
                      fluidRow(
                        column(4,actionButton(inputId="stop_button",label="Stop",
                                                class="btn-danger btn-large btn-block")),
                        column(4,actionButton(inputId="step_button",label="Step",
                                                class="btn-primary btn-large btn-block")),
                        column(4,actionButton(inputId="run_button",label="Run",
                                                class="btn-success btn-large btn-block"))
                      ))),
    conditionalPanel("input.dummyInput=='detail'",
                     fluidRow(plotOutput("detail_plot",
                                         hover="detail_plot_hover",
                                         height="600px")),
                     #The button to go back to the plot matrix view
                     fluidRow(
                       actionButton(inputId="back_button",
                                    label="",icon=icon(name="th"),
                                    class="btn-primary btn-large btn-block")
                     ),
                     #The text field in which information on the point that is hovered
                     #over is given.
                     fluidRow(
                       wellPanel(htmlOutput("tooltip"))
                     )
    )
  )
  return(ui)
}
#Creates a Shiny server to handle the logic of which plots are shown
makeServer <- function(dsc,dsd,step=NULL,delay=5000) {
  
  if(is.data.frame(dsd)) {
    points <- dsd
    initial_data_frame <- format_data_from_dsc(dsc,points=points)
    if(is.null(points[["class"]])) {
      number_of_dimensions <- ncol(points)
    } 
    else {
      number_of_dimensions <- ncol(points)-1
    }
  } else {
    initial_data_frame <- format_data_from_dsc(dsc)
    #Try to get the number of dimensions of the stream
    number_of_dimensions <- dsd[["d"]]
    #If that failed just take one point and find out how many dimensions the stream data has
    if(is.null(number_of_dimensions)) {
      number_of_dimensions <- ncol(get_points(dsd,1,cluster=F,class=F))
    }
  }
  

  server <- function(input,output,session){
    #A reactiveValues object to keep track of global application state. In this
    #case, we keep track of whether we are showing the clustering as a
    #scatterplot Matrix (display_mode=="matrix") or a detailed view in which two
    #dimensions are plotted against each other (display_mode=="detail") 
    #Additionally we are keeping track of the data frame that is currently being
    #shown as well as whether we are currently running the stream clustering
    #continuously.
    state <- reactiveValues(display_mode="matrix",current_data_frame=initial_data_frame,
                            running=F,
                            should_perform_step=F,
                            plot_was_recently_drawn=F)
    #When the "back" button (looks like a grid of squares) is pressed, the
    #display mode should be set to "matrix"
    observeEvent(input$back_button, {
      state$display_mode  <- "matrix"
    })  
    observeEvent(input$step_button, {
      state$should_perform_step <- T
    })
    observeEvent(input$run_button, {
      state$running <- T
    })
    observeEvent(input$stop_button, {
      state$running <- F
    })
    #If the current state of the app is running then start performing a step and
    #repeat this action after "delay"
    observe({
       if(state$running) {
         isolate({state$should_perform_step <- T})
         invalidateLater(delay,session)
       }
    })
    #Whenever a step should be performed, get new data from the
    #stream and push the data into the clusterer.
    observe({
      if(state$should_perform_step) {
        isolate({
          new_points <- get_points(dsd,step,class=T)
          withProgress({
            update(dsc,DSD_Memory(new_points[,1:(ncol(new_points)-1)]),step)
          },message="Updating the Clustering")
          res <- format_data_from_dsc(dsc,points=new_points)
          state$current_data_frame  <-  res
        })
        state$should_perform_step <- F
      }
    })
    #Always have the selected value in the dummyInput reflect the current state of the application
    #This expression is executed every time state$display_mode changes.
    observe({
      updateSelectInput(session=session,inputId="dummyInput",selected=state$display_mode)
    })
    #Keep track of the last plot that was clicked on in the scatterplot matrix.
    #Changes whenever the main plot is clicked on.
    #A helper function in helper.R determines which plot was being clicked on.
    last_plot_clicked_on_in_matrix <- reactive({
      c <- input$plot_matrix_click
      if(is.null(c)) {
        return(NULL)
      } else {
        state$display_mode <- "detail"
        return(from_coords_to_plot(x=c$x,y=c$y,domain=c$domain,number_of_dimensions=number_of_dimensions))
      }
    })
    #Make sure that last_plot_clicked_on_in_matrix is always current whenever a
    #click occurs. If it weren't for this observe block, its value would get 
    #recomputed only when the detail plot is being shown, which is not the behavior
    #we want because it keeps the detail plot from ever being shown.
    observe({
      last_plot_clicked_on_in_matrix()
    })
    #Draw the plot matrix from the current data frame
    output$plot_matrix <- renderPlot({
      withProgress({
        list_of_plots <- create_plot_matrix(state$current_data_frame)
        incProgress(message="Creating Scatterplot Matrix")
        plotmatrix <- make_plot_matrix(list_of_plots,ncol=number_of_dimensions)
        incProgress(amount = 0.3,message="Displaying Scatterplot Matrix")
        if(is.null(plotmatrix)) return()
        grid::grid.draw(plotmatrix)
      },message="Creating Scatterplots")
      isolate({state$plot_was_recently_drawn <- T})
    })
    output$detail_plot <- renderPlot({
      if(!is.null(last_plot_clicked_on_in_matrix())) {
        res <- state$current_data_frame %>%
          basic_plot_from_dataframe(last_plot_clicked_on_in_matrix()) %>%
          style_plot_for_detail()
        isolate({state$plot_was_recently_drawn <- T})
        res
      }
    })
    output$tooltip <- renderPrint({
      row <- row_from_two_values(dataframe=state$current_data_frame,
                                 hover_list=input$detail_plot_hover,
                                 area_around_cursor=0.05)
      res <- dataframe_row_to_html(row)
      cat(res)
    })
  }
  return(server)
}

