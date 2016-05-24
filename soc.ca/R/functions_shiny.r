#' Explore the cloud of individuals
#'
#' @param object a a soc.ca class object as created by \link{soc.mca} and 
#'   \link{soc.csa}
#' @param active Defines the active modalities in a data.frame with rows of individuals and columns of factors, without NA's' 
#' @param sup Defines the supplementary modalities in a data.frame with rows of individuals and columns of factors, without NA's 
#'
#' @return an html application
#' @export
#'
#' @examples
#' \dontrun{
#' example(soc.mca)
#' ind.explorer(result, active, sup)
#' }

ind.explorer     <- function(object, active, sup = NULL){
  
  # -----------
  # Creating data
  
  data              <- data.frame("Dimension" = object$coord.ind, active)
  indicator.matrix  <- object$indicator.matrix
  if(!is.null(sup)) indicator.matrix  <- cbind(indicator.matrix, indicator(sup))
  if(!is.null(sup)) data              <- cbind(data, sup)
  rownames(data)    <- object$names.ind
  
  # ------------
  # UI 
  
  ui  <- shinyUI(
    fluidPage(
      titlePanel(
        "Explore"
        ), 
      fluidRow(
        column(3, selectInput(inputId = "dim.plot.x", label = "Dim X-axis", choices = 1:object$nd, selected = 1, width = 100),
                  selectInput(inputId = "dim.plot.y", label = "Dim Y-axis", choices = 1:object$nd, selected = 2, width = 100),
                  selectInput("fill", "Fill", choices = colnames(indicator.matrix), width = 150),
                  checkboxInput(inputId = "ellipse", "Draw ellipse", value = FALSE),
                  checkboxInput(inputId = "density", "Draw density", value = FALSE),
                  checkboxInput(inputId = "labels", "Modality labels", value = TRUE)
               ),
        column(8, plotOutput("map.ind", click = "plot_click"))
        
      ),
      fluidRow(
        column(3),
        column(8, plotOutput("map.mod", height = 600))
      ),
      fluidRow(
        column(3),
        column(8, tableOutput("info"))
        )
  )
  )
  
  # ---------------
  # Server
  server <- shinyServer(function(input, output) {
    
    # Click
    output$map.ind <- renderPlot({
      
      dim.plot       <- c(as.numeric(input$dim.plot.x), as.numeric(input$dim.plot.y))
      fill.ind       <- subset(indicator.matrix, select = input$fill)
      fill.var       <- factor(fill.ind, levels = c(NA, 1), labels = c(input$fill))
      
     # if(input$csa) object         <- soc.csa(object, class.indicator = which(fill.ind == 1), sup = sup)
      
      p <- map.ind(object, dim = dim.plot, point.fill = fill.var) + coord_cartesian()
      p <- p + scale_fill_manual(values = c("black", "white"), na.value = "white")
      if(input$ellipse == TRUE) p <- map.ellipse(object, ca.plot = p, variable = fill.var)
      if(input$density == TRUE) p <- map.density(object, map = p)
      p
    })
    
    output$info <- renderTable({
      np <- nearPoints(data, input$plot_click, xvar = "Dimension.1", yvar = "Dimension.2", threshold = 10)
      if(is.null(input$plot_click)) np <- data[1,]
      t(np[, -(1:object$nd)])
    })
    
    output$map.mod <- renderPlot({
      dim.plot       <- c(as.numeric(input$dim.plot.x), as.numeric(input$dim.plot.y))
      np <- nearPoints(data, input$plot_click, xvar = "Dimension.1", yvar = "Dimension.2", threshold = 10)
      if(is.null(input$plot_click)) np <- data[1,]
      ind           <- which(object$names.ind %in% rownames(np))
      mod.names     <- colnames(indicator.matrix)
      mod.names.ind <- mod.names[colSums(indicator.matrix[ind, , drop = FALSE] == 1) > 0]
      active.mod.ind <- object$names.mod %in% mod.names.ind
      sup.mod.ind   <- object$names.sup %in% mod.names.ind
      map.title     <- paste(rownames(np), collapse = " & ")
      variable      <- c(object$variable[active.mod.ind], object$variable.sup[sup.mod.ind])
      
      map.select(object, dim = dim.plot, list.mod = active.mod.ind,
                 list.sup = sup.mod.ind, map.title = map.title,
                 label = input$labels, label.repel = TRUE,
                 label.fill = variable, point.color = variable)
    })})
  
  
  # ---------------
  # runApp
  
  app <- list(ui = ui, server = server)
  runApp(app)
}

# coords.ind.plot <- function(object, ind, dim = 1:9){
#              
# Name       <- object$names.ind[ind]  
# coords     <- data.frame(Name, Dimension = object$coord.ind[ind, dim])
# coords     <- melt(coords, id.vars = "Name")
# if(length(ind) == 1) coords$variable <- dim
# 
# ctr        <- data.frame(Name, Ctr = object$ctr.ind[ind, dim])
# ctr        <- melt(ctr, id.vars = "Name")
# gd         <- data.frame(coords, Contribution = ctr$value)
# gd$Name    <- as.factor(gd$Name)
# 
# levels(gd$variable) <- dim
# 
# p <- ggplot(gd, aes(x = variable, y = value, size = Contribution, fill = Name, group = Name)) 
# p <- p + geom_line(aes(color = Name), size = 0.5) + geom_point(shape = 21)
# p + xlab(label = "Dimension") + ylab(label = "Coordinate") + theme_minimal()
# }
