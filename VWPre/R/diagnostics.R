#' Plots diagnostic plots of subject/item variance.
#' 
#' \code{plot_var_app} calculates and plots within-subject/item standard deviation,
#' along with standardized by-subject/item means for a given interest area, within
#' a given time window.
#' 
#' @export
#' @import dplyr
#' @import lazyeval
#' @import ggplot2
#' @import shiny
#' 
#' @param data A data table object output by either \code{\link{bin_prop}}. 
#' \code{\link{transform_to_elogit}}, or \code{\link{create_binomial}}.
#' @examples
#' \dontrun{
#' library(VWPre)
#' # For plotting the grand average with the included theme
#' plot_var_app(data = dat) 
#' }
plot_var_app <- function(data = data) {
  dat <- data
  
  shiny::shinyApp(
    ui =  # Use a fluid Bootstrap layout
      shiny::fluidPage(    
        
        # Give the page a title
        shiny::titlePanel("Variability"),
        
        # Generate a row with a sidebar
        shiny::sidebarLayout(      
          
          # Define the sidebar with one input
          shiny::sidebarPanel(
            shiny::selectInput("type", "Group:", 
                        choices=c("Subjects", "Item")),
            shiny::selectInput("scale", "Input:", 
                                       choices=c("Proportions", "Empirical Logits")),
                           shiny::conditionalPanel(
                             condition = "input.scale == 'Proportions'",
                             shiny::selectInput(
                               "PCol", "Interest Areas", c(intersect(grep("_P",names(data), value=TRUE),
                                                                    grep("IA_",names(data),value=TRUE)), "Select"),
                               selected = "Select")),
                           shiny::conditionalPanel(
                             condition = "input.scale == 'Empirical Logits'",
                             shiny::selectInput(
                               "ECol", "Interest Areas", c(intersect(grep("_E",names(data), value=TRUE),
                                                                    grep("IA_",names(data),value=TRUE)), "Select"),
                               selected = "Select")),
            shiny::sliderInput("mintime", 
                        "Min Time:", 
                        value = min(dat$Time),
                        min = min(dat$Time), 
                        max = max(dat$Time)),
            shiny::sliderInput("maxtime", 
                        "Max Time:", 
                        value = max(dat$Time),
                        min = min(dat$Time), 
                        max = max(dat$Time))
          ),
          
          # Create a spot for the plot
          shiny::mainPanel(
            shiny::plotOutput("Plot")  
          )
          
        )
      ),
    
    server = function(input, output) {
	
		dat<- shiny::reactive({
        data <- data
        return(data)
		})
        Col <- shiny::reactive({
        if (input$scale=='Proportions') {
          col <- input$PCol
        } else if (input$scale=='Empirical Logits') {
          col <- input$ECol
        } 
        return(col)
		})
		
      output$Plot <- shiny::renderPlot({
        
        # Render a plot
		
		Col <- Col()
		
		if (Col == "Select") {
          warning("Please select an interest area.")
        } else {
		
        if (input$type == "Item") {
          
          dat1 <- dat() %>% filter(Time >= input$mintime & Time <= input$maxtime) %>%
            group_by(Item) %>%
            summarise_(Avg = interp(~mean(Col), Col = as.name(Col)), StDev = interp(~sd(Col), Col = as.name(Col))) %>%
            ungroup() %>%
            mutate(., Zscore = (Avg-mean(Avg))/sd(Avg))
          
          ggplot(dat1, aes(Item, Zscore)) + 
            geom_segment(aes(x = Item, y = 0, xend = Item, yend = Zscore)) +
            geom_point(aes(size = StDev), shape = 21, fill = "gray", alpha = 0.75) +
            geom_hline(yintercept=0) + 
            geom_hline(yintercept=2.5, color = "gray", linetype = 2) +
            geom_hline(yintercept=-2.5, color = "gray", linetype = 2) + 
            labs(y = "Z-score of looks") +
            scale_size(name="Within\nItem SD") +
            theme_bw() +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  axis.text.x  = element_text(angle=90, vjust=0.5, size=8)
            )
          
        } else {
          
          dat1 <- dat() %>% filter(Time >= input$mintime & Time <= input$maxtime) %>%
            group_by(Subject) %>%
            summarise_(Avg = interp(~mean(Col), Col = as.name(Col)), StDev = interp(~sd(Col), Col = as.name(Col))) %>%
            ungroup() %>%
            mutate(., Zscore = (Avg-mean(Avg))/sd(Avg))
          
          ggplot(dat1, aes(Subject, Zscore)) + 
            geom_segment(aes(x = Subject, y = 0, xend = Subject, yend = Zscore)) +
            geom_point(aes(size = StDev), shape = 21, fill = "gray", alpha = 0.75) +
            geom_hline(yintercept=0) + 
            geom_hline(yintercept=2.5, color = "gray", linetype = 2) +
            geom_hline(yintercept=-2.5, color = "gray", linetype = 2) + 
            labs(y = "Z-score of looks") +
            scale_size(name="Within\nSubject SD") +
            theme_bw() +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  axis.text.x  = element_text(angle=90, vjust=0.5, size=8)
            )
        }
        
        }
      })
    }
  )
}



#' Plots diagnostic average plots of subjects/items.
#' 
#' \code{plot_indiv_app} calculates and plots interest area averages for a 
#' given subject/item.
#' 
#' @export
#' @import dplyr
#' @import tidyr
#' @import lazyeval
#' @import ggplot2
#' @import shiny
#' 
#' @param data A data table object output by either \code{\link{bin_prop}}. 
#' \code{\link{transform_to_elogit}}, or \code{\link{create_binomial}}.
#' @examples
#' \dontrun{
#' library(VWPre)
#' # For plotting the grand average with the included theme
#' plot_indiv_app(data = dat)
#' } 
plot_indiv_app <- function(data = data) {
  shiny::shinyApp(
    ui =  # Use a fluid Bootstrap layout
      shiny::fluidPage(    
        
        # Create a spot for the plot
        shiny::mainPanel("",
                  shiny::fluidRow(
                    shiny::plotOutput("Indiv")#, width = "100%", height = "500px")
                  ),
                  
                  # Give the page a title
                  shiny::titlePanel("Individual Averages"),
                  
                  # Generate a row with a sidebar
                  shiny::hr(),
                  
                  shiny::fluidRow(
                    shiny::column(4, offset = 0,
                           shiny::selectInput("scale", "Scale:", 
                                       choices=c("Proportions", "Empirical Logits")),
                           shiny::conditionalPanel(
                             condition = "input.scale == 'Proportions'",
                             shiny::selectizeInput(
                               "PCols", "Interest Areas", intersect(grep("_P",names(data), value=TRUE),
                                                                    grep("IA_",names(data),value=TRUE)),
                               selected = NULL, multiple = TRUE,
                               options = list(placeholder = 'select interest areas'))),
                           shiny::conditionalPanel(
                             condition = "input.scale == 'Empirical Logits'",
                             shiny::selectizeInput(
                               "ECols", "Interest Areas", intersect(grep("_E",names(data), value=TRUE),
                                                                    grep("IA_",names(data),value=TRUE)),
                               selected = NULL, multiple = TRUE,
                               options = list(placeholder = 'select interest areas')))
                    ),
                    shiny::column(4, offset = 0,
                           shiny::selectInput("type", "Group:", 
                                       choices=c("Subjects", "Items")),
                           shiny::conditionalPanel(
                             condition = "input.type == 'Items'",
                             shiny::selectInput("item", "Plot:", 
                                         choices=unique(levels(data$Item)))),
                           shiny::conditionalPanel(
                             condition = "input.type == 'Subjects'",
                             shiny::selectInput("subj", "Individual:", 
                                         #selected = unique(levels(data$Subject))[3],
                                         choices=unique(levels(data$Subject))))
                    ),
                    shiny::column(4, offset = 0,
                           shiny::sliderInput("mintime", 
                                       "Min Time:", 
                                       value = min(data$Time),
                                       min = min(data$Time), 
                                       max = max(data$Time)),
                           shiny::sliderInput("maxtime", 
                                       "Max Time:", 
                                       value = max(data$Time),
                                       min = min(data$Time), 
                                       max = max(data$Time))
                    )
                    
                  )
        )
      ),
    
    server = function(input, output) {
      
      granddata<- shiny::reactive({
        data <- data
        #data <- as.data.frame(data)
        # data <- data[data$Subject==input$subj,]
        #data <- filter(data, Subject==input$indiv)
        return(data)
      })
      cols <- shiny::reactive({
        if (input$scale=='Proportions') {
          cols <- input$PCols 
        } else if (input$scale=='Empirical Logits') {
          cols <- input$ECols
        } 
        return(cols)
      })
      SCALE <- shiny::reactive({
        SCALE <- input$scale
        return(SCALE)
      })
      YLIM <- shiny::reactive({
        if (input$scale=='Proportions') {
          YLIM <- c(0,1) 
        } else if (input$scale=='Empirical Logits') {
          YLIM <- c(-4,4)
        } 
        return(YLIM)
      })
      SN <- shiny::reactive({
        SN <- c("Time", cols())
        return(SN)
      })
      INDIV <- shiny::reactive({
        if (input$type=='Subjects') {
          INDIV <- input$subj 
        } else if (input$type=='Items') {
          INDIV <- input$item 
        } 
        return(INDIV)
      })
      indivdata<- shiny::reactive({
        if (input$type=='Subjects') {
          data <- data[data$Subject==input$subj ,]
        } else if (input$type=='Items') {
          data <- data[data$Item==input$item ,]
        } 
        return(data)
      })
      
      # Render a plot
      output$Indiv <- shiny::renderPlot({
        
        ylim <- YLIM()
        sel_names <- SN()
        Cols <- cols()
        Ind <- INDIV()
        scale <- SCALE()
        
        if (is.null(Cols)) {
          warning("Please select interest areas.")
        } else {
          
          GrandAvg <- granddata() %>%
            select(match(sel_names,names(.))) %>%
            tidyr::gather_("IA", "VALUE", unique(Cols), na.rm = FALSE, convert = FALSE) %>%
            group_by_("IA", "Time") %>%
            summarise(mean = mean(VALUE, na.rm = T), se = sd(VALUE) / sqrt(length(VALUE))) %>%
            mutate(alpha = 0.25, group = "Grand Average")
          GrandAvg$group <- as.factor(GrandAvg$group)
          
          IndivAvg <- indivdata() %>% select(match(sel_names,names(.))) %>%
            tidyr::gather_("IA", "VALUE", unique(Cols), na.rm = FALSE, convert = FALSE) %>%
            group_by_("IA", "Time") %>%
            summarise(mean = mean(VALUE, na.rm = T), se = sd(VALUE) / sqrt(length(VALUE))) %>%
            mutate(alpha = 1, group = "Individual Average")
          IndivAvg$group <- as.factor(IndivAvg$group)
          
          Avg <- rbind(IndivAvg, GrandAvg)
          
          if (input$scale == "Empirical Logits") {
            ylim[1] <- min(Avg$mean) - 0.25
            ylim[2] <- max(Avg$mean) + 0.25
          } else {
            ylim <- c(0,1)
          }
          
          ggplot(Avg, aes(x = Time, y = mean, colour = IA)) + 
            geom_point(alpha=0.7) +
            geom_line(alpha=0.7) +
            geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .3, alpha=0.7) +
            facet_grid(. ~ group) +
            ylab(scale) +
            scale_x_continuous(limits = c(input$mintime, input$maxtime)) + 
            scale_y_continuous(limits = c(ylim[1], ylim[2])) +
            scale_colour_brewer(palette = "Set1") +
            theme_bw() +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank())
          
        }
        
        
      })
      
    }
  )
}