#' Interactive Backtest Graphics
#' 
#' This function takes a data frame and returns an interactive interface for the
#' backtest data. The interface contains drop-down menus for different
#' strategies, portfolios and instruments/sectors. The user can slice her data
#' set according to the strategies for selecting instruments or the portfolio
#' numbers of her individual portfolios. The user can also look at a specific
#' group of instruments according to the sector of instruments. These sectors
#' will be in the same dropdown menu as instruments. The elements in each of the
#' three dropdown menus are generated according to the input data set. The
#' dropdown menus for strategies and portfolio numbers will only contain
#' "Strategy Summary" or "Portfolio Summary" if there's no strategy or portfolio
#' information available in the input data set.
#' 
#' The user can also slice the data set according to the division of the data 
#' set into main strategies, and the division of any main strategy into small 
#' substrategies. The "substrategy" column of the data set should indicate the
#' division of main strategies, if there are any substrategies. Note that the
#' structure of substrategies should conform to main strategies. That is, a
#' single substrategy should be under only one main strategy.
#' 
#' Summary statistics are displayed in tables under the drop-down menu. The 
#' summary statistics contain data for the backtesting time horizon, profits, 
#' returns, volatility, market values, best/worst performers and biggest
#' drawdowns. Different interactive plots for the cumulative P&L, daily P&L, Net
#' Market Value (NMV), Gross Market Value (GMV) and number of contracts are
#' shown in the right panel.
#' 
#' The \code{backtestGraphics} function takes in a data set as well as the names
#' of essential columns in the data set, if these names are different from the 
#' default ones. The column names should be assigned to the corresponding 
#' variables so that the function can recognize these columns. If some essential
#' columns are missing, the function will try to fill in these columns with 
#' existing data. If the function fails to do so due to a lack of data, the 
#' function will return an error about the missing columns. If a column name is 
#' wrong, the function will treat that column as a missing column. The input 
#' data set has to contain either a name column or an ID column, a date column, 
#' either a net market value column or a number of contract column and a P&L 
#' column. The other columns either help the user slice the data set into 
#' different pieces to diagnose the data better, or can be calculated from the 
#' data inside the input data set.
#' 
#' @param x is a data frame that contains the necessary information.
#' @param trade.freq is the trading frequency of the data set, numeric type. The
#'   variable uses the number of dates between two trading dates to indicate 
#'   trading frequency. The default value of this variable is NULL, and in this 
#'   case the function will automatically calculate the trading period.
#' @param name.var is the column name of the instrument name column in the input
#'   data set, character type. The default value of this variable is "name". The
#'   user has to specify name.var if she passes in a data frame with a different
#'   column name for the instrument name. The input data set must contain either
#'   a name column or an ID column.
#' @param id.var is the column name of the instrument ID column, character type.
#'   The default value of this variable is "id". If the input data set labels
#'   the column of instrument ID's with some other column name, the user has to
#'   pass in the column name here. The input data set has to contain either a
#'   name column or an ID column.
#' @param date.var is the column name of the date column, character type. The 
#'   default value of this variable is "date". This column has to exist, and the
#'   column name has to be correct in order for the function to process the data
#'   set properly.
#'   
#' @param nmv.var is the column name of the "net market value" column, character
#'   type. The default value of this variable is "nmv". The input data set has
#'   to contain either a "net market value" column or a "number of contract" 
#'   column.
#' @param gmv.var The column name of the "gross market value" column if exists, 
#'   character type. The default value of this variable is "gmv". Such column
#'   will be automatically calculated from the "net market value" column if it
#'   does not exist.
#' @param pnl.var The column name of the profit-and-loss column, character type.
#'   The default value of this variable is "pnl". The data set has to contain 
#'   such column so that the function can function properly. If such column is 
#'   missing, the function will return an error indicating the problem.
#' @param contract.var The column name of the contract number column, character 
#'   type. The default value of this variable is "contract". If such column is
#'   missing, the function will instead use net market value as the contract
#'   number for each day.
#' @param capital.num The constant number of allocated capital for the whole 
#'   portfolio, numeric type. The default value of such variable is \code{NULL}.
#'   The function will use the number of allocated capital to calculate return 
#'   rates. If the allocated capital for the portfolio is not specified, the 
#'   function will use each day's gross market value to calculate the return 
#'   rate for each day.
#'   
#' @param sector.var is the column name of the sector column, character type.
#'   The default value of this variable is "sector". The sector column helps the
#'   user to group instruments into big groups according to industries. The
#'   function will still perform properly if the column name for sector is
#'   wrong.
#' @param strategy.var The column name of the strategy column, if any. Character
#'   type. The default value of this variable is "strategy". This column can be 
#'   missing from the input data set.
#' @param substrategy.var The column name of the substrategy column, if any. 
#'   Character type. The default value of this variable is "substrategy". This 
#'   column can be missing from the input data set.
#' @param portfolio.var The column name of the portfolio number column, if any. 
#'   Character type. The default value of this variable is "portfolio". This 
#'   column can be missing from the input data set.
#'   
#' @return a Shiny interface. The interactive shiny interface displays P&L and 
#'   cumulative P&L on one chart, and number of contracts, net market value and 
#'   gross market value on the other chart. Some summary statistics about return
#'   and performances are displayed on the left sidebar as tables.
#'   
#' @examples
#' \dontrun{
#' backtestGraphics(data = credit.bt)
#' }
#' 
#' @importFrom xts xts
#' @importFrom scales comma_format
#' @import shiny
#' @importFrom dplyr select filter mutate arrange group_by ungroup summarise
#'   left_join tbl_df
#' @import dygraphs
#' @importFrom stats sd
#'   
#' @export

backtestGraphics <- function(x,
                             trade.freq      = NULL,
                             name.var        = "name",   
                             id.var          = "id",
                             date.var        = "date",
                             nmv.var         = "nmv",
                             gmv.var         = "gmv",
                             pnl.var         = "pnl",
                             contract.var    = "contract",
                             capital.num     = NULL,
                             sector.var      = "sector",
                             strategy.var    = "strategy",
                             substrategy.var = "substrategy",
                             portfolio.var   = "portfolio"
                             ) {
  
  ## Initialize some variables that are fixed in the following helper functions
  ## so that CMD Check won't release notes
  
  nmv <- sector <- id <- gmv <- pnl <- contract <- nmv.temp <- name <- portfolio <- strategy <- 
    substrategy <- gmv.temp <- pnl.temp <- contract.temp <- NULL
  
  ## Check if the necessary columns are in the data set, and fill in the missing
  ## columns according to the existent columns. Then change the column names to
  ## some fixed names that are recognizable by all helper functions
  
  x <- cleanup_column(x,
                      name.var        = name.var,   
                      id.var          = id.var,
                      date.var        = date.var,
                      nmv.var         = nmv.var,
                      gmv.var         = gmv.var,
                      pnl.var         = pnl.var,
                      contract.var    = contract.var,
                      sector.var      = sector.var,
                      strategy.var    = strategy.var,
                      substrategy.var = substrategy.var,
                      portfolio.var   = portfolio.var)
  
  shinyApp(
    
    ######################################################
    ## The first part of the function is the user interface. On the side bar 
    ## panel, we create three dropdown menus for the user to slice her data frame,
    ## and all the elements in the dropdown menus are generated from the input 
    ## data frame. Below the dropdown menus is a button to conduct graphing and
    ## calculation, and then tables of summary statistics. The main panel
    ## contains the plots.
    ######################################################
    
    ui = pageWithSidebar(
      
      headerPanel(paste("Backtest Dashboard")),
      
      sidebarPanel(
        
        ## Create dropdown menu for strategy and portfolio. The dropdown menus 
        ## include strategies, portfolios and instruments, and the instrument
        ## menu has sectors as well if available. Substrategies have to conform
        ## to strategies, otherwise the slicing of data st might be wrong. 
        ## The strategy and the portfolio menus will only display "Strategy
        ## Summary" or "Portfolio Summary" if the strategy or the portfolio
        ## columns are not available because the "unique" function will return
        ## NULL.
        
        selectInput("strategy",
                    label    = "Strategies",
                    choices  = c("Strategy Summary", unique(x$strategy), unique(x$substrategy)),
                    selected = "Strategy Summary"),
        
        selectInput("portfolio",
                    label    = "Portfolios",
                    choices  = c("Portfolio Summary", unique(x$portfolio)),
                    selected = "Portfolio Summary"),
        
        ## If the number of unique id's in the data frame is more than 50, then
        ## create input-search textbox for portfolio summary, different sectors and commodities 
        
        selectizeInput("instrument", 
                       label   = "Instruments",
                       choices = c("Instrument Summary", unique(x$sector), unique(x$id)), 
                       options = list(maxOptions = 50)                      
        ),
        
        ## add fluidRow and columns so that space is more compactly used
        
        fluidRow(
          
          ## Add an actionButton so that Shiny will evaluate inputs only when
          ## user click this button this way Shiny won't eagerly evaluate when
          ## the user is still typing
          
          column(4,
                 
                 ## Unless submitButton is pressed, shiny won't evaluate inputs
                 
                 actionButton(inputId = "visualize", label = "Visualize!"),
                 
                 ## print additional blank row to place buttom in better position
                 
                 h3("   ")
          )
        ),
        
        ## Create texts and display summary statistics in summary tables
        
        tabsetPanel(
          tabPanel("Summary", p(" "), tableOutput("calc")),
          tabPanel("Detail",  p(" "), 
                   p("Top Three Drawdowns:"),    tableOutput("drawdown"),
                   p("Three Best Performers:"),  tableOutput("best3Performer"),
                   p("Three Worst Performers:"), tableOutput("worst3Performer")
          )
        ) 
      ),
      
      ## Create a column to embed the plots and the radio buttons to select
      ## plots. The graphs and the radio buttons are aligned at the center of
      ## the column.
      
      mainPanel(
        fluidRow(
          column(12, align = "center",
                 
                 ## Use dygraph to draw the first plot
                 
                 dygraphOutput("plot1"),
                 
                 ## Create radio buttons to choose between cumulative P&L and daily P&L
                 
                 div(style = "text-align: center",      
                     radioButtons("radio1", 
                                  label    = "", 
                                  choices  = list("Cumulative P&L" = "cumpnl", 
                                                  "Point-In-Time P&L" = "pnl"), 
                                  selected = "cumpnl", 
                                  inline   = TRUE )),
                 
                 ## Use dygraph to draw the second plot
                 
                 dygraphOutput("plot2"),
                 
                 ## Create radio buttons to choose among nmv, gmv and number of contracts
                 
                 div(style = "text-align: center", 
                     radioButtons("radio2", 
                                  label    = "", 
                                  choices  = list("NMV" = "nmv", "GMV" = "gmv",
                                                  "Number of Contracts" = "contract"),
                                  selected = "nmv", 
                                  inline   = TRUE)))
        )
      )
    ),
    
    #####################################################
    ## The second part of the backtestGraphics function is here. This part looks
    ## at the user's choices on the dropdown menus and the radio buttons. Then,
    ## the function calles a bunch of helper functions to slice the data frame
    ## into different pieces and calculates the summary statistics for the
    ## selected piece. Note that all actual operations are done outside the
    ## backtestGraphics function, and the function only calles helper functions
    ## to process data, and then arrange and format all outputs into tables. The
    ## graphing functions are also called in this part.
    #####################################################
    
    server <- function(input,output,session){
      
      ## set NULL to avoid NOTES in CMD Check... :3
      
      name <- portfolio <- strategy <- substrategy <- NULL
      
      ## Subtract data for each day's total net market value, total
      ## profit & loss and total contract number.
      
      intermediate <- eventReactive(
        input$visualize, withProgress(message = "Parsing Data", value = 0, {
          f <- slice_data(x, input, capital.num)
          
          ## Look at the dates in the data frame and figure out the trading frequency.
          
          if(is.null(trade.freq)){
            f$trade.freq <- trade_freq(f$x)
          } else {
            f$trade.freq <- trade.freq
          }
          
          return(f)
        })
      )
      
      ## Calculate all the summary statistics of the data frame
      
      valuefun <- eventReactive(
        input$visualize, withProgress(message = "Summarizing", value = 0, {
          x.intermediate <- intermediate()
          f <- stat_calculation(x.list = x.intermediate)
          return(f)     
        })
      )
      
      ## Create a table to combine all summary statistics together. All the 
      ## numbers are also formatted here with roundings and commas
      
      output$calc <- renderTable({
        
        x.intermediate <- intermediate()
        stat.value <- valuefun()
        
        items <- c("Start Date", 
                   "End Date",
                   "Allocated Capital",
                   "Average GMV ($)",
                   "Number of Instruments",
                   "Cumulative P&L ($)",
                   "Annualized P&L ($)",
                   "Annualized P&L Volatility ($)",
                   "Annualized Return on AC (%)",
                   "Annualized Volatility on AC (%)",
                   "Sharpe Ratio",
                   "Best Month ($)",
                   "Worst Month ($)"
        )
        
        datavalues <- c(as.character(stat.value$day1),
                        as.character(stat.value$day2),
                        comma_format(digits = 4)(capital.num),
                        comma_format(digits = 4)(stat.value$gmv.mean),
                        x.intermediate$instrument,
                        comma_format(digits = 4)(stat.value$pnl$pnl.cum),
                        comma_format(digits = 4)(stat.value$pnl$pnl.annualized),
                        comma_format(digits = 4)(stat.value$pnl$pnl.vol), 
                        comma_format(digits = 3)(stat.value$annualizedret),
                        comma_format(digits = 3)(stat.value$volret),
                        comma_format(digits = 3)(stat.value$pnl$pnl.sharpe),
                        paste(as.character(stat.value$performance$best.month),
                              " (", comma_format(digits = 4)(stat.value$performance$best.pnl),")",
                              sep = ""),
                        paste(as.character(stat.value$performance$worst.month),
                              " (",comma_format(digits = 4)(stat.value$performance$worst.pnl),")",
                              sep = "")
        )
        
        tab <- as.data.frame(cbind(items, datavalues))
        
        return(tab)
        
      }, align = "rlr",
      include.colnames = FALSE,
      include.rownames = FALSE)
      
      ## The summary table for three best performers. The interface returns the 
      ## summary of performers aabout the whole portfolio if individual
      ## instruments are selected, and returns the summary about a single sector
      ## if individual sectors are selected
      
      output$best3Performer <- renderTable(
        
        withProgress(message = "Rendering Best Performer Table", 
                     value = 0, {
                       x.intermediate <- intermediate()
                       
                       if (input$instrument == "Sector") {
                         tbl <- best_worst_three(x.intermediate$sect)$best.3
                       } else {
                         tbl <- best_worst_three(x.intermediate$x.temp)$best.3
                       }
                       
                       ## Format the numbers with commas and roundings
                       
                       tbl$pnl <- as.character(tbl$pnl)
                       for(i in 2:nrow(tbl)){
                         tbl$pnl[i] <- comma_format(digit = 4)(as.numeric(tbl$pnl[i]))
                       }
                       
                       return(tbl)
                     }), align = "rrr",
        include.colnames = FALSE,
        include.rownames = FALSE)
      
      ## The summary table for three worst performers
      
      output$worst3Performer <- renderTable(
        withProgress(message = "Redering Worst Performer Table",
                     value = 0, {
                       x.intermediate <- intermediate()
                       
                       ## use search rather than drop-down menu
                       
                       if (input$instrument == "Sector") {
                         tbl <- best_worst_three(x.intermediate$sect)$worst.3
                       } else {
                         tbl <- best_worst_three(x.intermediate$x.temp)$worst.3
                       }
                       
                       tbl$pnl <- as.character(tbl$pnl)
                       for(i in 2:nrow(tbl)){
                         tbl$pnl[i] <- comma_format(digit = 4)(as.numeric(tbl$pnl[i]))
                       }
                       
                       return(tbl)
                     }), align = "rrr",
        include.colnames = FALSE,
        include.rownames = FALSE)
      
      ## The summary table for top three drawdowns. If the number of drawdowns
      ## available is less than 3 (like a single drawdown dominates the whole
      ## time span), the table will return fewer drawdowns.
      
      output$drawdown <- renderTable(
        withProgress(message = "Rendering Worst Drawdowns",
                     value = 0, {
                       stat.value <- valuefun()
                       
                       tbl <- stat.value$pnl$pnl.drawdown
                       
                       tbl$pnl <- as.character(tbl$pnl)
                       for(i in 2:nrow(tbl)){
                         tbl$pnl[i] <- comma_format(digit = 4)(as.numeric(tbl$pnl[i]))
                       }
                       
                       colnames(tbl) <- c("Start Date","End Date", "P&L ($)")
                       return(tbl)
                     }), align = "rrrr", 
        include.colnames = FALSE,
        include.rownames = FALSE)
      
      ## Draw the first plot for cumulative P&L and daily P&L. The choice of
      ## plot depends on the user's selection. The user can switch between
      ## cumulative O&L and daily P&L by clicking the radio buttons below the
      ## plot. Only cumulative P&L among all the plots is drawn in a line, so
      ## we put it into another function
      
      output$plot1 <- renderDygraph(
        withProgress(message = "Rendering P&L Plots",
                     value = 0, {
                       x.intermediate <- intermediate()
                       
                       ## check if the filtered data is empty
                       
                       if((!is.data.frame(x.intermediate$x)) | dim(x.intermediate$x)[1] == 0){
                         stop("There does not exist any data that satisfies your selection. Please select a new combination of strategy + portfolio + instrument.")
                       }
                       
                       if(input$radio1 == "cumpnl") {
                         
                         p1 <- cumpnl_plot(x = x.intermediate$x)
                       } 
                       else {
                         p1 <- interactive_plot(x = x.intermediate$x, type = "pnl")
                       }
                       return(p1)
                     }))
      
      ## Draw the second plot for nmv, gmv and number of contracts
      
      output$plot2 <- renderDygraph(
        withProgress(message = "Rendering Market Value Plots",
                     value = 0, {
                       x.intermediate <- intermediate()
                       
                       ## check if the filtered data is empty
                       
                       if((!is.data.frame(x.intermediate$x)) | dim(x.intermediate$x)[1] == 0){
                         stop("There does not exist any data that satisfies your selection. Please select a new combination of strategy + portfolio + instrument.")
                       }
                       
                       p2 <- interactive_plot(x = x.intermediate$x, type = input$radio2)
                       
                       return(p2)
                     }))
    }
  )
}