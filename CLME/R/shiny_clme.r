#' Shiny GUI for CLME
#'
#' @description Opens a graphical user interface to run \pkg{CLME}, built from the \pkg{shiny} package.
#'
#' @rdname shiny_clme
#'
#' @param input  input from GUI.
#' @param output output to GUI.
#' 
#' @details 
#' Currently the GUI does not allow specification of custom orders for the alternative hypothesis. Future versions may enable this capability.
#' The data should be a CSV or table-delimited file with the first row being a header. Variables are identified using their column letter or number (e.g., 1 or A). Separate multiple variables with a comma (e.g., 1,2,4 or A,B,D), or select a range of variables with a dash (e.g., 1-4 or A-D). Set to 'None' (default) to indicate no covariates or random effects.
#' If group levels for the constrained effect are character, they may not be read in the proper order. An extra column may contain the ordered group levels (it may therefore have different length than the rest of the dataset).
#' 
#' @note
#' This function is primarily designed to be called by \code{\link{clme}}. 
#' 
#' By default, homogeneous variances are assumed for the residuals and (if included) random effects. Heterogeneity can be induced using the arguments \code{Nks} and \code{Qs}, which refer to the vectors \eqn{ (n_{1}, n_{2}, \ldots, n_{k}) }{(n1, n2 ,... , nk)} and \eqn{ (c_{1}, c_{2}, \ldots, c_{q}) }{(c1, c2 ,... , cq)}, respectively. See \code{\link{clme_em}} for further explanation of these values.
#' 
#' 
#' @examples
#' \dontrun{ shiny_clme() }
#' 
#' @import shiny
#' @export
#' 
shiny_clme <- function(){
  library("CLME")
  runApp(
    list(
      ui     = shinyUI_clme,
      server = shinyServer_clme
    )
  )
}




##############################################################################
##
## The user interface for the shiny app
##
##############################################################################
#  shinyUI(bootstrapPage())

#'
#' @rdname shiny_clme
#' 
#' @export
#' 

shinyUI_clme <- fluidPage( 
  titlePanel("Constrained Linear Mixed Effects"),
  sidebarLayout(
    sidebarPanel(
      
      ##
      ## Data input
      ##
      fileInput('file1', 'Data (choose CSV file)',
                accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv', '.xlsx')
      ),
      selectInput(inputId = "dlmt",
                  label = "Delimiter:",
                  choices=c("Comma-delimited", "Tab-delimited", "xlsx") 
      ), 
      ##
      ## Main controls
      ##
      hr(),
      selectInput(inputId = "solver",
                  label = "Type of Solver / fitting norm:",
                  choices=c("Least Squares (LS)", "Least Absolute Value (L1)",
                            "General LS (GLS)", "Asymmetrix LS" , 
                            "L1 approx", "Huber", "SILF", "Chebyshev",
                            "L-p Power Norm", "Quantile", "Poisson" ) 
      ),      
      selectInput(inputId = "tsfunc",
                  label = "Test Statistic:",
                  choices=c("LRT", "Williams") 
      ),
      checkboxGroupInput(inputId = "order",
                         label = "Order (select at least one):",
                         choices=c("Simple","Umbrella","Tree")
      ),
      checkboxInput(inputId = "decreasing",
                    label = "Decreasing (inverted) order:",
                    value=FALSE
      ),
      helpText("Identify columns of data"),
      helpText("Use column letters or numbers"),
      helpText("e.g., 1-3 or A-C or a list: 1,4,6"),
      textInput(inputId = "yy",
                   label = "Column of response:"),
      textInput(inputId = "p1",
                   label = "Column of constrained effect:"),
      textInput(inputId = "p2",
                   label = "Column(s) of Covariates:", value="None"),
      textInput(inputId = "q",
                   label = "Column(s) of random effects:", value="None"),

      numericInput(inputId = "nsim",
                   label = "Number of Bootstraps:",
                   min=0 , max=50000 , value=1000
      ),
      
      ##
      ## Action buttons
      ##
      hr(),
      helpText("Click to run model or update output."),
      actionButton(inputId = "compute1",
                   label = "Run model"
      ),
      actionButton(inputId = "compute2",
                   label = "Update plot"
      ),
      actionButton(inputId = "compute3",
                   label = "Update summary"
      ),
      
      
      ##
      ## Output controls
      ##
      hr(),
      checkboxInput(inputId = "outcheck",
                    label = "Format Output:",
                    value=FALSE
      ),
      conditionalPanel(condition = "input.outcheck",
                       checkboxInput(inputId = "plotci",
                                     label = "CI on Plot:",
                                     value=FALSE
                       )
      ),
      conditionalPanel(condition = "input.outcheck",
                       sliderInput(inputId = "alpha",
                                   label = "Alpha level:",
                                   min=0 , max=0.15 , value=0.05 , step=0.01
                       )
      ),
      conditionalPanel(condition = "input.outcheck",
                       checkboxInput(inputId = "makeFactor",
                                   label = "Force constrained effect to be factor",
                                   value=FALSE
                       )
      ),
      #br(),
      #conditionalPanel(condition = "input.outcheck",
      #                 helpText("Number of decimal places for:")
      #),
      #conditionalPanel(condition = "input.outcheck",
      #                 sliderInput(inputId = "digits",
      #                             label = "p-values:",
      #                             min=0 , max=8 , value=4 , step=1
      #                 )
      #),

      ##
      ## Extra parameters
      ##
      hr(),
      checkboxInput(inputId = "varssq",
                    label = "Heteroscedasticity:",
                    value=FALSE
      ),
      conditionalPanel(condition = "input.varssq",
                       textInput(inputId = "gfix",
                                    label = "Column of variance groups:")
      ),
      checkboxInput(inputId = "xlevel1",
                    label = "Define order of constrained groups:",
                    value=FALSE
      ),
      conditionalPanel(condition = "input.xlevel1",
                       textInput(inputId = "xlevels",
                                 label = "Column of ordered group levels:")
      ),
      ##
      ## Technical controls
      ##
      checkboxInput(inputId = "technical",
                    label = "Select Control Parameters:",
                    value=FALSE
      ),
      conditionalPanel(condition = "input.technical",
                       numericInput(inputId = "emiter",
                                    label = "Max EM Iterations:",
                                    min=10 , max=50000 , value=500
                       )
      ),
      conditionalPanel(condition = "input.technical",
                       numericInput(inputId = "mqiter",
                                    label = "Max MINQUE Iterations:",
                                    min=10 , max=50000 , value=500)
      ),
      conditionalPanel(condition = "input.technical",
                       numericInput(inputId = "emeps",
                                    label = "EM Convergence Criteria:",
                                    min=0 , max=50000 , value=0.0001)
      ),
      conditionalPanel(condition = "input.technical",
                       numericInput(inputId = "mqeps",
                                    label = "MINQUE Convergence Criteria:",
                                    min=0 , max=50000 , value=0.0001)
      ),
      conditionalPanel(condition = "input.technical",
                       numericInput(inputId = "ranseed",
                                    label = "Set RNG Seed:",
                                    min=0 , max=Inf , value=42)
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data Summary" , 
                 plotOutput(outputId = "fig0", height = "650px"),
                 tableOutput(outputId = "sum_table")
        ),
        tabPanel("Model Summary", 
                 verbatimTextOutput(outputId = "summary")
        ), 
        tabPanel("Model Plot"   ,
                 plotOutput(outputId = "fig1", height = "650px")
        ),
        tabPanel("Model Data"   ,
                 dataTableOutput(outputId = "datatbl")
        )
      )
    )
  ))





##############################################################################
##
## The server for the shiny app
##
##############################################################################

#'
#' @rdname shiny_clme
#' 
#' 
#' @import graphics
#' @importFrom openxlsx read.xlsx
#' @importFrom utils read.csv
#' 
#' @export
#' 
shinyServer_clme <- function(input, output) {
    
  clme_out <- reactive({
    
    compute1 <- input$compute1
    
    ## Put all the code to run CLME inside this
    if( compute1 > 0 ){
      isolate({
        
        file1   <- input$file1[4]
        solver  <- input$solver
        dlmt    <- input$dlmt
        #data1   <- as.matrix( read.csv( file=paste(file1) ) )
        
        if( dlmt=="Comma-delimited"){
          data1   <- read.csv( file=paste(file1) )
        }
        if( dlmt=="Tab-delimited" ){
          data1   <- read.csv( file=paste(file1) , sep="\t")
        }
        if( dlmt=="xlsx"){
          data1   <- read.xlsx( xlsxFile=paste(file1), colNames=TRUE )
        }
                
        yy      <- input$yy
        p1      <- input$p1
        p2      <- input$p2
        q       <- input$q
        nsim    <- input$nsim
        
        alpha      <- 0.05        
        makeFactor <- FALSE
        if( input$outcheck==TRUE ){
          makeFactor <- input$makeFactor
          alpha      <- input$alpha
        }
        
        if( p2 == '' | p2==0 ){ p2 <- "None" }
        if(  q == '' |  q==0 ){  q <- "None" }
        
        ##
        mySolver <- switch( EXPR=solver ,
                            "Least Squares (LS)"        = "LS",
                            "Least Absolute Value (L1)" = "L1",
                            "General LS (GLS)"          = "GLS",
                            "Asymmetrix LS"             = "asyLS",
                            "L1 approx"                 = "L1eps",
                            "Huber"                     = "huber",
                            "SILF"                      = "SILF",
                            "Chebyshev"                 = "chebyshev",
                            "L-p Power Norm"            = "Lp",
                            "Quantile"                  = "quantile",
                            "Poisson"                   = "poisson")       
        
                
        ## Get the indexes for X2 and U
        parse_idx <- function( arg1 ){
          arg1 <- toupper(arg1)
          if( arg1=="NONE" ){
            p4 <- 0
          } else{
            p2s <- strsplit(arg1, ",")[[1]]
            if( length(p2s) > 1 ){
              p3 <- vector()
              for( ii in 1:length(p2s) ){
                next1 <- sapply( p2s[ii],  FUN=function(x){ strsplit(x, "-") }  )[[1]]
                
                if( length(next1)==1 ){
                  if(  next1 %in% LETTERS ){
                    p3 <- append( p3 , which(next1 == LETTERS) )
                  } else{
                    p3 <- append( p3 , as.numeric(next1) )
                    
                  }
                } else{
                  next1a <- next1[1]
                  next1b <- next1[2]
                  if( next1a%in%LETTERS & next1b%in%LETTERS ){
                    n1a <- which(next1a == LETTERS)
                    n1b <- which(next1b == LETTERS)
                    p3 <- append( p3 , n1a:n1b )
                  } else{
                    p3 <- append( p3 , as.numeric(next1a):as.numeric(next1b) )
                  }
                }
              }
            } else{
              next1 <- sapply( p2s,  FUN=function(x){ strsplit(x, "-") }  )[[1]]
              
              if( length(next1)==1 ){
                if(  next1 %in% LETTERS ){
                  p3 <-  which(next1 == LETTERS)
                } else{
                  p3 <-as.numeric(next1)
                }
              } else{
                next1a <- next1[1]
                next1b <- next1[2]
                if( next1a%in%LETTERS & next1b%in%LETTERS ){
                  n1a <- which(next1a == LETTERS)
                  n1b <- which(next1b == LETTERS)
                  p3  <- n1a:n1b
                } else{
                  p3 <- as.numeric(next1a):as.numeric(next1b)
                }
              }
            }
              p4 <- as.numeric(p3)
          }
          p4
        }
        
        idx_yy <- parse_idx(yy)
        yn <- colnames(data1)[idx_yy]
        
        idx_x1 <- parse_idx(p1)
        x1n <- colnames(data1)[idx_x1]
        
        cov <- ran <- FALSE
        if( p2 != "None" ){
          idx_x2 <- parse_idx(p2)
          x2n <- colnames(data1)[idx_x2]
          cov <- TRUE
        }
        if( q != "None" ){
          idx_u  <- parse_idx(q)
          uun <- colnames(data1)[idx_u]
          ran <- TRUE
        }
        
        ncon   <- length(idx_x1)
        
        
        ## Create the formula
        
        ## Reorder the levels of X1 if need be
        if( ncon==1 & input$xlevel1 ){
          xlev    <- parse_idx(input$xlevels)
          nlev    <- length(levels(as.factor( data1[,xlev] )))
          xlevels <- as.character( data1[1:nlev,xlev] )
          if( any(xlevels=="") ){
            xlevels <- xlevels[ -which(xlevels=="")]  
          }
          data1[,idx_x1] <- factor( data1[,idx_x1], levels=xlevels )
        } else if( ncon==1 & makeFactor==TRUE ){
          data1[,idx_x1] <- factor( data1[,idx_x1] )
        }
        
        ## Build the formula
        x1f <- paste( x1n , collapse=" + ")
        
        if( !cov & !ran ){
          ## No extra terms
          formula <- formula( paste( yn, "~", x1f ) )
        } else if( cov & !ran ){
          ## Yes covariates, no random effects
          x2f     <- paste( x2n , collapse=" + ")
          formula <- formula( paste( yn, "~", x1f, "+", x2f ) )
        } else if( !cov & ran ){
          ## No covariates, yes random effects
          uf      <-  paste( "(1|", uun , ")", collapse=" + " )
          formula <- formula( paste( yn, "~", x1f, "+", uf ) )
        } else if( cov & ran ){
          ## Covariates and random effects
          x2f     <- paste( x2n , collapse=" + ")
          uf      <-  paste( "(1|", uun , ")", collapse=" + " )
          formula <- formula( paste( yn, "~", x1f, "+", x2f, "+", uf ) )
        }
                
        
        ## Input control arguments
        tsf <- lrt.stat
        if( input$tsfunc=="Williams" ){
          tsf <- w.stat
        }
        
        constraints <- list()
        constraints$order      <- tolower(input$order)
        constraints$decreasing <- input$decreasing
        
        rep.idx <- grep( "tree" , constraints$order )
        constraints$order[rep.idx] <- "simple.tree"
        
        gfix <- NULL
        if( input$varssq ){
          idx_grp <- parse_idx(input$gfix)
          gfix <- data1[,idx_grp]
        }
        
        emiter <- 500
        mqiter <- 500
        emeps  <- 0.00001
        mqeps  <- 0.00001
        seedvl <- NULL
                
        if( input$technical==TRUE ){
          emiter <- input$emiter
          mqiter <- input$mqiter
          emeps  <- input$emeps
          mqeps  <- input$mqeps
          seedvl <- input$ranseed
        }
        
        ## Run the model
        data2 <- as.data.frame( data1 )
        
        withProgress(message = 'Status:', value = 1, {
          setProgress(1/2, detail = paste("Computing"))
          
          if( ncon==1 & input$xlevel1 ){
            clme.results <- summary( clme(
              formula=formula, data=data2, gfix=gfix, constraints=constraints, 
              verbose = c(FALSE, FALSE, FALSE), tsf=tsf, tsf.ind = w.stat.ind,
              mySolver=mySolver, ncon=ncon, levels=list(idx_x1, xlevels), 
              em.eps=emeps, em.iter=emiter, mq.eps=mqeps, mq.iter=mqiter
            ), alpha=alpha, seed=seedvl, nsim=nsim )
          } else{
            clme.results <- summary( clme(
              formula=formula, data=data2, gfix=gfix, constraints=constraints, 
              verbose = c(FALSE, FALSE, FALSE), tsf=tsf, tsf.ind = w.stat.ind,
              mySolver=mySolver, ncon=ncon,
              em.eps=emeps, em.iter=emiter, mq.eps=mqeps, mq.iter=mqiter
            ), alpha=alpha, seed=seedvl, nsim=nsim )
          }
          
        })
        
        clme.results
        
      })
      
      
    }   
  })
  
  ##
  ## Boxplot of the data
  ##
  output$fig0 <- renderPlot({
    clme_out <- clme_out()
    if( length(clme_out)>1 ){
      dframe   <- clme_out$dframe
      if( is.factor( dframe[,2] ) ){
        if( length(levels(dframe[,2]))==clme_out$P1 ){        
          boxplot( dframe[,1] ~ dframe[,2],
                   xlab=colnames(dframe)[2], ylab=colnames(dframe)[1]  )
        }
        #if( (ncol(dframe)-1) >= clme_out$P1 ){
        #  xx1 <- apply( dframe[,2:(clme_out$P1+1)], 2, FUN=function(x){ (max(x)==1)*(min(x)==0) }  )
        #  
        #}
      } else{
        xx <- yy <- c(0,1)  
        plot( xx, yy, xlab="", ylab="", xaxt='n', yaxt='n', main="", col=0 )
        legend( "top", inset=0.45, box.lwd=0,
                legend="Cannot detect appropriate plot type.\n Try making constrained variable a factor.")
      }
    } else{
      plot( 1:5 , 1:5 , col="white", xaxt='n', yaxt='n', xlab="", ylab="", frame.plot=FALSE )
    }
  })
  
  
  output$sum_table <- renderTable({
    clme_out <- clme_out()
    
    if( length(clme_out)>1 ){
      dframe   <- clme_out$dframe
      
      funq1 <- function(x) quantile(x, 0.25 )
      funq3 <- function(x) quantile(x, 0.75 )
      
      xbar <- aggregate( dframe[,1] , by=list(dframe[,2]), FUN="mean")
      ndat <- aggregate( dframe[,1] , by=list(dframe[,2]), FUN="length")
      stdv <- aggregate( dframe[,1] , by=list(dframe[,2]), FUN="sd")
      minx <- aggregate( dframe[,1] , by=list(dframe[,2]), FUN="min")
      maxx <- aggregate( dframe[,1] , by=list(dframe[,2]), FUN="max")
      medn <- aggregate( dframe[,1] , by=list(dframe[,2]), FUN="median")
      qrt1 <- aggregate( dframe[,1] , by=list(dframe[,2]), FUN=funq1)
      qrt3 <- aggregate( dframe[,1] , by=list(dframe[,2]), FUN=funq3)
      
      tbl1 <- cbind( ndat, xbar[,2], stdv[,2], minx[,2], qrt1[,2], medn[,2], qrt3[,2], maxx[,2] )
      colnames(tbl1) <- c("Groups", "N", "Mean", "Std", "Min", "Q1", "Med", "Q3", "Max")
      
      format(tbl1, digits=3)
      
    } else{
      tbl1 <- matrix( "No results yet" , ncol=1, nrow=1)
      format(tbl1, digits=3)      
    }

  })
  
  
  ##
  ## Summarize the model
  ##
  output$fig1 <- renderPlot({
    
    clme_out <- clme_out()
    if( length(clme_out)>1 ){
      compute2 <- input$compute2
      if( compute2 > 0 ){
        
        isolate({
          ciwd  <- FALSE
          alpha <- 0.05
          if( input$outcheck==TRUE ){
            ciwd  <- input$plotci
            alpha <- input$alpha
          }
          plot( clme_out , ci=ciwd , alpha=alpha)
        })
      }      
    } else{
      plot( 1:5 , 1:5 , col="white", xaxt='n', yaxt='n', xlab="", ylab="", frame.plot=FALSE )
    }

  })
  
  output$summary <- renderPrint({
    
    clme_out <- clme_out()
    if( length(clme_out)>1 ){
      compute3 <- input$compute3
      if( compute3 > 0 ){
        isolate({
          clme_out
        })
      }
    } else{
      print( "Model has not yet been run."  )
    }
    
  })
  

  output$datatbl <-  renderDataTable({
    clme_out <- clme_out()
    if( length(clme_out) > 1 ){
      clme_out$dframe
    }
  })
  
  
  ## To check the output of various things when I'm updating
  #output$etc <-  renderPrint({
  #  print( clme_out() )
  #})
  
  
}








