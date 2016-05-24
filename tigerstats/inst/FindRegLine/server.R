library(shiny)
library(magrittr)
library(DT)

# implement dean attali's local storage suggestion:
outputDir <- "scores"
dir.create(outputDir, showWarnings = FALSE)

saveData <- function(data) {
  # Create a unique file name
  fileName <- sprintf("%s_%s.csv", as.integer(Sys.time()), digest::digest(data))
  # Write the file to the local system
  write.csv(
    x = data,
    file = file.path(outputDir, fileName), 
    row.names = FALSE, quote = TRUE
  )
}

loadData <- function() {
  # Note:  function returns NULL if no scores yet
  # Read all the files into a list
  files <- list.files(outputDir, full.names = TRUE)
  data <- lapply(files, read.csv, na.strings = c("","NA"),
                 stringsAsFactors = FALSE) 
  # Concatenate all data together into one data.frame
  data <- do.call(rbind, data)
  data
}


# bounds for intercept
lowa <- -5
higha <- 5
# bounds for slope
lowb <- -2
highb <- 2
# error sd
sigma <- 3
n <- 10 # number of points in in cloud
x <- 1:10 # x-values

# read in players records
scoreCount <- length(list.files(outputDir))
if ( scoreCount > 0 ) {
  leaders <- loadData()
  } else leaders <- data.frame()

# Define server logic for FindRegLine
function(input, output, session) {
  
  # users should start in different places
  set.seed(as.numeric(Sys.time()))
  
  #initiate status values
  rv <- reactiveValues(
    beginning = TRUE,
    playing = FALSE,
    reporting = FALSE,
    leaders = leaders
  )
  
  #set up:
  player_rank <- NULL
  ta <- round(runif(1, min = lowa, max = higha), 2)
  tb <- round(runif(1, min = lowb, max = highb), 2)
  y <- ta+tb*x+rnorm(n,mean=0,sd=sigma)
  #SS for the regression line:
  mod <- lm(y~x)
  ess <- sum((resid(mod))^2)
  #determine nice limits for plot (and a slider):
  reg.slope <- coef(mod)[2]
  reg.int <- coef(mod)[1]
  #Find range of y-intercepts of lines through
  #points on scatterplot having slope = reg.slope
  int.min <- min(y-reg.slope*x)
  int.max <- max(y-reg.slope*x)
  int.band <- (int.max-int.min)/2
  #Expand this range, and make sure it includes 0:
  int.mid <- (int.max+int.min)/2
  lowa.slider <- floor(min(c(int.mid-1.2*int.band,-1,min(y)-1)))
  higha.slider <- ceiling(max(c(int.mid+1.2*int.band,1,max(y)+1)))
  #plot limits reflect this range, too:
  ymin <- lowa.slider
  ymax <- higha.slider
  y.mean <- mean(y)
  your.y <- rep(y.mean,n)
  #SS for the line initially placed (a= mean(y),b=0):
  total.ss <- sum((y-mean(y))^2)
  your.ss <- total.ss
  turns <- 0
  close <- 100
  score <- close + turns

  # initiate a slider info:
  rvSlider <- reactiveValues(
    lowa.slider = lowa.slider,
    higha.slider = higha.slider,
    y.mean = y.mean
  )
  
  #initiate bslider info
  rvbSlider <- reactiveValues(
    value = 0
  )
  
  rvGraph <- reactiveValues(
    ymin = ymin,
    ymax = ymax
  )

  #make the a and b sliders
  output$aslider <- renderUI({
        sliderInput("a",min=rvSlider$lowa.slider,max=rvSlider$higha.slider,
                label="y-Intercept",step=0.01,value=rvSlider$y.mean)
        })
  
  output$bslider <- renderUI({
        input$reset
        sliderInput("b",min=2*lowb,max=2*highb,label="Slope",
                    step=0.01,value=0)
          })
  
  observeEvent(input$updateBoard,
               {
                  temp <- loadData()
                  if (! is.null(temp)) {
                    rv$leaders <- temp
                    } else rv$leaders <- data.frame()
               }
               )
  
  observeEvent(input$submit,
               {
                 rv$beginning <- FALSE
                 rv$playing <- TRUE
                 turns <<- turns + 1
                 your.y <<- input$a+input$b*x
                 your.ss <<- sum((y-your.y)^2)
                 close <<- 100*(your.ss-ess)/(total.ss-ess)
                 score <<- turns+close
               })
  
  observeEvent(input$enditall,
               {
                 rv$reporting <- TRUE
                 rv$playing <- FALSE
                 if (input$player != "") {
                   # make record for the game just ended
                   name <- input$player
                   lastScore <- score
                   time <- Sys.time()
                   # update the leader board in case others are playing.
                   # first, get the latest data:
                   leaders <- loadData()
                   # compute rank of player:
                   if ( is.null(leaders) ) rank <- 1
                   if ( ! is.null(leaders) ) {
                     # sort it:
                     leaders <- leaders[order(leaders$score),]
                     leScore <- max(which(lastScore >= leaders$score))
                     # the above will return -Inf if our player has best score
                     if (!is.infinite(leScore)) {
                      rank <- leScore + 1
                      } else {
                        rank <- 1
                      }
                     }
                   # store rank for reporting:
                   player_rank <<- rank
                   # add this game to the board
                   game <- data.frame(name = name, score = lastScore, 
                                      time = time)
                   saveData(game)
                   # update board so user will see his/her name right away:
                   rv$leaders <- loadData()
                 }
               })
  
  observeEvent(input$reset,
              {
                rv$beginning <- TRUE
                rv$reporting <- FALSE
                rv$playing <- FALSE
                #set up again:
                ta <<- round(runif(1, min = lowa, max = higha), 2)
                tb <<- round(runif(1, min = lowb, max = highb), 2)
                y <<- ta+tb*x+rnorm(n,mean=0,sd=sigma)
                #SS for the regression line:
                mod <<- lm(y~x)
                ess <<- sum((resid(mod))^2)
                #determine nice limits for plot (and a slider):
                reg.slope <<- coef(mod)[2]
                reg.int <<- coef(mod)[1]
                #Find range of y-intercepts of lines through
                #points on scatterplot having slope = reg.slope
                int.min <<- min(y-reg.slope*x)
                int.max <<- max(y-reg.slope*x)
                int.band <<- (int.max-int.min)/2
                #Expand this range, and make sure it includes 0:
                int.mid <<- (int.max+int.min)/2
                rvSlider$lowa.slider <- floor(min(c(int.mid-1.2*int.band,-1,min(y)-1)))
                rvSlider$higha.slider <- ceiling(max(c(int.mid+1.2*int.band,1,max(y)+1)))
                #plot limits reflect this range, too:
                rvGraph$ymin <- rvSlider$lowa.slider
                rvGraph$ymax <- rvSlider$higha.slider
                y.mean <<- mean(y)
                rvSlider$y.mean <- y.mean # so a slider sets correctly
                your.y <<- rep(y.mean,n)
                #SS for the line initially placed (a= mean(y),b=0):
                total.ss <<- sum((y-mean(y))^2)
                your.ss <<- total.ss
                turns <<- 0
                close <<- 100
                score <<- close + turns
              })
  
  output$beginning <- reactive({
    rv$beginning
  })
  
  output$playing <- reactive({
    rv$playing
  })
  
  output$reporting <- reactive({
    rv$reporting
  })
  
  outputOptions(output,"beginning", suspendWhenHidden = FALSE)
  outputOptions(output,"playing", suspendWhenHidden = FALSE)
  outputOptions(output,"reporting", suspendWhenHidden = FALSE)
  
  
 output$gamecloud <- renderPlot({
   input$submit #just in case user decides not to change line
   plot(x,y,pch=16,col="blue",ylim=c(rvGraph$ymin,rvGraph$ymax),
        xlim=c(0,n))
   points(0,0,cex=0.8,pch=16,col="green")
   abline(input$a,input$b)
   abline(0,0,lty=2,col="green")
   lines(x=c(0,0),y=c(rvGraph$ymin,rvGraph$ymax),lty=2,col="green")
   current.y <- input$a + input$b * x
   for(i in 1:n)  {
     lines(x=c(x[i],x[i]),y=c(current.y[i],y[i]))
   }
 })
 
 output$finalcloud <- renderPlot({
   input$enditall
   ymin <- rvGraph$ymin
   ymax <- rvGraph$ymax
   plot(x,y,pch=16,col="blue",ylim=c(ymin,ymax),
        xlim=c(0,n))
   points(0,0,cex=0.8,pch=16,col="green")
   isolate(abline(input$a,input$b))
   abline(0,0,lty=2,col="green")
   lines(x=c(0,0),y=c(ymin,ymax),lty=2,col="green")
   coefs <- coef(mod)
   abline(coefs,col="red",lwd=3)
 })
 
 output$score <- renderTable({
   input$submit
   input$reset
   tab <- rbind(your.ss,ess,close,turns,round(score,3))
   colnames(tab) <- "Report"
   rownames(tab) <- c("Your ESS",
              "Reg Line's ESS",
              "Closeness Measure",
              "Turns So Far","Score So Far")
   tab
 })
 
 output$rank <- reactive({
   input$enditall
   paste0("<h3>Your rank for this game is: ", player_rank,"</h3>")
 })

# before DT: 
#  output$leaders <- renderDataTable({
#    input$enditall
#    input$updateBoard
#    leaders <<- read.csv(file = "leaders.csv", 
#             header = TRUE, stringsAsFactors = FALSE)
#    leaders[order(leaders$score),]
#  })
#  
# outputOptions(output, "leaders", suspendWhenHidden = FALSE)
 
# try with DT
 
  observeEvent(input$updateBoard, {
    if ( ! rv$beginning ) {
      temp <- loadData()
      if ( ! is.null(temp) ) {
        rv$leaders <- temp[order(temp$score), ]
        }
      }
  })
  
  output$leaders <- DT::renderDataTable({
    leaders <- rv$leaders
    if ( nrow(leaders ) > 0 ) {
      leaders <- leaders[order(leaders$score), ]
      leaders$rank <- 1:nrow(leaders)
      datatable(leaders, rownames = FALSE, 
                caption = "Here is the Leader Board:",
                filter = "bottom") %>%
        formatStyle(
          'name',
          backgroundColor = styleEqual(input$player, 'lightblue'))
      }
  })
 
 output$revelation <- renderTable({
   if (input$enditall > 0) {
     coefs <- round(coef(mod),2)
     tab <- rbind(c(input$a,input$b),coefs)
     rownames(tab) <- c("Your Line","Regression Line")
     colnames(tab) <- c("y-Intercept","Slope")
     tab
     }
 })
 
 output$downloadData <- downloadHandler(
   filename = function() {
     "leaderboard.csv"
   },
   content = function(file) {
     sep <- ","
     # Write to a file specified by the 'file' argument
     write.table(leaders, file, sep = sep,
                 row.names = FALSE)
   }
 )
  
}

