  options(warn = -1)  # remove warnings 
	# create year list from 1994 to current year
	yr <- as.character(c(1994:as.numeric(format(Sys.Date(), "%Y"))))

	# message for choosing working directory
	tryCatch({
		setwd(jchoose.dir(default = getwd(), caption = "Choose working directory"))
	}, error = function(e) {
		msg0 <- paste0("All files will be downloaded into current working directory: ", getwd(), "\n You can change working directory from home page")
		err <- tcltk::tkmessageBox(message = msg0, icon = "info")
	})

	## import functions from edgar_shiny_functions.R
	source(system.file("shiny-edgar/edgarapp/edgar_shiny_functions.R", package = "edgar"), local = TRUE)

	# read positive and negative dictionary
	neg.words <- utils::read.csv(system.file("data/negwords.csv", package = "edgar"))
	pos.words <- utils::read.csv(system.file("data/poswords.csv", package = "edgar"))
	neg.words <- neg.words$WORDS
	pos.words <- pos.words$WORDS


shinyServer(function(input, output, session) {
    
    # Event observer for Change working directory button
    observeEvent(input$change.wd, {
        tryCatch({
            setwd(jchoose.dir(default = getwd(), caption = "Choose working directory"))
        }, error = function(e) {
            msg01 <- paste0("All files will be downloaded into current working directory: ", 
							getwd(), "\n You can change working directory from home page")
            err <- tcltk::tkmessageBox(message = msg01, icon = "info")
        })
    })
    
    # Event observer for Get Master Index Tab
    observeEvent(input$down.master, {
        if (input$getmasteryears > 0) {
            dir.create("Master Index")
            year.array <- unlist(strsplit(input$getmasteryears, " "))
            if ("ALL" %in% year.array) {
                year.array <- as.character(c(1994:as.numeric(format(Sys.Date(), "%Y"))))
            }
            
            output$master.download.status <- renderTable({
                withProgress(message = "Downloading Index:", value = 0.01, {
                  GetMasterIndexShiny(year.array)
                })
            }, escape = FALSE)
            
        }
    })
    
    # Event observer for Filing Info Tab
    observeEvent(input$getinfo, {
        if (input$year.getinfo > 0) {
            year <- input$year.getinfo
            if (file.exists(paste0("Master Index/", year, "master.Rda"))) {
                withProgress(message = "Please Wait", {
                  load(paste0("Master Index/", year, "master.Rda"))
                  data <- year.master
                })
                
                if (nrow(data) > 0) {
                  datasetInput <- reactive({
                    CreateLinkShiny(data)
                  })
                  
                  output$master.table <- renderDataTable({
                    datasetInput()
                  }, escape = FALSE)
                  
                  withProgress(message = "Please Wait", {
                    output$tb <- renderUI({
                      dataTableOutput("master.table")
                    })
                  })
                } else {
                  output$tb <- renderText({
                    ("Rda file is corrupted. Please re-download the master Index file for the selected year using 'Get Master Index' tab.")
                  })
                }
            } else {
                output$tb <- renderText({
                  ("Rda file not found in the current directory. Please download the master Index file for the selected year using 'Get Master Index' tab.")
                })
            }
        }
    })
    
    
    # Event handler for download Filing Tab
    output$stat.table.Filing1 <- renderDataTable({
        input$downFilings
        year3 <- isolate(input$year3)
        cik3 <- gsub("(^[[:space:]]+|[[:space:]]+$)", "", isolate(input$cik3))
        ftype3 <- gsub("(^[[:space:]]+|[[:space:]]+$)", "", isolate(input$ftype3))
        
        if (year3 > 0 & cik3 > 0 & ftype3 > 0) {
            withProgress(message = "Please Wait", {
                isolate(DownloadFilingsShiny(year3, cik3, ftype3))
            })
        }
    }, escape = FALSE)
    
    output$stat.table.Filing <- renderUI({
        dataTableOutput("stat.table.Filing1")
    })
    
    # Event handler for daily info tab
    observeEvent(input$get.dailyinfo, {
        withProgress(message = "Please Wait", {
            
            if (input$dailydate > 0) {
                dir.create("Daily Index")
                if (as.numeric(Sys.Date()) > input$dailydate) {
                  date <- as.character(input$dailydate)
                  date <- unlist(strsplit(date, "-"))
                  year <- date[1]
                  month <- date[2]
                  day <- date[3]
                  
                  data <- GetDailyInfoShiny(day, month, year)
                  
                  if (data != 0) {
                    output$daily.table <- renderDataTable({
                      withProgress(message = "Please Wait", {
                        data
                      })
                    }, escape = FALSE)
                    
                    output$daily.tb <- renderUI({
                      dataTableOutput("daily.table")
                    })
                  } else {
                    output$daily.table <- renderText({
                      "Server error or daily index not available for this date"
                    })
                    
                    output$daily.tb <- renderUI({
                      textOutput("daily.table")
                    })
                  }
                } else {
                  output$daily.table <- renderText({
                    "Please select the appropriate date."
                  })
                  
                  output$daily.tb <- renderUI({
                    textOutput("daily.table")
                  })
                }
            }
        })
    })
    
    # Event observer for sentiment analysis tab
    observeEvent(input$get.sentianalysis, {
        withProgress(message = "Please Wait", {
            filpath <- jchoose.files(default = getwd(), caption = "Select 10-K file", multi = FALSE)
            
            if (grepl("10-K", filpath) && grepl(".txt", filpath)) {
                word_frq <- GetwordfrqShiny(filpath)
                words <- unlist(word_frq$WORD)
                neg.word.table <- word_frq[words %in% neg.words, ]
                pos.word.table <- word_frq[words %in% pos.words, ]
                
                polaritydata <- reactive({
                  data.frame(Polarity = c("Negative", "Positive"), 
                             Total_words = c(sum(neg.word.table$FREQUENCY), sum(pos.word.table$FREQUENCY)))  
                })
                
                
                
                output$negwordcloud <- renderPlot({
				  # Wordcloud for nagative words
                  wordcloud::wordcloud(words = neg.word.table$WORD, freq = neg.word.table$FREQUENCY, 
									   scale = c(4, 0.8), max.words = Inf, 
									   random.order = F, colors = RColorBrewer::brewer.pal(8, "Dark2"))
                })
                
                output$neghist <- renderPlot({
				    # Extracting top 15 negative words
					neghistdata <- neg.word.table[1:15, ]
					# Negative words histogram
					ggplot2::qplot(neghistdata$WORD, weight = neghistdata$FREQUENCY, data = neghistdata, geom = "bar", xlab = "Word", ylab = "Frequency") + 
					ggplot2::coord_flip() + 
					ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1, vjust = 0.5)) + 
					ggplot2::geom_bar(fill = "#FF0000") + 
					ggplot2::theme(axis.text = ggplot2::element_text(size = 11, face = "bold"), axis.title = ggplot2::element_text(size = 13, face = "bold")) + 
					ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = "black")) + 
					ggplot2::theme(plot.title = ggplot2::element_text(lineheight = 3, face = "bold", color = "black", size = 15))
                })
                
                output$poswordcloud <- renderPlot({
				   # Wordcloud for positive words
					wordcloud::wordcloud(words = pos.word.table$WORD, freq = pos.word.table$FREQUENCY, 
				                       scale = c(4, 0.8), max.words = Inf, 
                                       random.order = F, colors = RColorBrewer::brewer.pal(8, "Dark2"))
                })
                
                output$poshist <- renderPlot({
				    # Extracting top 15 positive words
					poshistdata <- pos.word.table[1:15, ]
					# Positive words histogram
					ggplot2::qplot(poshistdata$WORD, weight = poshistdata$FREQUENCY, data = poshistdata, geom = "bar", xlab = "Word", ylab = "Frequency") + 
					ggplot2::coord_flip() + 
					ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1, vjust = 0.5)) + 
					ggplot2::geom_bar(fill = "#339900") + 
					ggplot2::theme(axis.text = ggplot2::element_text(size = 11, face = "bold"), axis.title = ggplot2::element_text(size = 13, face = "bold")) + 
					ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = "black")) + 
					ggplot2::theme(plot.title = ggplot2::element_text(lineheight = 3, face = "bold", color = "black", size = 15))
                })
				
				        
                output$polarity.hist <- renderPlot({
					# Polarity histogram
					ggplot2::ggplot(data = polaritydata(), ggplot2::aes(x = Polarity, y = Total_words)) + 
					ggplot2::geom_bar(fill = c("#FF0000", "#339900"), stat = "identity", width = 0.4, position = ggplot2::position_dodge(width = 0.5)) + 
					ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 0.5)) + 
					ggplot2::theme(axis.text = ggplot2::element_text(size = 14, face = "bold"), axis.title = ggplot2::element_text(size = 16, face = "bold")) + 
					ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = "black")) + 
					ggplot2::xlab("Polarity") + 
					ggplot2::ylab("Frequency") + 
					ggplot2::theme(plot.title = ggplot2::element_text(lineheight = 3, face = "bold", color = "black", size = 18))
                })
                
                output$negwordcloud1 <- renderUI({
                  plotOutput("negwordcloud")
                })
                output$neghist1 <- renderUI({
                  plotOutput("neghist")
                })
                output$poswordcloud1 <- renderUI({
                  plotOutput("poswordcloud")
                })
                
                output$poshist1 <- renderUI({
                  plotOutput("poshist")
                })
                output$polarity.hist1 <- renderUI({
                  plotOutput("polarity.hist")
                })
            } else {
                output$senti.error <- renderText({
                  "Error: Please select 10-K Filing"
                })
                
                output$negwordcloud1 <- renderUI({
                  NULL
                })
                output$neghist1 <- renderUI({
                  NULL
                })
                
                output$poswordcloud1 <- renderUI({
                  NULL
                })
                output$poshist1 <- renderUI({
                  NULL
                })
                
                output$polarity.hist1 <- renderUI({
                  NULL
                })
                
                msg33 <- "Please select 10-K file only.."
                err <- tcltk::tkmessageBox(message = msg33, icon = "error")
            }
        })
    })
}) 
