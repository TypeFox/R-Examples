library(shiny)
shinyServer(function(input, output, session) {

  # ulfs mail
  # 1 need aus mail: chisq test statistic
  # 2 odds ratios
  
  # reactive dataset object ####
  # this object is basically a list containing
  # - infile
  # - names
  # - cleaned file for display
  dataset <- reactive({
    if(is.null(input$file1)) return(NULL)
    upfile <- read.csv2(input$file1$datapath,
                        sep = input$sep,
                        dec = input$dec,
                        quote = input$quote,
                        header = input$header)
    # wow this is super super ugly, 
    # took me three hours to fix. 
    # NA are not nice for render output,
    # but my function does not work without NA.
    # might need an NA string for my read in... 
    no_na <- upfile
    no_na[is.na(no_na)] = ''
    
    react_list <- list()
    initial_names <- names(upfile)
    react_list$nms <- initial_names
    react_list$in_file <- upfile
    react_list$clean_file <- no_na
    react_list
  })


  # intermediate reactive object to compute results ####
  # select the relevant columns
  # get the number of questions
  # split it by experimental condition
  # compute contigency tables for all points in time
  # this is only needed for Chisq stuff 
  intermediate <- reactive({
    if(is.null(input$file1)) return(NULL)
    intermed <- dataset()$in_file[,c(input$quest_cols,input$cond_col)]
    intermed$drop_out <- extract_drop_out_from_df(intermed,input$quest_cols)
    
    le <- length(input$quest_cols) 
    
    sp <- split(intermed,factor(intermed[,input$cond_col]))
    
    out_li <- lapply(sp,function(x){
      
      grd <- expand.grid(1:le,0)
      # count drop out positions, like table() 
      pos_count <- plyr::count(x$drop_out)
      names(grd) <- names(pos_count)
      grd[pos_count$x,'freq'] <- pos_count$freq
      
      cs <- cumsum(grd$freq)
      cs[-length(cs)]
      
    })
    
    participants <- table(intermed[,input$cond_col])
    participant_stats <- list()
    participant_stats$participants <- participants
    participant_stats$out_li <- out_li
    participant_stats$dropout <- 
      unlist(lapply(out_li,'[',input$select_question))
    participant_stats$remain <- 
      participants - participant_stats$dropout

    participant_stats
  })

  # rendered UI components #### 
  output$choose_questions <- renderUI({
    if(is.null(dataset()$nms)) return(NULL)
#     checkboxGroupInput("quest_cols", "Choose question columns", 
#                        choices  = dataset()$nms,
#                        selected = dataset()$nms[-c(1,length(dataset()$nms))])
    selectInput('quest_cols', 'Choose question columns',
                selected = dataset()$nms[-c(1,length(dataset()$nms))],
                multiple=TRUE,
                choices  = dataset()$nms,
                selectize=FALSE)
  })
  
  output$choose_condition <- renderUI({    
    if(is.null(dataset()$nms)) return(NULL)
    selectInput('cond_col',
                'Choose the column that holds the experimental conditions',
                choices = c('None',dataset()$nms),
                selected = 'None',
                multiple = F
                )
  })
  
  
 output$choose_question <- renderUI({    
   if(is.null(dataset()$nms)) return(NULL)
   sliderInput('select_question',
              'Select Question to be tested',
              min = 1,
              max = length(input$quest_cols)-1,
              step = 1,
              value = 1
  )
})



# dynamically rendertable , ifelse did not work for some reason, 
# however now the div class wrapping the table depends on the 
# p.value... thus we can change the color of the table background 
# depending on the p-value
output$color_table <- renderUI({  
  pval <- chisq.test(as.table(cbind(intermediate()$remain,
                                   intermediate()$dropout)))$p.value
  
  # this is a hack just to have no error message when there is no p-value
  pval <- ifelse(is.na(pval),1,pval)
  
  
  if(pval < .05){
    tags$div(class='alertme',tableOutput('chisq_table'))
  }  else {
    tags$div(class='noalert',tableOutput('chisq_table'))
  }
  
})

# end of dynamic UI components ####
output$out_table <- renderTable({
    head(dataset()$clean_file)
  })
  
# Chisq results for the selected question  
output$chisq_table <- renderTable({
    data.frame(question = as.character(input$select_question),
               p.value =  chisq.test(as.table(cbind(intermediate()$remain,
                                                    intermediate()$dropout)))$p.value
                )
  },include.rownames=F)
  

# The contigency table ####
output$contingency <- renderTable({
    share_remain = data.frame(intermediate()$remain)$Freq / data.frame(intermediate()$participants)$Freq
    
    data.frame(
      condition = data.frame(intermediate()$participants)$Var1,
      N = data.frame(intermediate()$participants)$Freq,
      dropout = as.character(intermediate()$dropout),
      remain = as.character(data.frame(intermediate()$remain)$Freq)
      )
  },include.rownames = F)

# create an odds ratio table
output$or_table <- renderTable({
  share_remain = data.frame(intermediate()$remain)$Freq / data.frame(intermediate()$participants)$Freq
  
  odds <- get_odds(share_remain)
  condition <- data.frame(intermediate()$participants)$Var1
  
  m <- combn(odds,2)
  or <- m[1,]/m[2,]
  nms <- combn(as.character(condition),
               2,paste0,collapse=' vs. ')
  
  d <- data.frame(t(or))
  colnames(d) <- nms
  d
  
#   data.frame(t(data.frame(
#     condition = nms,
#     OR = round(or,digits = 2),row.names = NULL
#     )))
  
},include.rownames = F)




# Create the Chart  
output$mychart <- renderLineChart({
    if(is.null(input$file1)) return(NULL)
    if(is.null(input$quest_cols)) return(NULL)
    if(is.null(input$cond_col)) return(NULL)
#     
    intermed <- dataset()$in_file[,c(input$quest_cols,input$cond_col)]
      
    intermed$drop_out <- extract_drop_out_from_df(intermed,input$quest_cols)
    le <- length(input$quest_cols)
    
    xl <- lapply(split(intermed,factor(intermed[,input$cond_col])),
                 function(x) compute_shares_remain(x,x$drop_out,le))
    df <- as.data.frame(xl)
    df$total <- compute_shares_remain(intermed,intermed$drop_out,le)
    names(df) <- gsub('X','',names(df))
    df
        
    })
})



# my little debugger... 
#   output$test <- renderText({
#     intermed <- dataset()$in_file[,c(input$quest_cols,input$cond_col)]
#     intermed$drop_out <- extract_drop_out_from_df(intermed,input$quest_cols)
#     intermed$drop_out
#     
#   })



