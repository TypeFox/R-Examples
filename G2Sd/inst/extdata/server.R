library(shiny)
library(G2Sd)
library(xlsx)



grancompat <- function(x)
{
  x <- as.data.frame(x)
  n.sieve <- nrow(x)
  n.sample <- ncol(x)
  if (!is.data.frame(x)) 
    stop("dataframe expected.")
  if (any(x < 0))
    stop("negative entries in dataframe.")
  if (any(x > 300))
    warning("Some high values are present.", call. = FALSE,immediate.=TRUE)
  if (n.sieve!=29)
  {
    cat("Compatibility progress.... \n \n")
    
    ref_sieve=c(25000,20000,16000,12500,10000,8000,6300,5000,4000,2500,
                2000,1600,1250,1000,800,630,500,400,315,250,200,160,125,
                100,80,63,50,40,0)      
    
    init_df <- as.data.frame(matrix(data=0,ncol=n.sample,nrow=length(ref_sieve)));colnames(init_df) <- colnames(x)
    row.names(init_df) <-ref_sieve
    
    #     column.integration <- .mgrep(row.names(x),as.character(ref_sieve),FUN="which") ## error
    if (any(is.na(pmatch(row.names(x),ref_sieve))))
      stop("Incorrect sieve values.")
    
    else 
      #       init_df[column.integration,] <- x  ##error
    {for (sieve in row.names(x))
      init_df[sieve,] <- x[sieve,]}
    
  }
  else init_df <- x
  return(init_df)
}


shinyServer(function(input, output) {
  output$contents <- renderTable({
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    bdd_gran <- read.xlsx(inFile$datapath, sheetIndex=input$sheetindex, Header=input$header)
    row.names(bdd_gran) <- bdd_gran[,1];bdd_gran <- bdd_gran[,-1]
#     ; bdd_gran <- bdd_gran[,-dim(bdd_gran)[2]]
    bdd_gran <- round(bdd_gran,3);bdd_gran <- grancompat(bdd_gran)
    
    output$stat <- renderTable({ 
      
      if (input$gran=='Statistics') {
        result <-granstat(bdd_gran,statistic=input$method,aggr=FALSE)$stat
        
        output$downloadData <- downloadHandler(
          filename = "output.xlsx",
          content = function(file){
            fname <-paste0(input$filename,".xlsx") 
            wb <- createWorkbook()
            sheet  <- createSheet(wb, sheetName="Statistics")
            addDataFrame(result, sheet)
            saveWorkbook(wb, fname)
            #           
          })
      }
      if (input$gran=='Index') {
        result <-granstat(bdd_gran,statistic=input$method,aggr=FALSE)$index
        
        output$downloadData <- downloadHandler(
          filename = "output.xlsx",
          content = function(file){
            wb <- createWorkbook()
            sheet  <- createSheet(wb, sheetName="Index")
            addDataFrame(result, sheet)
            saveWorkbook(wb, file)
            #           
          })
        #           
        
      }
      if (input$gran=='Texture') {
        result <-granstat(bdd_gran,statistic=input$method,aggr=FALSE)$sedim
        
        output$downloadData <- downloadHandler(
          filename = "output.xlsx",
          content = function(file){
            wb <- createWorkbook()
            sheet  <- createSheet(wb, sheetName="Texture")
            addDataFrame(result, sheet)
            saveWorkbook(wb, file)
            #           
            #           
          })
      }
      if (input$gran=='All') {
        result <-granstat(bdd_gran,statistic=input$method,aggr=TRUE)
        
        output$downloadData <- downloadHandler(
          filename = "output.xlsx",
          content = function(file){
            wb <- createWorkbook()
            sheet  <- createSheet(wb, sheetName="Statistics")
            addDataFrame(granstat(bdd_gran,statistic=input$method,aggr=FALSE)$stat, sheet)
            sheet  <- createSheet(wb, sheetName="Index")
            addDataFrame(granstat(bdd_gran,statistic=input$method,aggr=FALSE)$index, sheet)
            sheet  <- createSheet(wb, sheetName="Texture")
            addDataFrame(granstat(bdd_gran,statistic=input$method,aggr=FALSE)$sedim, sheet)
            saveWorkbook(wb, file)
            #           
            #           
          })
      }
      
      result 
      
    })
    
    
    output$plot1 <-renderPlot({
      
      par(mfrow=c(length(c(input$from:input$to)),1))  
      for (i in c(input$from:input$to))
        granplot(bdd_gran,xc=i,hist=input$hist,cum=input$cum,main=names(bdd_gran)[i],cexname=1.2)
      
      # output$downloadplot1 <- downloadHandler(
      #   filename = "histo_output.png",
      #   content = function(file) {
      #     png(file)})
      
    },height="auto",width=600) 
    
    output$plot2 <-renderPlot({
      
      grandistrib(bdd_gran,scale=input$distritype,xlab = "")
      
    },height="auto",width=800) 
    
    
    bdd_gran
    
    
    
    #   output$stat <- renderDataTable({
    #     
    #     result_gran <- granstat(input$gran,statistic=method,aggr=aggr)
    #     
    #   })
  })
})


