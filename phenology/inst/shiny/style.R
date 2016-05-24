options(rstudio.markdownToHTML =
          function(inputFile, outputFile) {     
            require(markdown)
            markdownToHTML(inputFile, outputFile, stylesheet='Shiny-Phenology.css')  
          }
)