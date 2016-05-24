## creates a notebook interface tohandle plots
setClass("gGraphicsNotebooktcltk",
         representation=representation(
           width="numeric",height="numeric"
           ),
         contains="gNotebooktcltk",
         prototype=prototype(new("gNotebooktcltk"))
         )
setMethod(".ggraphicsnotebook",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   width=dpi*6, height=dpi*6,dpi=75,
                   container = NULL,
                   ...) {
            ## ... passed onto gnotebook


            
            force(toolkit)
            
            return(glabel("No ggraphics available in gWidgetstcltk", container=container)@widget)
            
          })
          
