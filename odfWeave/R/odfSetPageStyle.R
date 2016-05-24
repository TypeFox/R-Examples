odfSetPageStyle <- function(style="Standard")
{
   # Allocate a unique name for a new paragraph style
   styleNameEnv <- get('styleNameEnv', pos=.odfEnv, inherits=FALSE)
   paragraphStyleName <- uniqueName('P', styleNameEnv)

   # Create the new paragraph style that we need, and stash it
   # in the "New Style Environment" so it can be put into the
   # document during post processing.
   newStyle <- makeSetPageStyle(paragraphStyleName, style,
                                family='paragraph', type='common',
                                prevstyle='Standard')

   newStyleEnv <- get('newStyleEnv', pos=.odfEnv, inherits=FALSE)
   assign(paragraphStyleName, newStyle, pos=newStyleEnv)

   # Generate the XML for the new paragraph that uses the newly
   # generated paragraph style that will use the specified page style.
   x <- list(text=sprintf('<text:p text:style-name="%s"/>', paragraphStyleName))
   return(structure(x, class='odfPageBreak'))
}
