#################################################################################
#
# Note:
#
#   This file is currently unused.  The functions in newStyleGen.R are
#   being used instead.  This will be replaced by newStyleGen.R eventually.
#
#################################################################################

#################################################################################
# This file contains functions to write out style infomration into either
# content.xml or styles.xml. At the time of creation of these functions, there
# was no good resource to tell which styles wen into which file (besides trial 
# and error). 


# tagattr is a function that will write out style parameters, such as
#
#   cat(tagattr("style:vertical-align", "top"))
#
# creates
#
#   style:vertical-align="top"

"tagattr" <- function(name, val) paste(c(name, '="', val, '"'), collapse="")



"element" <- function(tag, atts, content="") {
   #tag:  character string, name the element
   #atts: character, attributes of a tag in the form name="value"
	#   empty items are discarded 
   #content: character string, text of an element
   #value: character string, and properly-formed element
   atts <- paste(atts[atts!=""], collapse=" ")
   if (all(content == "")) {
      close <- "/>"
      endtag <- ""
   } else {
      close <- ">"
      endtag <- paste(c('</', tag, '>'), collapse="")
   }
   paste(c('<', tag, ' ', atts, close, content, endtag), collapse="" )
}


# odfStyleGen is what is called to create the styles. x is a list of
# styles via getStyleDefs and type is the destination file.

"odfStyleGen" <- function(x, type = "styles") {
   out <- ""
   if(is.null(x)) return(out)
   
   "has" <- function(x) !is.null(x) && x != ""
   
   styles <-  unlist(lapply(x, function(x) x$type))
   
   # for styles that go into styles.xml...
   if(type == "styles") {
   
      # get any paragraph styles
      
      paraStyles <- styles[styles == "Paragraph"]
      
      for(i in seq(along = paraStyles)) {
         thisStyle <- x[[names(paraStyles)[i]]]

         # for each paragraph style, write out appropriate xml
         if(has(thisStyle$fontType)) {
            fontText <- ""
            if(length(grep("bold", thisStyle$fontType)))
               fontText <- c(fontText, tagattr("fo:font-weight", "bold"))
            if(length(grep("italic", thisStyle$fontType)))
               fontText <- c(fontText, tagattr("fo:font-style", "italic"))
            if(length(grep("underline", thisStyle$fontType)))
               fontText <- c(
                  fontText,
                  tagattr("style:text-underline-style", "solid"),
                  tagattr("style:text-underline-width", "auto"),
                  tagattr("style:text-underline-color", "font-color"))
            if(length(grep("shadow", thisStyle$fontType)) == 1)
               fontText <- c(fontText, tagattr("fo:text-shadow", "1pt 1pt"))
            if(length(grep("superscript", thisStyle$fontType)) == 1)
               fontText <- c(
                  fontText, tagattr("style:text-position", "super 58%"))
            if(length(grep("subscript", thisStyle$fontType)) == 1)
               fontText <- c(
                  fontText, tagattr("style:text-position", "sub 58%"))
         } else fontText <- ""

         style_style <- 'style:style'
         style_style_attr <- c(
            tagattr('style:name', paste(names(paraStyles)[i], collapse=" ")),
            tagattr('style:family', 'paragraph'),
            if(has(thisStyle$parentStyleName))
               tagattr('style:parent-style-name', thisStyle$parentStyleName)
         )
         style_text <- 'style:text-properties'
         style_text_attr <- c(
            if(has(thisStyle$fontColor))
               tagattr('fo:color', thisStyle$fontColor),
            if(has(thisStyle$fontSize))
               tagattr('fo:font-size', thisStyle$fontSize),
            fontText,
            if(has(thisStyle$fontName))
               tagattr('style:font-name', thisStyle$fontName))
         style_paragraph <- 'style:paragraph-properties'
         style_paragraph_attr <- c(
            tagattr('fo:text-align', thisStyle$textAlign))
         out <- paste(
            c(
            out,
            element(style_style, style_style_attr,
               c(
                  element(
                     style_text,
                     style_text_attr
                  ),
                  if(has(thisStyle$textAlign))
                     element(style_paragraph, style_paragraph_attr)
                  )
               )
            ),
            collapse="\n"
         )
      }

   } else if (type == "page") {
      pageStyles <- styles[styles == "Page"]
      for (i in seq(along = pageStyles))
      {
         thisStyle <- x[[names(pageStyles)[i]]]
         style_style <- 'style:page-layout'
         style_style_attr <- c(
            tagattr("style:name", paste("pm", i + 100, sep="")))
         page_properties <- 'style:page-layout-properties'
         page_properties_attr <- c(
            if(has(thisStyle$marginLeft))
               tagattr("fo:margin-left", thisStyle$marginLeft),
            if(has(thisStyle$marginRight))
               tagattr("fo:margin-right", thisStyle$marginRight),
            if(has(thisStyle$marginTop))
               tagattr("fo:margin-top", thisStyle$marginTop),
            if(has(thisStyle$marginBottom))
               tagattr("fo:margin-bottom", thisStyle$marginBottom),
            if(has(thisStyle$pageWidth))
               tagattr("fo:page-width", thisStyle$pageWidth),
            if(has(thisStyle$pageHeight))
               tagattr("fo:page-height", thisStyle$pageHeight),
            if(has(thisStyle$numFormat))
               tagattr("style:num-format", thisStyle$numFormat),
            if(has(thisStyle$printOrientation))
               tagattr("style:print-orientation", thisStyle$printOrientation)
         )
         style_header_style <- 'style:header-style'
         style_header_style_attr <- c()
         style_footer_style <- 'style:footer-style'
         style_footer_style_attr <- c()
         out <- paste(
            c(
               out,
               element(style_style, style_style_attr,
                  c(
                     element(page_properties, page_properties_attr),
                     element(style_header_style, style_header_style_attr),
                     element(style_footer_style, style_footer_style_attr)))
            ),
            collapse="\n")
      }
   } else if (type == "master") {
      pageStyles <- styles[styles == "Page"]
      for (i in seq(along = pageStyles))
      {
         thisStyle <- x[[names(pageStyles)[i]]]
         style_style <- 'style:master-page'
         style_style_attr <- c(
            tagattr("style:name", paste(names(pageStyles)[i], collapse=" ")),
            tagattr("style:page-layout-name", paste("pm", i + 100, sep="")))
         out <- paste(
            c(
               out,
               element(style_style, style_style_attr)
            ),
            collapse="\n")
      }
   } else if (type == "content") {
      tableStyles <- styles[styles == "Table"]
      for(i in seq(along = tableStyles)) {
         thisStyle <- x[[names(tableStyles)[i]]]
         style_style <- 'style:style'
         style_style_attr <- c(
            tagattr("style:name", paste(names(tableStyles)[i], collapse=" ")),
            tagattr("style:family", "table"))
         table_properties <- 'style:table-properties'
         table_properties_attr <- c(
            if(has(thisStyle$marginLeft))
               tagattr("fo:margin-left", thisStyle$marginLeft),
            if(has(thisStyle$marginRight))
               tagattr("fo:margin-right", thisStyle$marginRight),
            if(has(thisStyle$marginTop))
               tagattr("fo:margin-top", thisStyle$marginTop),
            if(has(thisStyle$marginBottom))
               tagattr("fo:margin-bottom", thisStyle$marginBottom),
            if(has(thisStyle$align))
               tagattr("table:align", thisStyle$align))
         out <- paste(
            c(
               out,
               element(style_style, style_style_attr,
                  element(table_properties, table_properties_attr))
            ),
            collapse="\n")
      }

      cellStyles <- styles[styles == "Table Cell"]
      for(i in seq(along = cellStyles))
      {
         thisStyle <- x[[names((cellStyles)[i])]]
         style_style <- "style:style"
         style_style_attr <- c(
            tagattr("style:name", paste(names(cellStyles)[i], collapse=" ")),
            tagattr("style:family", "table-cell"))
         table_cell <- "style:table-cell-properties"
         table_cell_attr <- c(
            if(has(thisStyle$verticalAlign))
               tagattr("style:vertical-align", thisStyle$verticalAlign),
            if(has(thisStyle$leftBorder))
               tagattr("fo:border-left", thisStyle$leftBorder),
            if(has(thisStyle$rightBorder))
               tagattr("fo:border-right", thisStyle$rightBorder),
            if(has(thisStyle$topBorder))
               tagattr("fo:border-top", thisStyle$topBorder),
            if(has(thisStyle$bottomBorder))
               tagattr("fo:border-bottom", thisStyle$bottomBorder),
            if(has(thisStyle$padding))
               tagattr("fo:padding", thisStyle$padding),
            if(has(thisStyle$backgroundColor))
               tagattr("fo:background-color", thisStyle$backgroundColor)
         )
         out <- paste(
            c(
               out,
               element(
                  style_style,
                  style_style_attr,
                  element(
                     table_cell,
                     table_cell_attr))),
            collapse="\n")
      }
      
      
      figFrameStyles <- styles[styles == "Figure Frame"]
      for(i in seq(along = figFrameStyles)) {
         thisStyle <- x[[names(figFrameStyles)[i]]]
         
         style_style <- 'style:style'
         style_style_attr <- c(
            tagattr("style:name", paste(names(figFrameStyles)[i], collapse=" ")),
            tagattr("style:family", "graphic"),
            tagattr("style:parent-style-name", "Frame"))
            
         figFrame_properties <- 'style:graphic-properties'
         
         figFrame_properties_attr <- c(
         
            if(has(thisStyle$verticalPosition))
               tagattr("style:vertical-pos", thisStyle$verticalPosition),
            if(has(thisStyle$verticalRelates))
               tagattr("style:vertical-rel", thisStyle$verticalRelates),               
         
            if(has(thisStyle$horizontalPosition))
               tagattr("style:horizontal-pos", thisStyle$horizontalPosition),
            if(has(thisStyle$horizontalRelates))
               tagattr("style:horizontal-rel", thisStyle$horizontalRelates),                   
               
            if(has(thisStyle$backgroundColor))
               tagattr("fo:background-color", thisStyle$backgroundColor),

            if(has(thisStyle$wrap) & !(thisStyle$frameAnchor %in% c("char", "as-char")))
               tagattr("style:wrap", thisStyle$wrap),
                
            if(has(thisStyle$leftBorder))
               tagattr("fo:border-left", thisStyle$leftBorder),
            if(has(thisStyle$rightBorder))
               tagattr("fo:border-right", thisStyle$rightBorder),
            if(has(thisStyle$topBorder))
               tagattr("fo:border-top", thisStyle$topBorder),
            if(has(thisStyle$bottomBorder))
               tagattr("fo:border-bottom", thisStyle$bottomBorder),
            if(has(thisStyle$padding))
               tagattr("fo:padding", thisStyle$padding)
            )
         out <- paste(
            c(
               out,
               element(style_style, style_style_attr,
                  element(figFrame_properties, figFrame_properties_attr))
            ),
            collapse="\n")
      }      
      
      bulletStyles <- styles[styles == "Bullet List"]
      for(i in seq(along = bulletStyles))
      {
         thisStyle <- x[[names(bulletStyles)[i]]]      
         startParaTag <- "  <style:style "
         listParaName <- tagattr("style:name", paste(names(bulletStyles)[i], "Paragraph", sep=""))
         listParaStyleFamily <- 'style:family="paragraph"'
         listParaParentStyle <- if(has(thisStyle$paraStyle)) tagattr("style:parent-style-name", thisStyle$paraStyle) else tagattr("style:parent-style-name", "Standard")
         listParaStyleName <- tagattr("style:list-style-name", names(bulletStyles)[i])
         endParaListTag <- "/>\n"
         listPart1 <- c(
            '  <text:list-style',
            tagattr("style:name", names(bulletStyles)[i]),
            ">\n") 
         listPart2 <- c(
            '   <text:list-level-style-bullet text:level=\"1\" text:style-name=\"Bullet_20_Symbols\" style:num-suffix=\".\"',
            if(has(thisStyle$bulletChar)) tagattr("text:bullet-char", thisStyle$bulletChar),
            ">\n")
         listPart3 <- c(
            '    <style:list-level-properties ',
            if(has(thisStyle$spaceBefore)) tagattr("text:space-before", thisStyle$spaceBefore),
            if(has(thisStyle$minLabelWidth)) tagattr("text:min-label-width", thisStyle$minLabelWidth),
            "/>\n")
         listPart4 <- c(
            '    <style:text-properties style:font-name=\"StarSymbol\"/>',
            '   </text:list-level-style-bullet>',
            '  </text:list-style>\n')
         out <- c(out,
            paste(
               paste(startParaTag, listParaName, listParaStyleFamily, listParaParentStyle, listParaStyleName, endParaListTag),
               paste(listPart1, collapse = " "),
               paste(listPart2, collapse = " "),
               paste(listPart3, collapse = " "),
               paste(listPart4, collapse = "\n"),
               collapse = "\n"))
               
      }      
   } else {
      stop('illegal type specified')
   }
   gsub(">", ">\n", out)
}
