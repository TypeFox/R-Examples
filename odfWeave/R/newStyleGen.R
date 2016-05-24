# This is intended to be the new version of odfStyleGen.
# It is essentially the same, except that instead of creating
# strings containing XML markup, it creates XML node objects.
# This is necessary because of the changes that were made
# to XML processing.  These XML node objects are added to
# the document object, and then the whole thing is written
# out to an XML file.
#
# This will be renamed to odfStyleGen after Max has signed
# off on this whole change.
#

newStyleGen <- function(x, type="styles")
{
   nonempty <- function(x) !is.null(x) && x != ""

   getFontAttrs <- function(styleName, thisStyle)
   {
      c(
         if (nonempty(thisStyle$fontType))
         {
            c(
               if (length(grep("bold", thisStyle$fontType)))
                  c("fo:font-weight"="bold"),
               if (length(grep("italic", thisStyle$fontType)))
                  c("fo:font-style"="italic"),
               if (length(grep("underline", thisStyle$fontType)))
                  c("style:text-underline-style"="solid",
                    "style:text-underline-width"="auto",
                    "style:text-underline-color"="font-color"),
               if (length(grep("shadow", thisStyle$fontType)) == 1)
                  c("fo:text-shadow"="1pt 1pt"),
               if (length(grep("superscript", thisStyle$fontType)) == 1)
                  c("style:text-position"="super 58%"),
               if (length(grep("subscript", thisStyle$fontType)) == 1)
                  c("style:text-position"="sub 58%")
            )
         }
      )
   }

   getStyleAttrs <- function(styleName, thisStyle)
   {
      c(
         "style:name"=styleName,
         "style:family"="paragraph",
         if (nonempty(thisStyle$parentStyleName))
            c("style:parent-style-name"=thisStyle$parentStyleName)
      )
   }

   getTextPropertyAttrs <- function(styleName, thisStyle)
   {
      c(
         if (nonempty(thisStyle$fontColor))
            c("fo:color"=thisStyle$fontColor),
         if (nonempty(thisStyle$fontSize))
            c("fo:font-size"=thisStyle$fontSize),
         getFontAttrs(styleName, thisStyle),
         if (nonempty(thisStyle$fontName))
            c("style:font-name"=thisStyle$fontName)
      )
   }

   getParagraphPropertyAttrs <- function(styleName, thisStyle)
   {
      c(
         if (nonempty(thisStyle$textAlign))
            c("fo:text-align"=thisStyle$textAlign)
      )
   }

   getPageLayoutAttrs <- function(styleName, thisStyle, idx)
   {
      c("style:name"=paste("pm", idx + 100, sep=""))
   }

   getPagePropertiesAttrs <- function(styleName, thisStyle)
   {
      c(
         if (nonempty(thisStyle$marginLeft))
            c("fo:margin-left"=thisStyle$marginLeft),
         if (nonempty(thisStyle$marginRight))
            c("fo:margin-right"=thisStyle$marginRight),
         if (nonempty(thisStyle$marginTop))
            c("fo:margin-top"=thisStyle$marginTop),
         if (nonempty(thisStyle$marginBottom))
            c("fo:margin-bottom"=thisStyle$marginBottom),
         if (nonempty(thisStyle$pageWidth))
            c("fo:page-width"=thisStyle$pageWidth),
         if (nonempty(thisStyle$pageHeight))
            c("fo:page-height"=thisStyle$pageHeight),
         if (nonempty(thisStyle$numFormat))
            c("style:num-format"=thisStyle$numFormat),
         if (nonempty(thisStyle$printOrientation))
            c("style:print-orientation"=thisStyle$printOrientation)
      )
   }

   getHeaderStyleAttrs <- function(styleName, thisStyle)
   {
      character()
   }

   getFooterStyleAttrs <- function(styleName, thisStyle)
   {
      character()
   }

   getMasterPageAttrs <- function(styleName, thisStyle, idx)
   {
      c("style:name"=styleName,
        "style:page-layout-name"=paste("pm", idx + 100, sep=""))
   }

   getTableStyleAttrs <- function(styleName, thisStyle)
   {
      c("style:name"=styleName, "style:family"="table")
   }

   getTablePropertiesAttrs <- function(styleName, thisStyle)
   {
      c(
         if (nonempty(thisStyle$marginLeft))
            c("fo:margin-left"=thisStyle$marginLeft),
         if (nonempty(thisStyle$marginRight))
            c("fo:margin-right"=thisStyle$marginRight),
         if (nonempty(thisStyle$marginTop))
            c("fo:margin-top"=thisStyle$marginTop),
         if (nonempty(thisStyle$marginBottom))
            c("fo:margin-bottom"=thisStyle$marginBottom),
         if (nonempty(thisStyle$align))
            c("table:align"=thisStyle$align)
      )
   }

   getTableCellStyleAttrs <- function(styleName, thisStyle)
   {
      c("style:name"=styleName, "style:family"="table-cell")
   }

   getTableCellPropertiesAttrs <- function(styleName, thisStyle)
   {
      c(
         if (nonempty(thisStyle$verticalAlign))
            c("style:vertical-align"=thisStyle$verticalAlign),
         if (nonempty(thisStyle$leftBorder))
            c("fo:border-left"=thisStyle$leftBorder),
         if (nonempty(thisStyle$rightBorder))
            c("fo:border-right"=thisStyle$rightBorder),
         if (nonempty(thisStyle$topBorder))
            c("fo:border-top"=thisStyle$topBorder),
         if (nonempty(thisStyle$bottomBorder))
            c("fo:border-bottom"=thisStyle$bottomBorder),
         if (nonempty(thisStyle$padding))
            c("fo:padding"=thisStyle$padding),
         if (nonempty(thisStyle$backgroundColor))
            c("fo:background-color"=thisStyle$backgroundColor)
      )
   }

   getFigureFrameStyleAttrs <- function(styleName, thisStyle)
   {
      c("style:name"=styleName,
        "style:family"="graphic",
        "style:parent-style-name"="Frame")
   }

   getFigureFramePropertiesAttrs <- function(styleName, thisStyle)
   {
      c(
         if (nonempty(thisStyle$verticalPosition))
            c("style:vertical-pos"=thisStyle$verticalPosition),
         if (nonempty(thisStyle$verticalRelates))
            c("style:vertical-rel"=thisStyle$verticalRelates),

         if (nonempty(thisStyle$horizontalPosition))
            c("style:horizontal-pos"=thisStyle$horizontalPosition),
         if (nonempty(thisStyle$horizontalRelates))
            c("style:horizontal-rel"=thisStyle$horizontalRelates),

         if (nonempty(thisStyle$backgroundColor))
            c("fo:background-color"=thisStyle$backgroundColor),

         if (nonempty(thisStyle$wrap) && !(thisStyle$frameAnchor %in% c("char", "as-char")))
            c("style:wrap"=thisStyle$wrap),

         if (nonempty(thisStyle$leftBorder))
            c("fo:border-left"=thisStyle$leftBorder),
         if (nonempty(thisStyle$rightBorder))
            c("fo:border-right"=thisStyle$rightBorder),
         if (nonempty(thisStyle$topBorder))
            c("fo:border-top"=thisStyle$topBorder),
         if (nonempty(thisStyle$bottomBorder))
            c("fo:border-bottom"=thisStyle$bottomBorder),
         if (nonempty(thisStyle$padding))
            c("fo:padding"=thisStyle$padding)
      )
   }

   getBulletStyleAttrs <- function(styleName, thisStyle)
   {
      c("style:name"=paste(styleName, "Paragraph", sep=""),  # XXX verify
        "style:family"="paragraph",
        if (nonempty(thisStyle$paraStyle))
        {
           c("style:parent-style-name"=thisStyle$paraStyle)
        } else {
           c("style:parent-style-name"="Standard")
        },
        c("style:list-style-name"=styleName)
      )
   }

   getListStyleAttrs <- function(styleName, thisStyle)
   {
      c("style:name"=styleName)
   }

   getListLevelStyleAttrs <- function(styleName, thisStyle)
   {
      c("text:level"="1",
        "text:style-name"="Bullet_20_Symbols",
        "style:num-suffix"=".",
        if (nonempty(thisStyle$bulletChar))
           c("text:bullet-char"=thisStyle$bulletChar)
      )
   }

   getListLevelPropertiesAttrs <- function(styleName, thisStyle)
   {
      c(
         if (nonempty(thisStyle$spaceBefore))
            c("text:space-before"=thisStyle$spaceBefore),
         if (nonempty(thisStyle$minLabelWidth))
            c("text:min-label-width"=thisStyle$minLabelWidth)
      )
   }

   getTextPropertiesAttrs <- function(styleName, thisStyle)
   {
      c("style:font-name"="StarSymbol")
   }

   styles <- unlist(lapply(x, function(x) x$type))
   styleNames <- names(x)

   if (type == "styles")
   {
      # These styles will go into the 'office:styles' section
      # of styles.xml
      paragraphFun <- function(idx)
      {
         thisStyle <- x[[idx]]
         styleName <- styleNames[idx]

         pattrs <- getParagraphPropertyAttrs(styleName, thisStyle)
         children <- c(
            if (!is.null(pattrs))
               list(xmlNode('style:paragraph-properties',
                            attrs=pattrs)),
            list(xmlNode('style:text-properties',
                         attrs=getTextPropertyAttrs(styleName, thisStyle)))
         )
         xmlNode('style:style',
                 attrs=getStyleAttrs(styleName, thisStyle),
                 .children=children)
      }
      idx <- which(styles == 'Paragraph')
      lapply(idx, paragraphFun)
   } else if (type == "page") {
      # These styles will go into the 'office:automatic-styles' section
      # of styles.xml
      pageFun <- function(idx)
      {
         thisStyle <- x[[idx]]
         styleName <- styleNames[idx]
         children <- c(
            list(xmlNode('style:page-layout-properties',
                         attrs=getPagePropertiesAttrs(styleName, thisStyle))),
            list(xmlNode('style:header-style',
                         attrs=getHeaderStyleAttrs(styleName, thisStyle))),
            list(xmlNode('style:footer-style',
                         attrs=getFooterStyleAttrs(styleName, thisStyle)))
         )
         xmlNode('style:page-layout',
                 attrs=getPageLayoutAttrs(styleName, thisStyle, idx),
                 .children=children)
      }
      idx <- which(styles == 'Page')
      lapply(idx, pageFun)
   } else if (type == "master") {
      # These styles will go into the 'office:master-styles' section
      # of styles.xml
      masterFun <- function(idx)
      {
         thisStyle <- x[[idx]]
         styleName <- styleNames[idx]
         xmlNode('style:master-page',
                 attrs=getMasterPageAttrs(styleName, thisStyle, idx))
      }
      idx <- which(styles == 'Page')
      lapply(idx, masterFun)
   } else if (type == "content") {
      # These styles will go into the 'office:automatic-styles' section
      # of content.xml

      # Get the list of table styles
      tableFun <- function(idx)
      {
         thisStyle <- x[[idx]]
         styleName <- styleNames[idx]
         children <- list(xmlNode('style:table-properties',
                                  attrs=getTablePropertiesAttrs(styleName, thisStyle)))
         xmlNode('style:style',
                 attrs=getTableStyleAttrs(styleName, thisStyle),
                 .children=children)
      }
      idx <- which(styles == 'Table')
      tableStyleNodes <- lapply(idx, tableFun)

      # Get the list of table cell styles
      tableCellFun <- function(idx)
      {
         thisStyle <- x[[idx]]
         styleName <- styleNames[idx]
         children <- list(xmlNode('style:table-cell-properties',
                                  attrs=getTableCellPropertiesAttrs(styleName, thisStyle)))
         xmlNode('style:style',
                 attrs=getTableCellStyleAttrs(styleName, thisStyle),
                 .children=children)
      }
      idx <- which(styles == 'Table Cell')
      tableCellStyleNodes <- lapply(idx, tableCellFun)

      # Get the list of figure frame/graphic styles
      figureFrameFun <- function(idx)
      {
         thisStyle <- x[[idx]]
         styleName <- styleNames[idx]
         children <- list(xmlNode('style:graphic-properties',
                                  attrs=getFigureFramePropertiesAttrs(styleName, thisStyle)))
         xmlNode('style:style',
                 attrs=getFigureFrameStyleAttrs(styleName, thisStyle),
                 .children=children)
      }
      idx <- which(styles == 'Figure Frame')
      figureFrameStyleNodes <- lapply(idx, figureFrameFun)

      # Get the list of bullet styles used for lists
      bulletFun <- function(idx)
      {
         thisStyle <- x[[idx]]
         styleName <- styleNames[idx]
         xmlNode('style:style',
                 attrs=getBulletStyleAttrs(styleName, thisStyle))
      }
      idx <- which(styles == 'Bullet List')
      bulletStyleNodes <- lapply(idx, bulletFun)

      # Get the list of list styles
      listFun <- function(idx)
      {
         thisStyle <- x[[idx]]
         styleName <- styleNames[idx]
         children <-
            list(xmlNode('text:list-level-style-bullet',
                         attrs=getListLevelStyleAttrs(styleName, thisStyle),
                         .children=list(xmlNode('style:list-level-properties',
                                               attrs=getListLevelPropertiesAttrs(styleName, thisStyle)))))
         xmlNode('text:list-style',
                 attrs=getListStyleAttrs(styleName, thisStyle),
                 .children=children)
      }
      idx <- which(styles == 'Bullet List')
      listStyleNodes <- lapply(idx, listFun)

      # Make a list with the one style used to support the odfPageBreak function
      breakchildren <- xmlNode('style:paragraph-properties',
                               attrs=c('fo:break-before'='page'))
      breaknode <- xmlNode('style:style',
                           attrs=c('style:parent-style-name'='Standard',
                                   'style:name'='OdfPageBreak',
                                   'style:family'='paragraph'),
                           .children=list(breakchildren))
      odfPageBreak <- list(odfPageBreak=breaknode)

      # Concatenate them all together into one big list to return
      c(tableStyleNodes, tableCellStyleNodes, figureFrameStyleNodes,
        bulletStyleNodes, listStyleNodes, odfPageBreak)
   } else {
     stop("illegal type specified: ", type)
   }
}
