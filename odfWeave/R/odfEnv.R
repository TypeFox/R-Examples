## the idea for this strucure came form the lattice package

.odfEnv <- new.env()

## Bullet character:
bullet <- "\342\227\217"
Encoding(bullet) <- "UTF-8"

## todo: style values should be validated

assign(
       "picPath",
       NA,
       pos = .odfEnv)

assign(
       "imageDefs",
       list(

            type = "png",
            device = if(capabilities("png")) "png" else "bitmap",
            plotHeight = if(capabilities("png")) 480 else 480/72,
            plotWidth = if(capabilities("png")) 480 else 480/72,
            dispHeight = 5,
            dispWidth = 5,
            args = list()
            ),
       pos = .odfEnv)


assign(
       "odfStyles",
       list(
            paragraph = "ArialNormal",
            input = "ttRed",
            output = "ttBlue",
            table = "RTable1",
            cell = "noBorder",
            header = "lowerBorder",
            cellText = "ArialCentered",
            headerText = "ArialCenteredBold",
            bullet = "Rbullet",
            figureFrame = "basicFigFrame",
            page = "RlandscapePage"
            ),
       pos = .odfEnv)

assign(
       "styleDefs",
       list(

            ## There are distinct text and paragraph properties that can be set (with some
            ## overlapping elements, such as background color). To specifiy text styles,
            ## these two sets of properties are specified under the "Paragraph" style here.
            ArialCenteredBold = list(
              type = "Paragraph",
              parentStyleName = "",
              textAlign = "center",
              fontName = "Arial",
              fontSize = "12pt",
              fontType = "bold",
              fontColor = "#000000"),
            ArialNormal = list(
              type = "Paragraph",
              parentStyleName = "",
              textAlign = "left",
              fontName = "Arial",
              fontSize = "12pt",
              fontType = "normal",
              fontColor = "#000000"),
            ArialCentered = list(
              type = "Paragraph",
              parentStyleName = "",
              textAlign = "center",
              fontName = "Arial",
              fontSize = "12pt",
              fontType = "normal",
              fontColor = "#000000"),
            ArialHighlight = list(
              type = "Paragraph",
              parentStyleName = "",
              textAlign = "center",
              fontName = "Arial",
              fontSize = "12pt",
              fontType = "bold",
              fontColor = "#ff0000"),
            ttBlue = list(
              type = "Paragraph",
              parentStyleName = "",
              textAlign = "left",
              fontName = "Courier New",
              fontSize = "10pt",
              fontType = "normal",
              fontColor = "#000080"),
            ttRed = list(
              type = "Paragraph",
              parentStyleName = "",
              textAlign = "left",
              fontName = "Courier New",
              fontSize = "10pt",
              fontType = "normal",
              fontColor = "#800000"),

            ## Cell specifications are also allowed to include text and paragraph
            ## properties. The "Table Cell" style will not include these properties

            noBorder = list(
              type = "Table Cell",
              backgroundColor="transparent",
              padding = "0.0382in",
              verticalAlign = "automatic",
              padding = "0.0382in",
              leftBorder = "none",
              rightBorder = "none",
              topBorder = "none",
              bottomBorder = "none"),

            lowerBorder = list(
              type = "Table Cell",
              backgroundColor="#FFFFFF",
              padding = "0.0382in",
              verticalAlign = "automatic",
              leftBorder = "none",
              rightBorder = "none",
              topBorder = "none",
              bottomBorder = "0.0007in solid #000000"),



            RTable1 = list(
              type = "Table",
              ##background style
              marginLeft = "0.05in",
              marginRight = "0.05in",
              marginTop = "0.05in",
              marginBottom = "0.05in",
              align = "margins"),

            Rbullet = list(
              type = "Bullet List",
              paraStyle = "ArialNormal",
              bulletChar= bullet,
              spaceBefore="0.25in",
              minLabelWidth="0.25in"),

            basicFigFrame = list(
              type = "Figure Frame",
              verticalPosition = "from-top",
              verticalRelates = "paragraph",
              horizontalPosition = "center",
              horizontalRelates = "paragraph",
              frameAnchor = "paragraph",
              imageAnchor = "paragraph",         
              wrap = "none",
              backgroundColor="transparent",
              padding = "0.02in",
              leftBorder = "0.0008in solid #ffffff",
              rightBorder = "0.0008in solid #ffffff",
              topBorder = "0.0008in solid #ffffff",
              bottomBorder = "0.0008in solid #ffffff"),

            RlandscapePage = list(
              type = "Page",
              printOrientation = "landscape",
              numFormat = "1",
              pageWidth = "11in",
              pageHeight = "8.5in",
              marginLeft = "1.25in",
              marginRight = "1.25in",
              marginTop = "1in",
              marginBottom = "1in")),
       
       pos = .odfEnv)
