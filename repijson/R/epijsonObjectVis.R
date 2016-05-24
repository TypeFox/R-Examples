if(getRversion() >= "2.15.1")  utils::globalVariables(c("x","y","record"))

#'  to view the structure of epijson objects and/or schema

#' allowing multiple attribute boxes with different labels

#' @param attribMeta label(s) for Metadata attributes
#' @param attribRecord label(s) for Record attributes
#' @param attribEvent label(s) for Event attributes
#' @param labelObject label for the overall object
#' @param labelMeta label for Metadata attributes
#' @param labelRecord label for Record attributes
#' @param labelEvent label for Event attributes
#' @param colAll optional single box colour to overide all other col args, e.g. 'grey'
#' @param colObject object box colour
#' @param colMeta metadata box colour
#' @param colRecord record box colour
#' @param colEvent event box colour
#' @param colAttrib attribute boxes colour
#' @param textSize size of labels default=4
#'
#' @return a ggplot object
#' @examples
#' #this gives the base schema
#' epijsonObjectVis()
#' #settin single box colour and increasing text size
#' epijsonObjectVis(colAll ='grey', textSize=7)
#' #this gives a diagram with named attributes
#' epijsonObjectVis( attribMeta = c("attribute: disease","attribute: data provider"),
#'                    attribRecord = c("attribute: gender","attribute: date of birth"),
#'                    attribEvent = c("attribute: recorder","attribute: test used") )
#' epijsonObjectVis( attribMeta = c("a"),
#'                    attribRecord = c("b","c"),
#'                    attribEvent = c("d","e","f") )
#' #repijson objects                    
#' epijsonObjectVis( attribMeta = 'ejAttribute',
#'                   attribRecord = 'ejAttribute',
#'                   attribEvent = 'ejAttribute',
#'                   labelObject = 'ejObject',
#'                   labelMeta = 'ejMetadata',
#'                   labelRecord = 'ejRecord',
#'                   labelEvent = 'ejEvent')                    
#' #repijson objects and constructors
#' epijsonObjectVis( attribMeta = 'ejAttribute create_ejAttribute()',
#'                    attribRecord = 'ejAttribute create_ejAttribute()',
#'                    attribEvent = 'ejAttribute create_ejAttribute()',
#'                    labelObject = 'repijson R objects and constructors : ejObject create_ejObject()',
#'                    labelMeta = 'ejMetadata create_ejMetadata()',
#'                    labelRecord = 'ejRecord create_ejRecord()',
#'                    labelEvent = 'ejEvent create_ejEvent()')
#' @export
epijsonObjectVis <- function( attribMeta = 'attributes [name, type, value, units]',
                               attribRecord = 'attributes [name, type, value, units]',
                               attribEvent = 'attributes [name, type, value, units]',
                               labelObject = 'Diagram of EpiJSON file structure',
                               labelMeta = 'metadata',
                               labelRecord = 'records [id]',
                               labelEvent = 'events [id, name, date, location]',
                               colAll = NULL,
                               colObject = "gray60",
                               colMeta = "gray60",
                               colRecord = "blue",
                               colEvent = "red",
                               colAttrib = "purple",
                               textSize = 4
                               )

{


  #count number of attributes
  nattM <- length(attribMeta)
  nattR <- length(attribRecord)
  nattE <- length(attribEvent)
  #todo sort this
  #assume that if num atts is 1, show multiple sheets, otherwise just show indiv sheets
  #(this gives misleading pic if there is indeed just one attribute)
  sheetsAttM <- sheetsAttR <- sheetsAttE <- 2
  if (nattM>1) sheetsAttM <- 0
  if (nattR>1) sheetsAttR <- 0
  if (nattE>1) sheetsAttE <- 0

  #set all box colours to a single one if colAll is set
  if (!is.null(colAll)) colObject <- colMeta <- colRecord <- colEvent <- colAttrib <- colAll
  

  xatt <- 0.65 #length att box
  yatt <- 0.06 #height att box (excl 'sheets')
  ysheets <- xsheets <- 0.02
  #could us these spacings below too
  xgap <- 0.04
  ygap <- 0.02
  #todo below it's not quite consistent whether ylab includes gap or not
  ylab <- 0.04  #spacing for label text, includes gap could base on font size

  #spacing for containers
  yM <- ylab + nattM*(yatt+ygap) + sheetsAttM*ysheets
  yE <- ylab + nattE*(yatt+ygap) + sheetsAttE*ysheets #+ ygap
#   #yR for Record (not records)
#   yR <- ylab + nattR*(yatt+ygap) + sheetsAttR*ysheets + yE + ysheets + ygap + ygap #todo this last ygap is a fudge because addSheets is 2 when it should be 1 for sizing
#   yRs <- ylab+yR+ygap+ysheets+ygap
  #remove outer records box
  yRs <- ylab + nattR*(yatt+ygap)+ sheetsAttR*ysheets + yE + ysheets + ygap + ygap #todo this last ygap is a fudge because addSheets is 2 when it should be 1 for sizing
  

  yAll <- ylab + ygap + yM + ygap + yRs + 2*ysheets + ygap

  #cat("yatt:",yatt," yM:",yM," yE:",yE,"yR:",yR,"\n")

  #ymax set to total height of plot
  xmin=0.01; xmax=0.99; ymin=0.01; ymax=ymin+yAll

  #todo : gg use could be refactored with pipes %>%

  gg <- box2(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, label=labelObject, colour=colObject, size=textSize)

  ## Metadata
  #work down, set ymax first, then ymin from ymax
  xmin=xmin+xgap; xmax=xmax-xgap; ymax=ymax-(ylab+ygap); ymin=ymax-yM;
  gg <- box2(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, label=labelMeta, gg=gg, colour=colMeta, size=textSize)
  #saving metadata dims to use later
  mxmin=xmin; mxmax=xmax; mymin=ymin; mymax=ymax

#   ##Records
#   #set ymax from ymin of Metadata
#   xmin=xmin; xmax=xmax; ymax=ymin-ygap; ymin=ymax-yRs
#   gg <- box2(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, label='Records', gg=gg, colour=colRecords, size=textSize)
# 
#   ##Record
#   #record is simpler because it just fits in records
#   xmin=xmin+xgap; xmax=xmax-xgap; ymax=ymax-(ylab+ygap); ymin=ymin+(ygap+ysheets);
#   gg <- box2(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, label=labelRecord, gg=gg, addSheets=2, colour=colRecord, size=textSize)
  
  #replacing record and records with just records
  ##Records
  #set ymax from ymin of Metadata
  xmin=xmin; xmax=xmax; ymax=ymin-ygap; ymin=ymax-yRs
  gg <- box2(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, label=labelRecord, gg=gg, addSheets=2, colour=colRecord, size=textSize)
  
  
  
  #saving record dims to use later
  rxmin=xmin; rxmax=xmax; rymin=ymin; rymax=ymax

  #Metadata attributes
  for(i in 1:nattM)
  {
    xmin=mxmin+xgap; xmax=xmin+xatt; ymax=mymax-i*(yatt); ymin=ymax-yatt;
    gg <- box2(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, label=attribMeta[i], gg=gg, addSheets=sheetsAttM, colour=colAttrib, size=textSize)
  }
  #Record attributes
  for(i in 1:nattR)
  {
    xmin=rxmin+xgap; xmax=xmin+xatt; ymax=rymax-i*(yatt); ymin=ymax-yatt;
    gg <- box2(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, label=attribRecord[i], gg=gg, addSheets=sheetsAttR, colour=colAttrib, size=textSize)
  }

  #Event needs to go below the last attribute of Record
  xmin=rxmin+xgap; xmax=rxmax-(xgap+xsheets); ymax=ymin-ygap-(sheetsAttR*ysheets); ymin=ymax-yE
  gg <- box2(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, label=labelEvent, gg=gg, addSheets=2, colour=colEvent, size=textSize)
  #saving event dims to use later
  exmin=xmin; exmax=xmax; eymin=ymin; eymax=ymax

  #Event attributes
  for(i in 1:nattE)
  {
    xmin=exmin+xgap; xmax=xmin+xatt; ymax=eymax-i*(yatt); ymin=ymax-yatt;
    gg <- box2(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, label=attribEvent[i], gg=gg, addSheets=sheetsAttE, colour=colAttrib, size=textSize)
  }

  ## remove axes
  gg <- gg + scale_x_continuous("", breaks=NULL)  + scale_y_continuous("", breaks=NULL)


  print(gg)

  invisible(gg)
}


#' generic labelled box function

#' @param xmin 0 to 1
#' @param xmax 0 to 1
#' @param ymin 0 to 1
#' @param ymax 0 to 1
#' @param label to print at top of box
#' @param colour border of box
#' @param fill fill colour for box NA=none
#' @param gg ggplot object, if passed box is added if NULL new ggplot object created
#' @param addSheets how many sheets to add behind box to indicate array
#' @param size font size for label
#' @param print TRUE/FALSE whether to print the ggplot object
#'
#' @return a ggplot object
#'
box2 <- function(xmin=0.01,
                 xmax=0.99,
                 ymin=0.01,
                 ymax=0.99,
                 label='box2',
                 colour="gray60",
                 fill="white",
                 gg=NULL,
                 addSheets=0,
                 size=4,
                 print=FALSE){

  #require(ggplot2)
  if(is.null(gg)) gg <- ggplot()
  #or for easier blank background plotting
  #require(cowplot) #cowplot on github
  #if(is.null(gg)) gg <- ggdraw()

  boxes <- data.frame(
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    label = label
  )

  #add label positions (todo allow args)
  boxes$x <- boxes$xmin + 0.01
  boxes$y <- boxes$ymax - 0.04

  #adding sheets (blank boxes behind)
  if (addSheets>0)
  {
    #copy coords
    sheets <- boxes
    for(sheetNum in addSheets:1)
    {
      #cat(sheet)
      sheets$xmin<-boxes$xmin+sheetNum*0.01
      sheets$xmax<-boxes$xmax+sheetNum*0.01
      sheets$ymin<-boxes$ymin-sheetNum*0.01
      sheets$ymax<-boxes$ymax-sheetNum*0.01
      gg <- gg + geom_rect(data=sheets, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),colour=colour, fill=fill)
    }
  }

  gg <- gg + geom_rect(data=boxes, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), colour=colour, fill=fill)

  gg <- gg + geom_text(data=boxes, aes(x=x, y=y, label=label, hjust=0, vjust=0),size=size)

  if(print) print(gg)

  invisible(gg)

}




#ggsave("filename.pdf")





