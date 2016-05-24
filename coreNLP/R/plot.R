#' Print Dependencies and POS of Annotation
#'
#' @importFrom   plotrix boxed.labels draw.arc
#' @method plot  annotation
#' @param x      an annotation object
#' @param y      sentence id
#' @param ...    other arguments passed to plot
#'
#' @export
plot.annotation = function(x,y=1L,...) {
  anno = x
  tkn = getToken(anno)
  dep = getDependency(anno)

  if (is.null(tkn) || is.null(dep))
    stop("Cannot plot without tokens and dependencies!")

  tkn = tkn[tkn$sentence == y,,drop=FALSE]
  dep = dep[dep$sentence == y,,drop=FALSE]

  if (nrow(tkn) < 1 || nrow(dep) < 1)
    stop("Plot needs at least 2 tokens and dependencies.")

  pos = c("DT","CC","CD","EX","FW","POS","PDT","MD","LS","IN","RP","SYM",
        "TO","UH","WDT","WRB",".","JJ","JJR","JJS","NN","NNS","NNP","NNPS",
        "RB","RBR","RBS","VB","VBD","VBG","VBN","VBP","VBZ","WP","WP$",
        "PRP","PRP$")
  posCol = c(rep("#DB9D85",17L),rep("#B1AF64",3L),rep("#6DBC86",4L),
             rep("#39BDBC",3L), rep("#87AEDF",6L),rep("#CD99D8",4L))

  tkn$CharacterOffsetBegin = as.numeric(tkn$CharacterOffsetBegin)
  tkn$CharacterOffsetEnd = as.numeric(tkn$CharacterOffsetEnd)
  tkn$CharacterOffsetBegin = tkn$CharacterOffsetBegin - min(tkn$CharacterOffsetBegin) + 1
  tkn$CharacterOffsetEnd = tkn$CharacterOffsetEnd - min(tkn$CharacterOffsetEnd) + 1
  dep = dep[dep$type != "root",,drop=FALSE]
  dep$governorIdx = as.numeric(as.character(dep$governorIdx))
  dep$dependentIdx = as.numeric(as.character(dep$dependentIdx))
  dep$type = as.character(dep$type)
  mend = max(as.numeric(tkn$CharacterOffsetEnd))

  theseSpace = tkn$CharacterOffsetBegin[-1] - tkn$CharacterOffsetEnd[-nrow(tkn)]
  thisText = tkn$word
  thisText[-length(thisText)][theseSpace > 0] = paste0(thisText[-length(thisText)][theseSpace > 0]," ")
  thisText = paste(thisText,collapse="")

  ASP = 1/2
  LIM = 10
  FACTOR =
  par(mar=c(0,0,0,0))
  plot(0,0,col="white",xlim=c(0,LIM),ylim=c(-2,LIM*0.95),asp=ASP,axes=FALSE,...)
  abline(h=0)

  wd = strwidth(thisText,family="mono")
  CEX = LIM / wd
  ht = strheight(thisText,family="mono",cex=CEX)

  if (3 * ht > 2) {
    CEX = CEX * 1 / ht * 2/3
    ht = strheight(thisText,family="mono",cex=CEX)
  }

  chBegin = tkn$CharacterOffsetBegin / mend * wd * CEX
  chEnd = tkn$CharacterOffsetEnd / mend * wd * CEX
  chMid = (chBegin + chEnd) / 2
  n = length(chBegin)

  index = match(tkn$POS,pos)
  cols = rep("#DB9D85",length(tkn$POS))
  cols[!is.na(index)] = posCol[index[!is.na(index)]]
  text(chBegin,-1.5*ht,tkn$word,pos=4L,offset=0,family="mono",cex=CEX)
  plotrix::boxed.labels(chMid,-3*ht,tkn$POS,family="mono",cex=CEX*0.8,col="black",bg=cols)
  points(chMid,rep(0,length(chBegin)),pch=19,cex=1)

  g = chMid[dep$governorIdx]
  d = chMid[dep$dependentIdx]

  plotrix::draw.arc((g+d)/2,rep(0,length(g)),radius=abs(g-d)/2,angle1=0,angle2=pi)
  text((g+d)/2,abs(g-d)/2 / ASP,dep$type,family="mono",cex=0.7,pos=3L,offset=0.2,col="blue")
  arrow = rep(">", length(g))
  arrow[g > d] = "<"
  text((g+d)/2,abs(g-d)/2 / ASP,arrow,cex=1,col="black")
}