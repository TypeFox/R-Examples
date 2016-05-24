#' Box-arrow diagrams for multi-state models.
#' 
#' Automated plotting of the states and transitions that characterize a multi
#' states model.
#' 
#' 
#' @param x An object of class \code{Hist}.
#' @param nrow the number of graphic rows
#' @param ncol the number of graphic columns
#' @param stateLabels Vector of names to appear in the boxes (states).
#' Defaults to attr(x,"state.names").  The boxes can also be individually
#' labeled by smart arguments of the form \code{box3.label="diseased"}, see
#' examples.
#' @param arrowLabels Vector of labels to appear in the boxes (states). One for
#' each arrow.  The arrows can also be individually labeled by smart arguments
#' of the form \code{arrow1.label=paste(expression(eta(s,u)))}, see examples.
#' @param arrowLabelStyle Either "symbolic" for automated symbolic arrow
#' labels, or "count" for arrow labels that reflect the number of transitions
#' in the data.
#' @param arrowLabelSymbol Symbol for automated symbolic arrow labels. Defaults
#' to "lambda".
#' @param changeArrowLabelSide A vector of mode logical (TRUE,FALSE) one for
#' each arrow to change the side of the arrow on which the label is placed.
#' @param tagBoxes Logical. If TRUE the boxes are numbered in the upper left
#' corner. The size can be controlled with smart argument boxtags.cex. The
#' default is boxtags.cex=1.28.
#' @param startCountZero Control states numbers for symbolic arrow labels and
#' box tags.
#' @param oneFitsAll If \code{FALSE} then boxes have individual size, depending
#' on the size of the label, otherwise all boxes have the same size dependent
#' on the largest label.
#' @param margin Set the figure margin via \code{par(mar=margin)}. Less than 4
#' values are repeated.
#' @param cex Initial cex value for the state and the arrow \code{labels}.
#' @param verbose If TRUE echo various things.
#' @param \dots Smart control of arguments for the subroutines text (box
#' label), rect (box), arrows, text (arrow label). Thus the three dots can be
#' used to draw individual boxes with individual labels, arrows and arrow
#' labels. E.g. arrow2.label="any label" changes the label of the second arrow.
#' See examples.
#' @note Use the functionality of the unix program `dot'
#' http://www.graphviz.org/About.php via R package Rgraphviz to obtain more
#' complex graphs.
#' @author Thomas A Gerds \email{tag@@biostat.ku.dk}
#' @seealso \code{\link{Hist}}\code{\link{SmartControl}}
#' @keywords survival
##' @examples
##' 
##' 
##' ## A simple survival model
##' 
##' SurvFrame <- data.frame(time=1:10,status=c(0,1,1,0,0,1,0,0,1,0))
##' SurvHist <- with(SurvFrame,Hist(time,status))
##' plot(SurvHist)
##' plot(SurvHist,box2.col=2,box2.label="experienced\nR user")
##' plot(SurvHist,
##'      box2.col=2,
##'      box1.label="newby",
##'      box2.label="experienced\nR user",
##'      oneFitsAll=FALSE,
##'      arrow1.length=.5,
##'      arrow1.label="",
##'      arrow1.lwd=4)
##' 
##' ## change the cex of all box labels:
##' plot(SurvHist,
##'      box2.col=2,
##'      box1.label="newby",
##'      box2.label="experienced\nR user",
##'      oneFitsAll=FALSE,
##'      arrow1.length=.5,
##'      arrow1.label="",
##'      arrow1.lwd=4,
##'      label.cex=1)
##' 
##' ## change the cex of single box labels:
##' plot(SurvHist,
##'      box2.col=2,
##'      box1.label="newby",
##'      box2.label="experienced\nR user",
##'      oneFitsAll=FALSE,
##'      arrow1.length=.5,
##'      arrow1.label="",
##'      arrow1.lwd=4,
##'      label1.cex=1,
##'      label2.cex=2)
##' 
##' 
##' ## The pbc data set from the survival package
##' library(survival)
##' data(pbc)
##' plot(with(pbc,Hist(time,status)),
##'      stateLabels=c("randomized","transplant","dead"),
##'      arrowLabelStyle="count")
##' 
##' ## two competing risks
##' comprisk.model <- data.frame(time=1:3,status=1:3)
##' CRHist <- with(comprisk.model,Hist(time,status,cens.code=2))
##' plot(CRHist)
##' plot(CRHist,arrow1.label=paste(expression(eta(s,u))))
##' 
##' plot(CRHist,box2.label="This\nis\nstate 2",arrow1.label=paste(expression(gamma[1](t))))
##' plot(CRHist,box3.label="Any\nLabel",arrow2.label="any\nlabel")
##' 
##' ## change the layout
##' plot(CRHist,
##'      box1.label="Alive",
##'      box2.label="Dead\n cause 1",
##'      box3.label="Dead\n cause 2",
##'      arrow1.label=paste(expression(gamma[1](t))),
##'      arrow2.label=paste(expression(eta[2](t))),
##'      box1.col=2,
##'      box2.col=3,
##'      box3.col=4,
##'      nrow=2,
##'      ncol=3,
##'      box1.row=1,
##'      box1.column=2,
##'      box2.row=2,
##'      box2.column=1,
##'      box3.row=2,
##'      box3.column=3)
##' 
##' ## more competing risks
##' comprisk.model2 <- data.frame(time=1:4,status=1:4)
##' CRHist2 <- with(comprisk.model2,Hist(time,status,cens.code=2))
##' plot(CRHist2,box1.row=2)
##' 
##' ## illness-death models
##' illness.death.frame <- data.frame(time=1:4,
##' 				  from=c("Disease\nfree",
##'                                       "Disease\nfree",
##'                                       "Diseased",
##'                                       "Disease\nfree"),
##' 				  to=c("0","Diseased","Dead","Dead"))
##' IDHist <- with(illness.death.frame,Hist(time,event=list(from,to)))
##' plot(IDHist)
##' 
##' ## illness-death with recovery
##' illness.death.frame2 <- data.frame(time=1:5,
##' from=c("Disease\nfree","Disease\nfree","Diseased","Diseased","Disease\nfree"),
##' to=c("0","Diseased","Disease\nfree","Dead","Dead"))
##' IDHist2 <- with(illness.death.frame2,Hist(time,event=list(from,to)))
##' plot(IDHist2)
##' 
##' ## 4 state models
##' x=data.frame(from=c(1,2,1,3,4),to=c(2,1,3,4,1),time=1:5)
##' y=with(x,Hist(time=time,event=list(from=from,to=to)))
##' plot(y)
##' 
##' ## moving the label of some arrows
##' 
##' d <- data.frame(time=1:5,from=c(1,1,1,2,2),to=c(2,3,4,3,4))
##' h <- with(d,Hist(time,event=list(from,to)))
##' plot(h,
##' tagBoxes=TRUE,
##' stateLabels=c("Remission\nwithout\nGvHD",
##'     "Remission\nwith\nGvHD",
##'     "Relapse",
##'     "Death\nwithout\nrelapse"),
##' arrowLabelSymbol='alpha',
##' arrowlabel3.x=35,
##' arrowlabel3.y=53,
##' arrowlabel4.y=54,
##' arrowlabel4.x=68)
##' 
##' ##'
#' @export 
plot.Hist <- function(x,
                      nrow,
                      ncol,
                      stateLabels,
                      arrowLabels,
                      arrowLabelStyle="symbolic",
                      arrowLabelSymbol='lambda',
                      changeArrowLabelSide,
                      tagBoxes=FALSE,
                      startCountZero=TRUE,                      
                      oneFitsAll,
                      margin,
                      cex,
                      verbose=FALSE,
                      ...){
    # {{{ margin 
    oldmar <- par()$mar
    oldoma <- par()$oma
    par(oma=c(0,0,0,0))
    oldxpd <- par()$xpd
    if (!missing(margin)){
        par(mar=rep(margin,length.out=4),xpd=TRUE)
    }
    else
        par(mar=c(0,0,0,0),xpd=TRUE)
    # }}}
    # {{{ find states 
    model.type <- attr(x,"model")
    states <- attr(x,"states")
    origStates <- states
    if (model.type!="multi.states"){ ## need an initial state
        states <- c("initial", states)
    }
    NS <- length(states)
    if (missing(stateLabels)){
        if (all(as.character(as.numeric(as.factor(origStates)))==origStates))  ## make nice state boxlabels if states are integers
            stateLabs <- switch(model.type,"survival"=paste(c("","Event"),states),"competing.risks"=paste(c("",rep("Cause",NS-1)),states),paste("State",states))
        else
            stateLabs <- states
    }
    else{
        if(length(stateLabels)==NS-1){
            stateLabs <- c("initial",stateLabels)
        }
        else{
            if (length(stateLabels)==NS){
                stateLabs <- stateLabels
            }
            else{
                stop("Wrong number of state names.")
            }
        }
    }
    ## forcedLabels
    thecall <- match.call(expand.dots=TRUE)
    labelhits <- match(paste("box",1:NS,".label",sep=""),names(thecall),nomatch=0)
    for (i in 1:NS){
        if (labelhits[i]!=0)
            ## may be language: thecall[[labelhits[i]]]
            ## if user specifies box2.label=c("Event 1")
            ## instead of box2.label="Event 1"
            stateLabs[i] <- eval(thecall[[labelhits[i]]])[1]
    }
    numstates <- as.numeric(as.character(factor(states,levels=states,labels=1:NS)))
    startCountZero <- TRUE 
    if (startCountZero)
        numstateLabels <- numstates-1
    else
        numstateLabels <- numstates
    # {{{  find transitions between the states

  ## first remove the censored lines from the transition matrix
  ## x <- x[x[,"status"]!=attr(x,"cens.code"),,drop=FALSE]
  x <- x[x[,"status"]!=0,,drop=FALSE]
  if (NROW(x)==0) stop("No uncensored transitions.")
  sumx <- summary(x,verbose=verbose)
  notCensored <- sumx$trans.frame$to!="unknown"
  sumx$trans.frame <- sumx$trans.frame[notCensored,]
  sumx$transitions <- sumx$transitions[notCensored]
  transitions <- sumx$trans.frame
  ordered.transitions <- unique(transitions)
  N <- NROW(ordered.transitions)
  # }}}

  # }}}
    # {{{ default layout: arranging the boxes

    state.types <- sumx$states
    state.types <- state.types[state.types>0]
    if (missing(nrow))
        if (model.type=="multi.states")
            nrow <- NS
        else
            if (ceiling(NS/2)==floor(NS/2))
                nrow <- NS-1
            else
                nrow <- NS
    if (missing(ncol))
        if (model.type=="multi.states")
            ncol <- NS
        else
            ncol <- 2
    ## placing boxes in rows and columns
    if (model.type=="multi.states"){
        adjustRowsInColumn <- rep(0,ncol)
        adjustColsInRow <- rep(0,nrow)
        box.col <- switch(as.character(NS),
                          "2"=c(1,ncol),
                          "3"=c(1,2,ncol),
                          "4"=c(1,1,ncol,ncol),
                          "5"=c(1,1,ceiling((ncol-1)/2),ncol,ncol),
                          "6"=c(1,3,3,5,6,6))
        box.row <- switch(as.character(NS),
                          "2"=c(1,1),
                          "3"=c(nrow,1,nrow),
                          "4"=c(1,nrow,1,nrow),
                          "5"=c(1,nrow,ceiling(nrow/2),1,nrow),
                          "6"=c(3,1,6,4,1,6))
    }
    else{ # survival or competing risks
        ## adjustRowsInColumn <- rep(1,ncol)
        ## adjustColsInRow <- rep(1,nrow)
        if (ceiling(NS/2)==floor(NS/2)){ ## equal number of states and unequal number of absorbing states
            box.col <- c(1,rep(ncol,NS-1))
            box.row <- c(NS/2,1:(NS-1))
        } else{ 
            box.col <- c(1,rep(ncol,NS-1))
            box.row <- c((NS+1)/2,(1:NS)[-(NS+1)/2])
        }
    }
    if (is.null(box.row) || is.null(box.col))
        stop("Please specify the layout for this ",NS," state model (")
    layoutDefaults <- data.frame(name=paste("box",1:NS,sep=""),
                                 row=box.row,
                                 column=box.col,
                                 stringsAsFactors=FALSE)
    layoutDefaultList <- lapply(1:NS,function(x)layoutDefaults[x,-1,drop=FALSE])
    names(layoutDefaultList) <- layoutDefaults$name
    layout <- SmartControl(list(...),
                           keys=c(layoutDefaults$name),
                           defaults=c(layoutDefaultList),
                           ignore.case=TRUE,
                           replaceDefaults=FALSE,
                           verbose=FALSE)

  # }}}
    # {{{ draw empty frame

  # plot
  Xlim <- 100
  Ylim <- 100
  plot(0,0,type="n",xlim=c(0,Xlim),ylim=c(0,Ylim),xlab="",ylab="",axes=FALSE)
  ## backGround(c(0,100),c(0,100),bg="yellow")

  # }}}
    # {{{ default values

  if (missing(cex))
    theCex <- 2
  else
    theCex <- cex
  if (found <- match("arrowLabel.cex",names(thecall),nomatch=0))
    arrowLabel.cex <- thecall[[found]]
  else
    arrowLabel.cex <- rep(theCex,N)
  ## boxes
  boxDefaults <- data.frame(name=paste("box",1:NS,sep=""),xpd=TRUE,stringsAsFactors=FALSE)
  ## box labels
  boxLabelDefaults <- data.frame(name=paste("label",1:NS,sep=""),stringsAsFactors=FALSE,label=stateLabs)
  ## arrows
  arrowDefaults <- data.frame(name=paste("arrow",1:N,sep=""),code=2,lwd=1,headoffset=strwidth("ab",cex=arrowLabel.cex),length=.13,stringsAsFactors=FALSE)
  arrowDefaults <- cbind(arrowDefaults,ordered.transitions)
  ## arrowlabels
  if (missing(changeArrowLabelSide))
    changeArrowLabelSide <- rep(FALSE,N)
  arrowlabelDefaults <- data.frame(name=paste("arrowlabel",1:N,sep=""),
                                   label=arrowLabelStyle,
                                   x=NA,
                                   y=NA,
                                   stringsAsFactors=FALSE,
                                   cex=arrowLabel.cex)
  arrowlabelDefaults <- cbind(arrowlabelDefaults,ordered.transitions)
  arrowlabelDefaults$numfrom <- factor(arrowlabelDefaults$from,levels=states,labels=numstateLabels)
  arrowlabelDefaults$numto <- factor(arrowlabelDefaults$to,levels=states,labels=numstateLabels)
  if (missing(arrowLabels)){
    arrowLabels <- NULL
  }
  arrowLabels.p <- TRUE
    if (length(arrowLabels)>0 &&is.logical(arrowLabels) && arrowLabels==FALSE){
    arrowLabels <- rep("",N)
    arrowLabels.p <- FALSE
  }
  else{
    if (length(arrowLabels)==0){
      arrowLabels <- lapply(1:N,function(i){
        bquote(paste(expression(.(as.name(arrowLabelSymbol))[.(paste(as.character(arrowlabelDefaults$numfrom[i]),
            as.character(arrowlabelDefaults$numto[i]),
            sep=""))](t))))
      })
    } else{
      stopifnot(length(arrowLabels)==N)
    }
  }
  arrowlabelhits <- match(paste("arrow",1:N,".label",sep=""),names(thecall),nomatch=0)
  for (i in 1:N){
    if (arrowlabelhits[i]!=0){
      arrowLabels[[i]] <- thecall[[arrowlabelhits[i]]]
    }
  }

  # }}}
    # {{{ compute box dimensions relative to cex of box labels

  ## to find the cex for the box labels, first initialize
  boxLabelCex <- rep(theCex,NS)
  ## then look for label.cex
  if (theLabelCex <- match("label.cex",names(thecall),nomatch=0)){
    boxLabelCex <- rep(thecall[[theLabelCex]],NS)
  }
  # finally adjust for box individual values 
  if (any(iLabelCex <- match(paste("label",1:NS,".cex",sep=""),names(thecall),nomatch=0))){
    for (i in 1:NS){
      if ((argi <- iLabelCex[i])!=0)
        boxLabelCex[i] <- thecall[[argi]]
    }
  }
  ## state.cex <- max(boxLabelCex)
  if (length(boxLabelCex)<length(stateLabs))
    boxLabelCex <- rep(boxLabelCex,length.out=length(stateLabs))
  state.width <- sapply(1:length(stateLabs),function(i){strwidth(stateLabs[i],cex=boxLabelCex[i])})
  state.height <- sapply(1:length(stateLabs),function(i){strheight(stateLabs[i],cex=boxLabelCex[i])})
  ## state.width <- sapply(stateLabs,strwidth,cex=boxLabelCex)
  ## state.height <- sapply(stateLabs,strheight,cex=boxLabelCex)

  if (missing(oneFitsAll))
    oneFitsAll <- length(unique(boxLabelCex))==1
  if (oneFitsAll==TRUE){
    max.width <- max(state.width)
    max.height <- max(state.height)
    ##     box.width <- max.width + xbox.rule * max.width
    ##     box.height <- max.height + ybox.rule * max.height
    box.width <- max.width + strwidth("ab",cex=max(boxLabelCex))
    box.height <- max.height + strwidth("ab",cex=max(boxLabelCex))
    ## really need to check this for each row:
    ##     if ((ncol * box.width) > Xlim) warning("The horizontal dimensions of the boxes are too big -- change layout or tune parameters `label.cex' and/or `xbox.rule'.")
    ##     if ((nrow * box.height) > Ylim) warning("The verticalf dimensions of the boxes are too big -- change layout or tune parameters `label.cex' and/or `ybox.rule'.")
  }
  else{
    box.width <- state.width + strwidth("ab",cex=boxLabelCex)
    box.height <- state.height + strwidth("ab",cex=boxLabelCex)
  }

  if (length(box.height)==1) box.height <- rep(box.height,NS)
  if (length(box.width)==1) box.width <- rep(box.width,NS)
  # }}}
    # {{{ arrange the boxes in the layout

  boxCol <- sapply(layout,function(x){x$column})
  if (any(boxCol>ncol)) ncol <- max(boxCol)
  boxRow <- sapply(layout,function(x){x$row})
  if (any(boxRow>ncol)) nrow <- max(boxRow)
  ybox.position <- numeric(NS)
  names(ybox.position) <- paste("box",numstates,sep="")
  # {{{y box positions
  for (x in 1:ncol){
      ## For each column find y positions for boxes
      boxesInColumn <- names(boxCol)[boxCol==x]
      boxesInColumnNumbers <- as.numeric(sapply(strsplit(boxesInColumn,"box"),function(x)x[[2]]))
      if (length(boxesInColumn)>0){
          ## if (adjustRowsInColumn[x]==1 && all(match(paste(boxesInColumn,"row",sep="."),names(thecall),nomatch=0)==0)){
          # adjust the y position of the boxes according to the number of boxes in column
          ## yPossible <- centerBoxes(Ylim,box.height[boxesInColumnNumbers],nrow,boxRow[boxesInColumn])
          ## for (b in 1:length(boxesInColumn))
          ## ybox.position[boxesInColumn[b]] <- yPossible[b]
          ## }
          ## else{
          yPossible <- centerBoxes(Ylim,box.height[boxesInColumnNumbers],nrow,boxRow[boxesInColumn])
          for (b in 1:length(boxesInColumn)){
              ybox.position[boxesInColumn[b]] <- yPossible[b]
              ## }
          }
      }
  }
  ## row 1 is on top but the y-axis starts at the button
  ## therefore need to transform
  ybox.position <- 100-(ybox.position+box.height)
  # }}}
  # {{{x box positions
  xbox.position <- numeric(NS)
  names(xbox.position) <- paste("box",numstates,sep="")
  for (x in 1:nrow){
      ## For each row find x positions for boxes
      boxesInRow <- names(boxRow)[boxRow==x]
      boxesInRowNumbers <- as.numeric(sapply(strsplit(boxesInRow,"box"),function(x)x[[2]]))
      if (length(boxesInRow)>0){
          ## if (adjustColsInRow[x]==1 && all(match(paste(boxesInRow,"row",sep="."),names(thecall),nomatch=0)==0)){
          # adjust the x position of the boxes according to the number of boxes in row
          ## xpossible <- centerBoxes(Ylim,box.height[boxesInRowNumbers],ncol,boxCol[boxesInRow])
          ## for (b in 1:length(boxesInRow))
          ## xbox.position[boxesInRow[b]] <- xpossible[b]
          ## }
          ## else{
          if (sum(box.width[boxesInRowNumbers])>Xlim)
              stop(paste("Sum of box widths in row",x,"exceed limit",Xlim))
          xpossible <- centerBoxes(Xlim,box.width[boxesInRowNumbers],ncol,boxCol[boxesInRow])
          ## if (any(xpossible<0)) browser()
          for (b in 1:length(boxesInRow)){
              xbox.position[boxesInRow[b]] <- xpossible[b]
          }
          ## }
      }
  }
  # }}}
  xtext.position <- xbox.position + (box.width - state.width)/2
  ytext.position <- ybox.position + (box.height - state.height)/2
  if (verbose){
      cat("\n\nBoxlabel data:\n\n")
      print(data.frame(stateLabs,
                       boxCol,
                       boxRow,
                       x.pos=round(xbox.position,2),
                       y.pos=round(ybox.position,2),
                       width=round(box.width,2),
                       label.width=round(state.width,2),
                       label.height=round(state.height,2),
                       boxLabelCex))
  }
  boxDefaults <- cbind(boxDefaults,xleft=xbox.position,ybottom=ybox.position,xright=xbox.position+box.width,ytop=ybox.position+box.height)
  boxLabelDefaults <- cbind(boxLabelDefaults,
                            x=xtext.position,
                            y=ytext.position,
                            cex=boxLabelCex)

  # }}}
    # {{{ compute arrow positions

  doubleArrow <- match(paste(arrowDefaults[,"to"],arrowDefaults[,"from"]),paste(arrowDefaults[,"from"],arrowDefaults[,"to"]),nomatch=0)
  arrowDefaults <- cbind(arrowDefaults,doubleArrow)
  arrowList <- for (trans in 1:N){
    from.state <- factor(ordered.transitions[trans,1],levels=states,labels=numstates)
    to.state <- factor(ordered.transitions[trans,2],levels=states,labels=numstates)
    ArrowPositions <- findArrow(Box1=c(round(xbox.position[from.state],4),round(ybox.position[from.state],4)),
                                Box2=c(round(xbox.position[to.state],4),round(ybox.position[to.state],4)),
                                Box1Dim=c(box.width[from.state],box.height[from.state]),
                                Box2Dim=c(box.width[to.state],box.height[to.state]),
                                verbose=FALSE)
    Len <- function(x){sqrt(sum(x^2))}
    from <- ArrowPositions$from
    to <- ArrowPositions$to
    ArrowDirection <- to-from
    ArrowDirection <- ArrowDirection/Len(ArrowDirection)
    ## perpendicular direction
    PerDir <- rev(ArrowDirection)*c(1,-1)/Len(ArrowDirection)
    ## shift double arrows
    dd <- arrowDefaults[trans,"doubleArrow"]
    if (dd!=0){
      dist <- strwidth(".",cex=arrowLabel.cex)
      arrowDefaults[trans,"headoffset"]+dist
      if (dd>trans){
        from <- from + sign(PerDir) * c(dist,dist)
        to <- to + sign(PerDir) * c(dist,dist)
      }
      else{
        from <- from + sign(PerDir) * c(dist,dist)
        to <- to + sign(PerDir) * c(dist,dist)
      }
    }
    # shift the start and end points of arrows by ArrowHeadOffset
    ArrowHeadOffset <- arrowDefaults[trans,"headoffset"]
    from <- from+sign(ArrowDirection)*c(ArrowHeadOffset,ArrowHeadOffset)*abs(ArrowDirection)
    to <- to-sign(ArrowDirection)*c(ArrowHeadOffset,ArrowHeadOffset)*abs(ArrowDirection)
    arrowDefaults[trans,"x0"] <- from[1]
    arrowDefaults[trans,"x1"] <- to[1]
    arrowDefaults[trans,"y0"] <- from[2]
    arrowDefaults[trans,"y1"] <- to[2]
    ## shift arrow label perpendicular (left) to arrow direction
    offset <- strwidth(".",cex=arrowLabel.cex)
    ArrowMid <- (to+from)/2
    ## points(x=ArrowMid[1],y=ArrowMid[2],col=3,pch=16)
    if (changeArrowLabelSide[trans]==TRUE)
    ArrowLabelPos <- ArrowMid - sign(PerDir) * c(offset,offset)
    else
    ArrowLabelPos <- ArrowMid + sign(PerDir) * c(offset,offset)
    try1 <- try(mode((arrowLabels[[trans]])[2])[[1]]=="call",silent=TRUE)
    ## try2 <- try(as.character(arrowLabels[[trans]])[[1]]=="paste",silent=TRUE)
    labIsCall <- (class(try1)!="try-error" && try1)
    ## labUsePaste <- (class(try2)!="try-error" && try2)
    if (labIsCall){ # symbolic label
      arrowLabels[[trans]] <- ((arrowLabels[[trans]])[2])[[1]][[2]]
    }
    ## relative label height
    lab <- arrowLabels[[trans]]
    labelHeight <- strheight(lab,cex=arrowlabelDefaults[trans,"cex"])
    ## relative label width 
    labelWidth <-  strwidth(lab,cex=arrowlabelDefaults[trans,"cex"])
    ## shift further according to label height and width in perpendicular direction
    if (changeArrowLabelSide[trans]==TRUE)
      ArrowLabelPos <- ArrowLabelPos-sign(PerDir)*c(labelWidth/2,labelHeight/2)
    else
    ArrowLabelPos <- ArrowLabelPos+sign(PerDir)*c(labelWidth/2,labelHeight/2)
    arrowlabelDefaults[trans,"x"] <- ArrowLabelPos[1] 
    arrowlabelDefaults[trans,"y"] <- ArrowLabelPos[2]
  }

  # }}}
    # {{{ Smart argument control

  boxDefaultList <- lapply(1:NS,function(x)boxDefaults[x,-1,drop=FALSE])
  names(boxDefaultList) <- boxDefaults$name
  boxLabelDefaultList <- lapply(1:NS,function(x)boxLabelDefaults[x,-1,drop=FALSE])
  names(boxLabelDefaultList) <- boxLabelDefaults$name
  arrowDefaultList <- lapply(1:N,function(x)arrowDefaults[x,-1,drop=FALSE])
  names(arrowDefaultList) <- as.character(arrowDefaults$name)
  arrowlabelDefaultList <- lapply(1:N,function(x)arrowlabelDefaults[x,-1,drop=FALSE])
  names(arrowlabelDefaultList) <- as.character(arrowlabelDefaults$name)
  boxTagsDefaultList <- list(labels=numstateLabels,cex=1.28,adj=c(-.5,1.43))
  smartArgs <- SmartControl(list(...),
                            keys=c(boxDefaults$name,
                                boxLabelDefaults$name,
                                as.character(arrowDefaults$name),
                                as.character(arrowlabelDefaults$name),
                                "boxtags"),
                            defaults=c(boxLabelDefaultList,arrowDefaultList,arrowlabelDefaultList,boxDefaultList,list("boxtags"=boxTagsDefaultList)),
                            ignore.case=TRUE,
                            replaceDefaults=FALSE,
                            verbose=verbose)
    
  # }}}
    # {{{  draw the boxes
  
  for (i in 1:NS) {
    suppressWarnings(do.call("rect",smartArgs[[paste("box",i,sep="")]]))
  }

  # }}}
    # {{{  label the boxes
  
  for (i in 1:NS) {
    suppressWarnings(do.call("text",c(list(adj=c(0,0)),smartArgs[[paste("label",i,sep="")]])))
  }

  # }}}
    # {{{  draw the arrows

  for (i in 1:N){
    suppressWarnings(do.call("arrows",c(smartArgs[[paste("arrow",i,sep="")]])))
  }

  # }}}
    # {{{ label the arrows
    if (verbose) arrowLabel.data <- NULL
    if (arrowLabels.p==TRUE){
        for (i in 1:N){
            labelList <- smartArgs[[paste("arrowlabel",i,sep="")]]
            if (verbose) arrowLabel.data <- rbind(arrowLabel.data,cbind("arrowLabel"=i,data.frame(labelList)))
            switch(labelList$label,"symbolic"={
                ## lab <- (arrowLabels[[i]])
                try1 <- try(mode((arrowLabels[[i]])[2])[[1]]=="call",silent=TRUE)
                ## try2 <- try(as.character(arrowLabels[[i]])[[1]]=="paste",silent=TRUE)
                labIsCall <- (class(try1)!="try-error" && try1)
                suppressWarnings(do.call("text",c(list(labels=bquote(arrowLabels[[i]])),labelList)))        
            }, "count"={
                tabTrans <- as.matrix(table(transitions))
                lab <- paste("n=",tabTrans[as.character(labelList$from),as.character(labelList$to)])
                suppressWarnings(do.call("text",c(list(labels=quote(lab)),labelList)))
            })
            ## suppressWarnings(do.call("text",c(list(adj=c(labelWidth/2,labelHeight/2),labels="label"),smartArgs[[paste("arrowlabel",i,sep="")]])))
        }
    }
    if (verbose) {
        cat("\n\nArrow label data:\n\n")
        print(arrowLabel.data)
    }
    # }}}
    # {{{  put numbers in the upper left corner of the boxes (if wanted)

  if (tagBoxes==TRUE){
    tagList <- smartArgs$boxtags
    nix <- lapply(1:NS,function(b) {
      lab <- tagList[b]
      text(x=xbox.position[b],
           y=ybox.position[b]+box.height,
           labels=tagList$labels[b],
           cex=tagList$cex,
           adj=tagList$adj)})
  }

  # }}}
    # {{{ reset margin
  par(mar=oldmar,xpd=oldxpd,oma=oldoma)
  # }}}
    if (verbose){
        cat("\nRelevel the factor 'event' in the dataset which defines the Hist object,\nto change the order of the boxes.\n")
    }
    invisible(smartArgs)
}


position.finder <- function(border,len,n){
## distribute the boxes of lenght len uniformly
## over [0,border]
 if (n==1)
    (border - len)/2
  else{
    seq(0,border-.5*len,len + (border-(n * len))/(n-1))
  }  
}

centerBoxes <- function(border,len,ncell,pos){
    ## box i has length len[i] and is centered in cell pos[i]
    ## return the position in [0,border] of the lower
    ## border of the boxes
    cellwidth <- border/ncell
    nboxes <- length(len)
    if ((luft <- border-sum(len))<0) stop("sum of box dimensions exceeds limits")
    if (nboxes>ncell) stop("too many boxes in one row")
    ## case: all boxes fit into given cell width
    ## if (all(len<cellwidth)){
    box.pos <- seq(from=0,to=border,by=cellwidth)[pos] + pmax(0,sapply(len,function(l) {(cellwidth - l)/2}))
    ## spread as far as possible
    boxPos <- sapply(1:length(box.pos),function(b){
        bp <- box.pos[b]
        if (ncell>1 && pos[b]==1) # at the left/lower border
            bp <- min(0,abs(box.pos[b]))
        if (ncell> 1 && pos[b]==ncell)# at the right/upper border
            bp <- max(border-len[b],box.pos[b])
        bp
    })
    ## }else{
    ## ## case: at least one box exceeds the cellwidth
    ## between <- luft/(nboxes-1)
    ## boxPos <- c(0,len[-nboxes]+between)
    ## }
    boxPos
}


## positionFinder <- function(border,len,n){
## distribute the whitespace between the boxes
## instead of the boxes
## wspace <- border-sum(len)
## if (n==1)
## (border - len)/2
## else{
## seq(0,border-.5*len,len + (border-(n * len))/(n-1))
## }  
## }
