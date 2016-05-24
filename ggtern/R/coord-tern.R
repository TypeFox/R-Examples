#' Ternary Coordinate System
#' 
#' \code{coord_tern} is a function which creates a transformation mechanism between the ternary system, and, the cartesian system.
#' It inherits from the fixed coordinate system, employing fixed ratio between x and y axes once transformed.
#' 
#' @section Aesthetics (Required in Each Layer):
#' \Sexpr[results=rd,stage=build]{ggtern:::rd_aesthetics("coord", "tern")}
#' 
#' Abovementioned limitations include the types of geometries which can be used (ie approved geometries), 
#' or modifications to required aesthetic mappings. One such essential patch is, for approved geometries previously 
#' requiring \code{x} and \code{y} coordinates, now require an additional \code{z} coordinate, and, 
#' \code{\link{geom_segment}} goes one step further in that it requires both an additional 
#' \code{z} and \code{zend} coordinate mappings. 
#' 
#' In essence, the required aesthetics are the product between what
#' is required of each 'layer' and what is required of the 'coordinate system'.
#' @param Tlim the range of T in the ternary space
#' @param Llim the range of L in the ternary space
#' @param Rlim the range of R in the ternary space
#' @inheritParams ggplot2:::coord_cartesian
#' @return \code{coord_tern} returns a CoordTern ggproto
#' @rdname coord-tern
#' @author Nicholas Hamilton
#' @export
coord_tern <- function(Tlim   = NULL, Llim   = NULL, Rlim   = NULL,expand=TRUE){
  all.coords <- c("x","y","z")
  all.scales <- c("T","L","R")
  mapping    <- sapply(all.scales,function(x){ getOption(sprintf('tern.default.%s',x))}  )
  if(length(intersect(all.coords,as.character(mapping))) != 3)
    stop("Options for T, L and R are x,y and z, must be assigned and NOT duplicated, alter the 'tern.default.X' (X=T,L or R) global option",call.=FALSE)
  ggproto(NULL, CoordTern,
          required_aes    = all.coords,
          required_scales = all.scales,
          mapping         = as.list(mapping),
          limits          = list(x = c(0,1), y = c(0,1)*.ratio(), T=Tlim, L=Llim, R=Rlim),
          ratio           = 1,
          expand          = expand,
          scales          = list(),
          labels_coord    = list(),
          theme           = theme_get()
  )
}

#' @rdname coord-tern
#' @export
CoordTern <- ggproto("CoordTern", CoordCartesian,
  aspect = function(self, ranges) { (diff(ranges$y.range) / diff(ranges$x.range)) * self$ratio },
  transform = function(self, data, scale_details){
    
    #Variables and functions.
    ix               = c('x','y')
    angle            = .theme.get.rotation(self)
    addOrigin        = function(d,ix,o){  d[,ix] = t(t(d[,ix,drop=FALSE]) + o); d }
    data             = tlr2xy(data,self)
    coordShift       = sapply(c('tern.axis.hshift','tern.axis.vshift'),function(x){ calc_element(x,self$theme) %||% 0})
    
    #For each coordinate group, conduct the re-centering, translation / rotation process
    for(group in .get.sets(ix,names(data)) ){
      ix.comb        = .combos(ix,group)
      
      #Run check to see if all the required columns are present within the group.
      #A required column is the cartesian product of the group and the required index, so it could be (xend,yend), (xstart,ystart), etc...
      #This shouldn't be a problem, since the same check is executed during the master transformation in the tlr2xy(...) operation above
      if(!all(ix.comb %in% names(data))) stop(sprintf("Problem with proposed tranlation/rotation, require columns (%s) wich are not present.",
                                                      joinCharacterSeries(setdiff(ix.comb,names(data)),lastWord='and') ))
      
      #Determine Origin
      xtrm           = tlr2xy(.get.tern.extremes(self,scale_details,FALSE),self)
      origin         = apply(xtrm[,ix],2,function(x){mean(x)}) ##ORIGIN TO BE CENTROID
      
      #Any multiple of 360 degrees can be disregarded
      if(angle %% 360 != 0){
        #Translation   
        data           = addOrigin(data,ix.comb,-origin); 
        xtrm           = addOrigin(xtrm,ix,     -origin)
        #Rotation  
        data[,ix.comb] = .rotation(data[,ix.comb],angle) 
        xtrm[,ix]      = .rotation(xtrm[,ix],     angle)
        #Translate Back
        data           = addOrigin(data,ix.comb,+origin)
        xtrm           = addOrigin(xtrm,ix,     +origin)
      }
      
      #Re-Center
      origin         = apply(xtrm[,ix],2,function(x){ mean(range(x))} ) ##ORIGIN TO BE MEDIAN
      target         = c(mean(scale_details$x.range),mean(scale_details$y.range)) + coordShift
      data           = addOrigin(data,ix.comb,(target - origin))
    }
    
    self$super$super$transform(data,scale_details)
  },
  render_axis_h = function(self,scale_details, theme) zeroGrob(), #not required
  render_axis_v = function(self,scale_details, theme) zeroGrob(), #not required
  render_fg     = function(self,scale_details, theme){
    items = list()
    extrm = .get.tern.extremes(self,scale_details)
    if(.theme.get.gridsontop(theme))
      items = .render.grids(self,extrm,scale_details,theme,items,onBackground = TRUE)
    else
      items = .render.grids(self,extrm,scale_details,theme,items,onBackground = FALSE)
    if(.theme.get.bordersontop(theme))
      items = .render.borders(self,extrm,theme,items)
    items = .render.arrows(self,extrm,scale_details,theme,items)
    items = .render.titles(self,extrm,scale_details,theme,items)
    gTree(children = do.call("gList", items))
  },
  render_bg     = function(self,scale_details, theme){
    items = list()
    extrm = .get.tern.extremes(self,scale_details)
    items = .render.background(self,extrm,theme,items)
    if(!.theme.get.gridsontop(theme))
      items = .render.grids(self,extrm,scale_details,theme,items,onBackground = TRUE)
    if(!.theme.get.bordersontop(theme))
      items = .render.borders(self,extrm,theme,items)
    gTree(children = do.call("gList",items))
  },
  train = function(self, scale_details) {
    train_cartesian <- function(scale_details,limits,name,continuousAmount) {
      if (self$expand) {
        expand <- ggint$expand_default(scale_details,continuous=continuousAmount)
      } else {
        expand <- c(0, 0)
      }
      if (is.null(limits)) {
        range <- scale_details$dimension(expand)
      } else {
        range <- range(scale_details$transform(limits))
        range <- expand_range(range, expand[1], expand[2])
      }
      
      out <- scale_details$break_info(range)
      names(out) <- paste(name, names(out), sep = ".")
      out
    }
    
    #Determine the epansiion amount
    expand.amount = calc_element('tern.panel.expand',theme=self$theme)
    
    #Adjust for rotation
    extremes        = .get.tern.extremes(self,list(x.range=self$limits$x,
                                                   y.range=self$limits$y))[,c('x','y')]
    #print(extremes)
    #print(self$limits)
    currentMidpoint = c(mean(self$limits$x),mean(self$limits$y))
    ternaryMidpoint = apply(extremes,2,function(x){mean(range(x))})
    shift           = 0*(currentMidpoint - ternaryMidpoint)
    
    #Execute the Training
    ret = c(
      train_cartesian(scale_details$x, self$limits$x - shift[1],"x",c(expand.amount,0) ),
      train_cartesian(scale_details$y, self$limits$y - shift[2],"y",c(expand.amount,0) )
    )
    
    #Done
    ret
  }
)

#----------------------------------------------------------------------------------
#Internals >>>> Ratio
#----------------------------------------------------------------------------------
.ratio = function(){ 0.5*tan(60*pi/180) }

#----------------------------------------------------------------------------------
#Internals >>>> Ternary Extremes
#----------------------------------------------------------------------------------
.get.tern.extremes <- function(self,scale_details,transform=TRUE){
  ex = .check.tern.extremes(self)
  if(transform) ex = self$transform(ex,scale_details)
  ex
}

.check.tern.extremes <- function(self,precision=7){
  ix           = self$required_scales
  scales       = self$scales[ix]
  limitsfunc   = Vectorize(function(r,c){ do.call(if(r == c)'max'else'min',list( scales[[ ix[c] ]]$limits %||% self$limits[[ ix[c] ]] %||% c(0,1) ))  })
  ex           = as.data.frame(outer(1:3,1:3,limitsfunc))
  colnames(ex) = as.character(self$mapping[ix])
  sums         = round(apply(ex,1,sum),precision)
  if(!all(sums == 1.0)){
    colnames(ex) = ix; ex$Sum = sums; print(ex)
    stop("Invalid Ternary Limits, Each Point Must Sum to Unity...",call.=FALSE) 
  }
  invisible(ex)
}

#----------------------------------------------------------------------------------
#Internals >>>> ANGLES
#Functions to determine the rotation angles for the various components
#----------------------------------------------------------------------------------
.get.angles             <- function(clockwise){ c(0,120,240) - as.logical(clockwise)*180 }
.get.angles.arrows      <- function(clockwise){ .get.angles(clockwise) + if(clockwise){-30}else{30} }
.get.angles.arrowmarker <- function(clockwise){ x = c(60,0,-60); if(clockwise){x}else{x[c(3,1:2)]}}
.get.angles.ticklabels  <- function(clockwise){ c(0,-60,60) }


#----------------------------------------------------------------------------------
#Internals >>>> Theme flags.
#----------------------------------------------------------------------------------
.theme.get.clockwise <- function(theme){
  clockwise = calc_element('tern.axis.clockwise',theme)
  ifthenelse(is.logical(clockwise),clockwise[1],getOption("tern.clockwise"))
}
.theme.get.gridsontop <- function(theme){
  ret = calc_element('tern.panel.grid.ontop',theme)
  ifthenelse(is.logical(ret),ret[1],FALSE)
}

.theme.get.bordersontop <- function(theme){
  ret = calc_element('tern.axis.line.ontop',theme)
  ifthenelse(is.logical(ret),ret[1],FALSE)
}

.theme.get.showtitles <- function(theme){
  showtitles = calc_element("tern.axis.title.show",theme=theme)
  ifthenelse(is.logical(showtitles),showtitles[1],getOption("tern.showtitles"))
}
.theme.get.showlabels <- function(theme){
  showlabels = calc_element("tern.axis.text.show",theme=theme)
  ifthenelse(is.logical(showlabels),showlabels[1],getOption("tern.showlabels"))
}
.theme.get.showgrid.major <- function(theme){
  showgrid   = calc_element("tern.panel.grid.major.show",theme=theme)
  ifthenelse(is.logical(showgrid),showgrid[1],getOption("tern.showgrid.major"))
}
.theme.get.showgrid.minor <- function(theme){
  showgrid   = calc_element("tern.panel.grid.minor.show",theme=theme)
  ifthenelse(is.logical(showgrid),showgrid[1],getOption("tern.showgrid.minor"))
}
.theme.get.outside       <- function(theme){
  outside     = calc_element("tern.axis.ticks.outside",theme=theme)
  ifthenelse(is.logical(outside),outside[1],getOption("tern.ticks.outside"))
}
.theme.get.showprimary   <- function(theme){
  showprimary = calc_element("tern.axis.ticks.primary.show",theme=theme)
  ifthenelse(is.logical(showprimary), showprimary[1],getOption("tern.ticks.showprimary"))
}
.theme.get.showsecondary <- function(theme){
  showsecondary = calc_element("tern.axis.ticks.secondary.show",theme=theme)
  ifthenelse(is.logical(showsecondary),showsecondary[1],getOption("tern.ticks.showsecondary"))
}
.theme.get.showarrows    <- function(theme){
  showarrows   = calc_element('tern.axis.arrow.show',theme=theme)
  ifthenelse(is.logical(showarrows),showarrows[1],getOption("tern.showarrows"))
}
.theme.get.label <- function(self,n,d=n,suffix=''){
  x      = function(n,s) sprintf("%s%s",n,s)
  if(!is.character(n)) return('')
  labels = self$labels_coord
  ix     = unique(c(x(n,suffix),x(self$mapping[[n]],suffix),n,self$mapping[[n]])) 
  id     = which(ix %in% names(labels))
  if(length(id) == 0) d else labels[[ix[id[1]]]]
}
.theme.get.rotation <- function(self){
  tryCatch({
    angle = calc_element('tern.panel.rotate',self$theme)
    return(ifthenelse(is.finite(angle),angle,0)[1])
  },error=function(e){ warning(e) })
  return(0)
}

#----------------------------------------------------------------------------------
#Internals >>>> Render Components
# -Backgrounds
# -Borders
# -Grids
# -Precession Arrows & Markers
# _Apex Titles
#----------------------------------------------------------------------------------
.render.background <- function(self,data.extreme,theme,items){
  grob  <- zeroGrob()
  tryCatch({
      e     <- calc_element('tern.plot.background',theme=theme,verbose=F)
      if(!identical(e,element_blank())){
        grob  <- polygonGrob( data.extreme$x, data.extreme$y, 
                              default.units = "npc",
                              id   = rep(1,nrow(data.extreme)),
                              gp   = gpar(  col  = e$colour,
                                            fill = alpha(e$fill,ifthenelse(!is.numeric(e$alpha),1,e$alpha)),
                                            lwd  = ifthenelse(!is.numeric(e$size),0,e$size)*find_global_tern(".pt"),
                                            lty  = e$linetype
                              )
        )
      }
  },error = function(e){  warning(e) })
  items[[length(items) + 1]] <- grob
  items
}

.render.borders    <- function(self,data.extreme,theme,items){
  
  .render.border = function(name,s,f,items){
    grob = zeroGrob()
    tryCatch({
      e = calc_element(name,theme=theme,verbose=FALSE)
      if(identical(e,element_blank()))return(items)
      grob = segmentsGrob(
          x0 = data.extreme$x[s], x1 = data.extreme$x[f], y0 = data.extreme$y[s], y1 = data.extreme$y[f],
          default.units="npc",
          gp = gpar(col     = e$colour, 
                    lty     = e$linetype,
                    lineend = e$lineend,
                    lwd     = e$size*find_global_tern(".pt"))
      )
    },error=function(e){
      warning(e)
    })
    items[[length(items) + 1]] <- grob
    items
  }
  
  #Render the Border
  .render.mainborder <- function(items){
    tryCatch({
      e  = calc_element('tern.plot.background',theme,verbose=F)
      if(identical(e,element_blank())) return(items)
      grob     <- polygonGrob(  x = data.extreme$x,
                                y = data.extreme$y,
                                default.units = "npc",
                                id   = rep(1,nrow(data.extreme)),
                                gp   = gpar(  col  = e$colour, 
                                              fill = NA, 
                                              lwd  = is.numericor(e$size,0)*find_global_tern(".pt"), 
                                              lty  = e$linetype)
      )
      items[[length(items) + 1]] = grob
    },error=function(e){})
    items
  }
  
  #Process the main border
  items <- .render.mainborder(items)
  
  #process the axes
  if(.theme.get.clockwise(theme)){
    items <- .render.border("tern.axis.line.T",2,1,items)
    items <- .render.border("tern.axis.line.L",3,2,items)
    items <- .render.border("tern.axis.line.R",1,3,items)
  }else{
    items <- .render.border("tern.axis.line.T",3,1,items)
    items <- .render.border("tern.axis.line.L",1,2,items)
    items <- .render.border("tern.axis.line.R",2,3,items)
  }
  
  items
}

.render.grids      <- function(self,data.extreme,scale_details,theme,items,onBackground=TRUE){
  #Process the flags.
  clockwise     <- .theme.get.clockwise(theme)
  outside       <- .theme.get.outside(theme)
  showprimary   <- .theme.get.showprimary(theme)
  showsecondary <- .theme.get.showsecondary(theme)
  showgrid.major<- .theme.get.showgrid.major(theme)
  showgrid.minor<- .theme.get.showgrid.minor(theme)
  showlabels    <- .theme.get.showlabels(theme)
  
  #major & minor ticklength
  tl.major <- tl.minor <- 0
  tryCatch({
    tl.major <- convertUnit(theme$tern.axis.ticks.length.major,"npc",valueOnly=T)
    tl.minor <- convertUnit(theme$tern.axis.ticks.length.minor,"npc",valueOnly=T)
  },error=function(e){  warning(e) })
  
  #Top, Left Right sequence.
  seq.tlr <- c("T","L","R")
  
  #ASSEMBLE THE GRID DATA.
  .get.data <- function(X,ix,major=TRUE,angle=0,angle.text=0){
    existing = data.frame()
    tryCatch({
      scale = self$scales[[X]]
      
      #DETERMINE THE BREAKS
      breaks <- if(major){ scale$breaks }else{ scale$minor_breaks }
      
      #BYPASS IF NECESSARY
      if(length(breaks) == 0) return(existing)
      
      labels <- ifthenelse(major,scale$labels,"")
      labels <- as.character(ifthenelse(identical(labels,waiver()),100*breaks,labels))
      
      #Assign new id.
      id <- (max(existing$ID,0) + 1)
      limits <- is.numericor(scale$limits,c(0,1))
      ix     <- min(ix,ifthenelse(major,length(tl.major),length(tl.minor)))
      majmin <- ifthenelse(major,"major","minor")  #Major or Minor Element Name part.
      
      #The new dataframe
      new            <- data.frame(ID = id,Scale=X,breaks,Labels=labels,Major=major)
      new            <- subset(new,breaks >= min(limits) & breaks <= max(limits))
      new$Prop       <- (new$breaks - min(limits)) / abs(diff(limits))
      new$TickLength <- ifthenelse(major,tl.major[ix],tl.minor[ix])
      new$NameText   <- paste0("tern.axis.text.",X)
      new$NameTicks  <- paste0("tern.axis.ticks.",majmin,".",X)
      new$NameGrid   <- paste0("tern.panel.grid.",majmin,".",X)
      new$Major      <- major
      
      ##Start and finish positions of scale.
      out       <- c("x","y")
      
      #Start indexes.
      ix.s <- which(seq.tlr == X);
      finish <- as.numeric(data.extreme[ix.s,out])
      
      #For Ticks
      ix.f <- ifthenelse(clockwise,if(ix.s == 3){1}else{ix.s+1},if(ix.s == 1){3}else{ix.s-1})
      start  <- as.numeric(data.extreme[ix.f,out])
      for(i in 1:length(out))
        new[,out[i]] <- new$Prop*(finish[i]-start[i]) + start[i]
      
      #FOR GRID
      ix.f <- ifthenelse(clockwise,if(ix.s == 1){3}else{ix.s-1},if(ix.s == 3){1}else{ix.s+1})
      start  <- as.numeric(data.extreme[ix.f,out])
      for(i in 1:length(out))
        new[,paste0(out[i],"end.grid")] <- new$Prop*(finish[i]-start[i]) + start[i]
      
      #The tick angles.
      new$Angle      <- angle + .theme.get.rotation(self)
      new$Angle.Text <- .valid.angle(angle.text)
      
      #Determine the tick finish positions for segments.
      new$xend <- cos(new$Angle*pi/180)*new$TickLength + new$x
      new$yend <- sin(new$Angle*pi/180)*new$TickLength + new$y
      
      #Determine the secondary tick start and finish positions.
      new$x.sec    <- new$xend.grid
      new$y.sec    <- new$yend.grid
      new$xend.sec <- cos((new$Angle+180)*pi/180)*new$TickLength + new$x.sec
      new$yend.sec <- sin((new$Angle+180)*pi/180)*new$TickLength + new$y.sec
      
      ##ADD TO EXISTING
      return(rbind(existing,new))
    },error=function(e){
      warning(e)
    })
    return(existing)
  }
  
  angles      <- .get.angles(clockwise) + (!outside)*180
  angles.text <- .get.angles.ticklabels(clockwise) + .theme.get.rotation(self)
  
  ##get the base data.
  d = merge(data.frame(ix=1:length(seq.tlr)),data.frame(major=c(T,F)),by=NULL)
  d = ddply(d,c('major','ix'),function(x){
    ix = x$ix[1]
    .get.data(X=seq.tlr[ix],ix=ix, major = x$major[1], angle = angles[ix], angle.text = angles.text[ix]);
  })
  if(empty(d)){ return(items) }
  
  #FUNCTION TO RENDER TICKS AND LABELS
  .render.ticks  <- function(name,items,d,primary=TRUE){
    grob = zeroGrob()
    tryCatch({  
      e <- calc_element(name,theme=theme,verbose=F)
      if(identical(e,element_blank()))return(items)
      grob     <- segmentsGrob(
        x0 = ifthenelse(!primary,d$x.sec,d$x), 
        x1 = ifthenelse(!primary,d$xend.sec,d$xend),
        y0 = ifthenelse(!primary,d$y.sec,d$y), 
        y1 = ifthenelse(!primary,d$yend.sec,d$yend),
        default.units="npc",
        gp = gpar(col     = e$colour, 
                  lty     = e$linetype,
                  lineend = e$lineend,
                  lwd     = e$size*find_global_tern(".pt"))
      )
    },error = function(e){
      warning(e)
    })
    items[[length(items) + 1]] <- grob
    return(items)
  }
  
  #Axis labels, ie labels next to ticks
  .render.labels <- function(name,items,d){ 
    grob = zeroGrob()
    tryCatch({
      d         = d[which(d$Labels != ''),]
      e         = calc_element(name,theme=theme,verbose=F)
      if(empty(d) || identical(e,element_blank())) return(items)
      xts       = ifthenelse(outside,d$x,d$xend); xtf = ifthenelse(outside,d$xend,d$x) 
      yts       = ifthenelse(outside,d$y,d$yend); ytf = ifthenelse(outside,d$yend,d$y) 
      a         = is.numericor(e$angle,0) + is.numericor(unique(d$Angle.Text)[1],0)
      dA        = a - atan2(ytf-yts,xtf-xts)*180/pi                  #DEGREES, Angle Difference between Ticks and Labels
      hj        = +cos((dA-180)*pi/180)*0.5 + is.numericor(e$hjust,0) #BACK TO RADIANS
      vj        = -sin((dA-180)*pi/180)*0.5 + is.numericor(e$vjust,0) #BACK TO RADIANS
      grob      = textGrob( label          = label_formatter(as.character(d$Labels)), 
                             x             = ifthenelse(showprimary || !outside,xtf,xts) + convertX(cos(pi*(d$Angle + (!outside)*180)/180)*unit(2,'pt'),'npc',valueOnly = T),
                             y             = ifthenelse(showprimary || !outside,ytf,yts) + convertY(sin(pi*(d$Angle + (!outside)*180)/180)*unit(2,'pt'),'npc',valueOnly = T),
                             default.units ="npc", 
                             hjust         = is.numericor(hj,0.5),
                             vjust         = is.numericor(vj,0.5),
                             rot           = a, 
                             gp            = gpar(col        = e$colour, 
                                                  fontsize   = e$size,
                                                  fontfamily = ifthenelse(is.character(e$family),e$family,"sans"), 
                                                  fontface   = e$face, 
                                                  lineheight = ifthenelse(is.numeric(e$lineheight),e$lineheight,1)))
    },error = function(e){  warning(e)  })
    items[[length(items) + 1]] <- grob
    return(items)
  }
  .render.grid   <- function(name,items,d,showgrid.major=TRUE,showgrid.minor=TRUE){
    if((unique(d$Major) & showgrid.major) | (!unique(d$Major) & showgrid.minor)){
      grob = zeroGrob()
      tryCatch({  
        e    = calc_element(name,theme=theme,verbose=F)
        if(identical(e,element_blank()))return(items)
        grob = segmentsGrob(
          x0 = d$x, 
          x1 = d$xend.grid,
          y0 = d$y, 
          y1 = d$yend.grid,
          default.units="npc",
          gp = gpar(col     = e$colour, 
                    lty     = e$linetype,
                    lineend = e$lineend,
                    lwd     = e$size*find_global_tern(".pt"))
        )
      },error = function(e){ 
        warning(e)
      })
      items[[length(items) + 1]] <- grob
    }
    return(items)
  }
  
  #PROCESS TICKS AND LABELS
  if(onBackground){
    if(showgrid.major | showgrid.minor){
      for(n in unique(d$NameGrid)){ 
        items <- .render.grid(  name=n,items=items,d=d[which(d$NameGrid  == n),], showgrid.major=showgrid.major,showgrid.minor=showgrid.minor)
      } 
    }
  }else{
    if(showprimary)
      for(n in unique(d$NameTicks)){items <- .render.ticks(name=n,items=items,d=d[which(d$NameTicks == n),],primary=TRUE)}
    if(showsecondary)
      for(n in unique(d$NameTicks)){items <- .render.ticks(name=n,items=items,d=d[which(d$NameTicks == n),],primary=FALSE)}
    if(showlabels)
      for(n in unique(d$NameText)){ items <- .render.labels(name=n,items=items,d=d[which(d$NameText  == n),])}
  }
  items
}

.render.titles     <- function(self,data.extreme,scale_details,theme,items){
  sidelength = sqrt( diff(data.extreme$x[1:2])^2 + diff(data.extreme$y[1:2])^2)
  if(!.theme.get.showtitles(theme)) return(items)
  .render.title = function(name,ix,items){
    grob = zeroGrob()
    tryCatch({
      e     <- calc_element(name,theme=theme,verbose=F)
      if(identical(e,element_blank()))return(items)
      
      ixc   <- c('x','y')
      point <- as.numeric(data.extreme[ix,ixc])
      base  <- as.numeric(apply(data.extreme[-ix,c('x','y')],2,mean))
      angle <- atan2((point[2]-base[2])*.ratio(),point[1]-base[1])
      n     <- regmatches(name,regexpr(".$",name)) 
      l     <- c(self$scales[[n]]$name,self$labels_coord[[n]],self$labels_coord[[ self$mapping[[n]] ]],n)
      x     <- data.extreme$x[ix] + 0.05*sidelength*cos(angle)
      y     <- data.extreme$y[ix] + 0.05*sidelength*sin(angle)
      grob  <- textGrob(label = label_formatter(l[1]), 
                       x = x, y = y,
                       hjust  = e$hjust - 0.5*cos(angle), 
                       vjust  = e$vjust - 0.5*sin(angle),
                       rot    = e$angle,
                       vp     = viewport(clip='inherit'),   #Change to off???
                       gp     = gpar(col        = e$colour, 
                                   fontsize   = e$size,
                                   fontfamily = ifthenelse(is.character(e$family),e$family,"sans"), 
                                   fontface   = e$face, 
                                   lineheight = e$lineheight))
    items[[length(items) + 1]] <- grob
    },error = function(e){
      warning(e)
    })
    items
  }
  
  #process the axes
  items <- .render.title("tern.axis.title.T",1,items)
  items <- .render.title("tern.axis.title.L",2,items)
  items <- .render.title("tern.axis.title.R",3,items)
}

.render.arrows     <- function(self,data.extreme,details,theme,items){
  if(!.theme.get.showarrows(theme))return(items)
  tryCatch({
    
    #clockwise or anticlockwise precession
    clockwise <- .theme.get.clockwise(theme)
    
    #The basic data.
    d.s <- data.extreme[ifthenelse(clockwise,c(2,3,1),c(3,1,2)),c('x','y')]
    d.f <- data.extreme[c(1,2,3),c('x','y')]
    rownames(d.s) <- rownames(d.f) #Correct rownames
    
    #Determine the length of the side
    sidelength = sqrt( diff(data.extreme$x[1:2])^2 + diff(data.extreme$y[1:2])^2)
    
    #arrow start and finish proportions
    arrowstart = calc_element('tern.axis.arrow.start', theme)
    arrowfinish= calc_element('tern.axis.arrow.finish',theme)
    
    #Ensure arrow start and finish length is 3.
    if(length(arrowstart) != 3 && length(arrowstart) >= 1)
      arrowstart <- rep(arrowstart[1],3)
    if(length(arrowfinish) != 3 && length(arrowfinish) >= 1)
      arrowfinish <- rep(arrowfinish[1],3)
    
    #Itterate over indexes 1:3
    for(i in c(1:3)){
      #Put in correct order.
      if(arrowfinish[i] < arrowstart[i]){
        warning(paste("Arrow size theme 'element tern.axis.arrow.finish[",i,"]' (",arrowfinish[i],") is < 'tern.axis.arrow.start[",i,"]' (",arrowstart[i],"), values will be swapped.",sep=""),call.=FALSE)
        #swapvalues
        tmp  = arrowstart[i]; arrowstart[i]  = arrowfinish[i]; arrowfinish[i] = tmp
      }
      #Check finish
      if(arrowfinish[i] > 1.0){
        warning(paste("Arrow size theme 'element tern.axis.arrow.finish[",i,"]' (",arrowfinish[i],") is > 1.0 and will be truncated",sep=""),call.=FALSE)
        arrowfinish[i] = 1.0
      }
      #Check start
      if(arrowstart[i] < 0.0){
        warning(paste("Arrow size theme 'element tern.axis.arrow.start[",i,"]' (",arrowstart[i],") is < 0.0 and will be truncated",sep=""),call.=FALSE)
        arrowstart[i] = 0.0
      }
    }
    
    #Cut down to relative proportion.
    dx = (d.f - d.s)
    d.f <- d.f - (1-arrowfinish)*dx
    d.s <- d.s +      arrowstart*dx
    d <- rbind(d.s,d.f)
    
    #Determine the start and end positions
    ixseq <- names(self$mapping)
    ixrow <- paste0("AT.",ixseq)
    ixcol <- c("x","y","xend","yend")
    ix    <- which(colnames(d) %in% ixcol[c(1:2)])
    d     <- cbind(d[1:3,ix],d[4:6,ix]);
    rownames(d) <- ixrow; colnames(d) <- ixcol
    
    #The arrow seperation in npc units.
    arrowsep      <- calc_element("tern.axis.arrow.sep",theme=theme,verbose=F)
    ticklength    <- max(calc_element("tern.axis.ticks.length.major",theme=theme,verbose=F),
                         calc_element("tern.axis.ticks.length.minor",theme=theme,verbose=F))
    
    #Ensure there are EXACTLY 3 values for each metric
    if(length(arrowsep)   != 3 && length(arrowsep)   >= 1){ arrowsep   = rep(arrowsep[1],3)   }
    if(length(ticklength) != 3 && length(ticklength) >= 1){ ticklength = rep(ticklength[1],3) }
    
    #Determine the Angles
    d[ixrow,"angle"]    <- .get.angles.arrows(clockwise)
    
    #get set of 3 arrowsep positions
    d[ixrow,"arrowsep"] <- arrowsep*sidelength
    
    #MOVE the Arrows Off the Axes.
    rotation = .theme.get.rotation(self)
    d[,ixcol[c(1,3)]]   <- d[,ixcol[c(1,3)]] + cos(pi*(d$angle + rotation)/180)*d$arrowsep #xcoordinates
    d[,ixcol[c(2,4)]]   <- d[,ixcol[c(2,4)]] + sin(pi*(d$angle + rotation)/180)*d$arrowsep #ycoorinates
    
    #Centerpoints, labels, arrowsuffix
    d$xmn = rowMeans(d[,ixcol[c(1,3)]])
    d$ymn = rowMeans(d[,ixcol[c(2,4)]])
    d$L   = unlist(lapply(ixseq,function(n){ .theme.get.label(self,n)  }))
    d$LA  = unlist(lapply(ixseq,function(n){ 
      #c(self$labels_coord[[ sprintf('%sarrow',n) ]],)[1]
      .theme.get.label(self,n,suffix='arrow')
    }))
    d$W   = unlist(lapply('W',function(n){  .theme.get.label(self,n,'')  }))
    d$A   = .get.angles.arrowmarker(clockwise)
    d$AL  = .valid.angle(d$A + .theme.get.rotation(self))
    
    ##Function to create new arrow grob
    .render.arrow <- function(name,ix,items){
      grob = zeroGrob()
      tryCatch({  
        e    = calc_element(name,theme=theme,verbose=F)
        if(identical(e,element_blank()))return(items)
        grob = segmentsGrob(x0 = d$x[ix], x1 = d$xend[ix], y0 = d$y[ix], y1 = d$yend[ix],
                            default.units ="npc",
                            arrow         = e$lineend,
                            gp            = gpar(col     = e$colour, 
                                                 lty     = e$linetype,
                                                 lineend = 'butt',
                                                 lwd     = e$size*find_global_tern(".pt"))
        )
      },error = function(e){
        warning(e)
      })
      items[[length(items) + 1]] <- grob
      return(items)
    }
    
    #Function to greate new label grob
    .render.label <- function(name,ix,items){
      grob = zeroGrob()
      tryCatch({  
        e    = calc_element(name,theme=theme,verbose=F)
        if(identical(e,element_blank()))return(items)
        dA   = e$angle + d$AL - 180*(atan2((d$yend - d$y)*.ratio(), d$xend - d$x)/pi + as.numeric(clockwise))
        grob = textGrob( label = arrow_label_formatter(d$LA[ix],d$W[ix]), 
                         x     = d$xmn[ix] + convertX(cos(pi*(d$angle[ix] + rotation)/180)*unit(2,'pt'),'npc',valueOnly = T), 
                         y     = d$ymn[ix] + convertY(sin(pi*(d$angle[ix] + rotation)/180)*unit(2,'pt'),'npc',valueOnly = T), 
                         hjust = e$hjust + 0.5*sin(dA[ix]*pi/180), 
                         vjust = e$vjust + 0.5*cos(dA[ix]*pi/180),
                         rot   = d$AL[ix] + e$angle, 
                         default.units="npc", 
                         gp   = gpar(col        = e$colour, 
                                     fontsize   = e$size,
                                     fontfamily = e$family, 
                                     fontface   = e$face, 
                                     lineheight = e$lineheight))
      },error = function(e){ warning(e) })
      items[[length(items) + 1]] <- grob
      items
    }
  
    #process the axes
    for(i in 1:length(ixseq)){
      items <- .render.arrow(paste0("tern.axis.arrow.",     ixseq[i]),i,items) #Arrows
      items <- .render.label(paste0("tern.axis.arrow.text.",ixseq[i]),i,items) #Markers
    }
  },error=function(e){ 
    message(e)
  })
  items
}

.rotation = function (xy, angle, degrees=TRUE) {
  if(degrees) angle = pi*angle/180
  xy <- as.matrix(xy)
  ca <- cos(angle)
  sa <- sin(angle)
  xy.rot <- xy %*% t(matrix(c( ca, sa, 
                              -sa, ca), 2, 2))
  return(xy.rot)
}

.valid.angle = function(x){
  if(length(x) > 1)
    return(sapply(x,.valid.angle))
  if(x >  90) x = .valid.angle(x - 180)
  if(x <=-90) x = .valid.angle(x + 180)
  x
}

