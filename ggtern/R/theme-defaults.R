#' Complete Themes
#' 
#' \code{ggtern} ships with a number of complete themes:
#' \itemize{
#'  \item{ Black and White Theme:}{ 
#'    \code{\link[=theme_tern_bw]{theme_bw(...)}}
#'  }
#'  \item{Minimal Theme:}{
#'    \code{\link[=theme_tern_minimal]{theme_minimal(...)}}
#'  }
#'  \item{Classic Theme:}{
#'    \code{\link[=theme_tern_classic]{theme_classic(...)}}
#'  }
#'  \item{Gray and White Theme:}{
#'    \code{\link[=theme_tern_gray]{theme_gray(...)}}
#'  }
#'  \item{Red, Green, Blue and White Theme:}{
#'    \code{\link[=theme_tern_rgbw]{theme_rgbw(...)}}
#'  }
#'  \item{Red, Green, Blue and Gray Theme:}{
#'    \code{\link[=theme_tern_rgbg]{theme_rgbg(...)}}
#'  }
#'  \item{Dark Theme:}{
#'    \code{\link[=theme_dark]{theme_dark(...)}}
#'  }
#'  \item{Darker Theme:}{
#'    \code{\link[=theme_darker]{theme_darker(...)}}
#'  }
#'  \item{Light Theme:}{
#'    \code{\link[=theme_light]{theme_light(...)}}
#'  }
#'  \item{Theme with Only Black Lines:}{
#'    \code{\link[=theme_linedraw]{theme_linedraw(...)}}
#'  }
#' }
#' @rdname theme_complete
#' @name theme_complete
#' @author Nicholas Hamilton
NULL

#' ggtern themes
#'
#' Themes set the general aspect of the plot such as the colour of the
#' background, gridlines, the size and colour of fonts.
#'
#' @param base_size base font size
#' @param base_family base font family
#'
#' @details \describe{
#'
#' \item{\code{theme_gray}}{
#' The signature ggplot2 theme with a grey background and white gridlines,
#' designed to put the data forward yet make comparisons easy.}
#'
#' \item{\code{theme_bw}}{
#' The classic dark-on-light ggplot2 theme. May work better for presentations
#' displayed with a projector.}
#'
#' \item{\code{theme_linedraw}}{
#' A theme with only black lines of various widths on white backgrounds,
#' reminiscent of a line drawings. Serves a purpose similar to \code{theme_bw}.
#' Note that this theme has some very thin lines (<< 1 pt) which some journals
#' may refuse.}
#'
#' \item{\code{theme_light}}{
#' A theme similar to \code{theme_linedraw} but with light grey lines and axes,
#' to direct more attention towards the data.}
#'
#' \item{\code{theme_dark}}{
#' The dark cousin of \code{theme_light}, with similar line sizes but a dark background. 
#' Useful to make thin coloured lines pop out.
#' }
#'
#' \item{\code{theme_darker}}{
#' A darker cousing to \code{theme_dark}, with a dark panel background.
#' }
#'
#' \item{\code{theme_minimal}}{
#' A minimalistic theme with no background annotations.
#' }
#'
#' \item{\code{theme_classic}}{
#' A classic-looking theme, with x and y axis lines and no gridlines.
#' }
#' 
#' \item{\code{theme_rgbw}}{
#'  A theme with white background, red, green and blue axes and gridlines
#' }
#' 
#' \item{\code{theme_rgbg}}{
#' A theme with grey background, red, green and blue axes and gridlines
#' }
#'
#' \item{\code{theme_void}}{ 
#' A completely empty theme.
#' }
#' 
#' \item{\code{theme_custom}}{
#' Theme with custom basic colours
#' }
#' }
#'
#' @examples
#' data(Feldspar)
#' p <- ggtern(Feldspar,aes(Ab,An,Or)) + 
#'      geom_point(aes(colour=T.C,size=P.Gpa)) + 
#'      facet_wrap(~Feldspar)
#' 
#' #Uncomment to run
#' p + theme_gray()
#' p + theme_rgbg()
#' p + theme_dark()
#' @aliases theme_tern_gray theme_tern_grey theme_grey theme_tern_bw theme_tern_classic theme_tern_rgbg theme_tern_rgbw theme_tern_minimal
#' @name ggtern_themes
#' @rdname ggtern_themes
NULL

#' @rdname ggtern_themes
#' @export
theme_gray  <- function(base_size = 12, base_family = ""){
  .theme_tern(base_size=base_size, base_family=base_family,
              tern.panel.background = 'white',
              tern.plot.background  = 'grey92',
              col.T                 = "gray50",
              col.L                 = "gray50",
              col.R                 = "gray50",
              col.axis.T            = "grey92",
              col.axis.L            = "grey92",
              col.axis.R            = "grey92",
              col.title.T           = "black",
              col.title.L           = "black",
              col.title.R           = "black",
              tern.axis.size        = 0.50,
              col.grid.minor        = 'gray99',
              #col.grid.minor        = 'white',
              ticklength.minor      = unit(0,"npc"),
              showarrow             = FALSE,
              col.arrow.text.T="black",col.arrow.text.L="black",col.arrow.text.R="black")
}
theme_tern_gray <- function(base_size = 12, base_family = ""){
  tern_dep("1.0.1.3","theme_tern_gray has been superceded by the ggplot2 standard theme_gray")
  theme_gray(base_size,base_family)
}
theme_tern_grey <- theme_tern_gray
theme_grey <- theme_gray

#' @rdname ggtern_themes
#' @export
theme_bw <- function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text                  = element_text(size   = rel(0.8)),
      axis.ticks                 = element_line(colour = "black"),
      legend.key                 = element_rect(colour = "grey80"),
      tern.panel.background      = element_rect(fill   = 'white', colour = 'white'),
      tern.plot.background       = element_rect(fill   = 'white', colour = 'white'),
      tern.axis.line             = element_line(colour = 'grey50',size=1.0),
      tern.panel.grid.major      = element_line(colour = "grey90"),
      tern.panel.grid.minor      = element_line(colour = "grey98", size = 0.25),
      strip.background           = element_rect(fill   = "grey80", colour = "grey50", size = 0.2)
    )
}
theme_tern_bw <- function(base_size = 12, base_family = ""){
  tern_dep("1.0.1.3","theme_tern_bw has been superceded by the ggplot2 standard theme_bw")
  theme_bw(base_size,base_family)
}


#' @rdname ggtern_themes
#' @export
theme_linedraw <- function(base_size = 12, base_family = "") {
  half_line <- base_size / 2
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      tern.axis.line    = element_line(colour = 'black', size = 1.0),
      legend.key        = element_rect(colour = "black", size = 0.25),
      strip.background  = element_rect(fill   = "black", colour = NA),
      strip.text.x      = element_text(
        colour = "white",
        margin = margin(t = half_line, b = half_line)
      ),
      strip.text.y      = element_text(
        colour = "white",
        angle = 90,
        margin = margin(l = half_line, r = half_line)
      )
    )
}


#' @rdname ggtern_themes
#' @export
theme_classic <- function(base_size = 12, base_family = ""){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
  theme(
    panel.grid.major = element_line(),
    panel.grid.minor = element_line(),
    strip.background = element_rect(colour = "black", size = 0.5),
    legend.key       = element_blank()
  )
}

theme_tern_classic <- function(base_size = 12, base_family = ""){
  tern_dep("1.0.1.3","theme_tern_classic has been superceded by the ggplot2 standard theme_classic")
  theme_classic(base_size,base_family)
}

#' @rdname ggtern_themes
#' @export
theme_void <- function(base_size = 12, base_family = "") {
  theme_bw() %+replace%
  theme(
    line =               element_blank(),
    rect =               element_blank(),
    text =               element_blank(),
    complete = TRUE
  )
}

#' A theme with grey background, red, green and blue axes and white gridlines
#' 
#' \code{theme_rgbg} is a theme with grey background, red, green and blue axes and gridlines
#' @rdname ggtern_themes
#' @export
theme_rgbg  <- function(base_size = 12, base_family = ""){
  theme_rgbw() %+replace%
    theme(
      tern.panel.background    = element_rect(fill='white'),
      tern.plot.background     = element_rect(fill='gray92'),
      tern.panel.grid.minor    = element_line(colour='white',size=0.25)
    )
}
theme_tern_rgbg <- function(base_size = 12, base_family = ""){
  tern_dep("1.0.1.3","theme_tern_rgbg has been superceded by theme_rgbg")
  theme_rgbg(base_size,base_family)
}
theme_rgb <- theme_rgbg

#' @rdname ggtern_themes
#' @export
theme_rgbw  <- function(base_size = 12, base_family = ""){
  .theme_tern(base_size=base_size, base_family=base_family,
              tern.plot.background ="white",
              tern.panel.background = 'white',
              col.T          ="darkred",
              col.L          ="darkblue",
              col.R          ="darkgreen",
              col.grid.T     ="darkred",
              col.grid.L     ="darkblue",
              col.grid.R     ="darkgreen",
              col.grid.minor ="gray90",grid.linetype=6,grid.linetype.minor=1,grid.major.size=0.25)
}
theme_tern_rgbw <- function(base_size = 12, base_family = ""){
  tern_dep("1.0.1.3","theme_tern_rgbw has been superceded by theme_rgbw")
  theme_rgbw(base_size,base_family)
}

#' @rdname ggtern_themes
#' @export
theme_minimal <- function(base_size = 12, base_family = "") {
  # Starts with theme_bw and then modify some parts
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      legend.background      = element_blank(),
      legend.key             = element_blank(),
      panel.background       = element_rect(),
      tern.panel.background  = element_rect(fill='white',colour=NA),
      panel.border           = element_blank(),
      strip.background       = element_blank(),
      plot.background        = element_blank(),
      axis.ticks             = element_blank(),
      axis.ticks.length      = unit(0, "lines")
    )
}


#' @rdname ggtern_themes
#' @export
theme_dark <- function(base_size = 12, base_family = "") {
  base = ggplot2::theme_dark()
  base %+replace%
    theme(
      tern.panel.background    = element_rect(fill   = "white", colour = NA),
      tern.plot.background     = base$panel.background,
      tern.panel.grid.major    = base$panel.grid.major,
      tern.panel.grid.major.T  = element_line(),
      tern.panel.grid.major.L  = element_line(),
      tern.panel.grid.major.R  = element_line(),
      tern.panel.grid.minor    = base$panel.grid.minor,
      tern.axis                = base$panel.grid.major,
      tern.axis.line           = element_line(),
      tern.axis.line.T         = element_line(),
      tern.axis.line.L         = element_line(),
      tern.axis.line.R         = element_line()
    )
}

#' @rdname ggtern_themes
#' @export
theme_darker <- function(base_size = 12, base_family = "") {
  base = theme_dark(base_size=base_size,base_family=base_family) 
  base %+replace%
    theme(tern.panel.background = element_rect(fill='grey75',colour=NULL),
          plot.background       = element_rect(fill='grey75',colour=NULL),
          legend.background     = element_rect(fill='grey50',colour=NULL))
}

#' @rdname ggtern_themes
#' @export
theme_light <- function(base_size = 12, base_family = "") {
  base = ggplot2::theme_light(base_size=base_size,base_family=base_family) 
  base %+replace%
    theme(
      tern.panel.background    = element_rect(fill   = "white", colour = NA),
      tern.plot.background     = base$panel.background,
      tern.panel.grid.major    = base$panel.grid.major,
      tern.panel.grid.major.T  = element_line(),
      tern.panel.grid.major.L  = element_line(),
      tern.panel.grid.major.R  = element_line(),
      tern.panel.grid.minor    = base$panel.grid.minor
    )
}


#Internals
#helper function
.theme_tern      <- function(base_size             = 12, 
                             base_family           = "",
                             base_ggplot2_theme    = "theme_gray",
                             tern.panel.background = NA,
                             tern.plot.background  = NA,
                             col.T               = "black",
                             col.L               = "black",
                             col.R               = "black",
                             col.grid.T          = "white",
                             col.grid.L          = "white",
                             col.grid.R          = "white",
                             col.grid.minor      = "grey95",
                             col.axis.T          = col.T,
                             col.axis.L          = col.L,
                             col.axis.R          = col.R,
                             col.arrow.T         = col.T,
                             col.arrow.L         = col.L,
                             col.arrow.R         = col.R,
                             col.title.T         = col.T,
                             col.title.L         = col.L,
                             col.title.R         = col.R,
                             col.arrow.text.T    = col.T,
                             col.arrow.text.L    = col.L,
                             col.arrow.text.R    = col.R,
                             #col.ticks.major     = "black",
                             tern.axis.size      = 0.5,
                             showarrow           = getOption("tern.showarrows"),
                             ticks.outside       = getOption("tern.ticks.outside"),
                             ticks.showsecondary = getOption("tern.ticks.showsecondary"),
                             ticks.showprimary   = getOption("tern.ticks.showprimary"),
                             grid.linetype       = 1,
                             grid.linetype.minor = grid.linetype,
                             grid.major.size     = tern.axis.size,
                             grid.minor.size     = grid.major.size / 2,
                             ticklength.major    = unit(0.010,"npc"),
                             ticklength.minor    = unit(0.005,"npc")){
  
  #TEXT SIZES
  size.base      <- max(base_size-4,2)
  size.text      <- max(base_size-2,4)
  size.title     <- max(base_size-0,6)
  
  get(base_ggplot2_theme,asNamespace("ggplot2"))(
    base_size=base_size,base_family=base_family) %+replace%
    theme(
      tern.panel.background         = element_rect(fill=tern.panel.background,color=NA),
      tern.plot.background          = element_rect(fill=tern.plot.background,color=NA),
      tern.axis.clockwise           = getOption("tern.clockwise"),
      tern.axis.arrow.show          = showarrow,
      tern.axis.title.show          = getOption("tern.showtitles"),
      tern.axis.text.show           = getOption("tern.showlabels"),
      tern.axis.arrow.start         = getOption("tern.arrow.start"),
      tern.axis.arrow.finish        = getOption("tern.arrow.finish"),
      tern.axis.hshift              = getOption("tern.hshift"),
      tern.axis.vshift              = getOption("tern.vshift"),
      tern.axis.arrow.sep           = as.numeric(getOption("tern.arrowsep")),
      
      tern.axis                     = element_line(size=tern.axis.size,linetype="solid"),
      tern.axis.line                = element_line(colour=.resolveCol(col.axis.T,col.axis.L,col.axis.R,0)),
      tern.axis.line.T              = element_line(colour=.resolveCol(col.axis.T,col.axis.L,col.axis.R,1)),
      tern.axis.line.L              = element_line(colour=.resolveCol(col.axis.T,col.axis.L,col.axis.R,2)),
      tern.axis.line.R              = element_line(colour=.resolveCol(col.axis.T,col.axis.L,col.axis.R,3)),
      
      tern.axis.arrow               = element_line(colour=.resolveCol(col.arrow.T,col.arrow.L,col.arrow.R,0),lineend=arrow(length=unit(2.5,"mm"))),
      tern.axis.arrow.T             = element_line(colour=.resolveCol(col.arrow.T,col.arrow.L,col.arrow.R,1)),
      tern.axis.arrow.L             = element_line(colour=.resolveCol(col.arrow.T,col.arrow.L,col.arrow.R,2)),
      tern.axis.arrow.R             = element_line(colour=.resolveCol(col.arrow.T,col.arrow.L,col.arrow.R,3)),
      
      tern.axis.text                = element_text(colour=.resolveCol(col.T,col.L,col.R,0),size=size.base,face="plain"),
      tern.axis.text.T              = element_text(colour=.resolveCol(col.T,col.L,col.R,1)),
      tern.axis.text.L              = element_text(colour=.resolveCol(col.T,col.L,col.R,2)),
      tern.axis.text.R              = element_text(colour=.resolveCol(col.T,col.L,col.R,3)),
      
      tern.axis.arrow.text          = element_text(colour=.resolveCol(col.arrow.text.T,col.arrow.text.L,col.arrow.text.R,0),size=size.text,face="plain"),
      tern.axis.arrow.text.T        = element_text(colour=.resolveCol(col.arrow.text.T,col.arrow.text.L,col.arrow.text.R,1)),
      tern.axis.arrow.text.L        = element_text(colour=.resolveCol(col.arrow.text.T,col.arrow.text.L,col.arrow.text.R,2)),
      tern.axis.arrow.text.R        = element_text(colour=.resolveCol(col.arrow.text.T,col.arrow.text.L,col.arrow.text.R,3)),
      
      tern.axis.title               = element_text(colour=.resolveCol(col.title.T,col.title.L,col.title.R,0),size=size.title,face="bold"),
      tern.axis.title.T             = element_text(colour=.resolveCol(col.title.T,col.title.L,col.title.R,1)),
      tern.axis.title.L             = element_text(colour=.resolveCol(col.title.T,col.title.L,col.title.R,2)),
      tern.axis.title.R             = element_text(colour=.resolveCol(col.title.T,col.title.L,col.title.R,3)),
      
      tern.panel.grid               = element_line(linetype=grid.linetype),
      tern.panel.grid.major         = element_line(colour = .resolveCol(col.grid.T,col.grid.L,col.grid.R,0),size=grid.major.size),
      tern.panel.grid.major.T       = element_line(colour = .resolveCol(col.grid.T,col.grid.L,col.grid.R,1)),
      tern.panel.grid.major.L       = element_line(colour = .resolveCol(col.grid.T,col.grid.L,col.grid.R,2)),
      tern.panel.grid.major.R       = element_line(colour = .resolveCol(col.grid.T,col.grid.L,col.grid.R,3)),
      tern.panel.grid.minor         = element_line(colour = col.grid.minor,size=grid.minor.size,linetype=grid.linetype.minor),
      
      tern.axis.ticks.outside       = ticks.outside,
      tern.axis.ticks.secondary.show= ticks.showsecondary,
      tern.axis.ticks.primary.show  = ticks.showprimary,
      
      tern.axis.ticks.length.major  = ticklength.major,
      tern.axis.ticks.length.minor  = ticklength.minor,
      
      tern.axis.ticks               = element_line(),
      tern.axis.ticks.major         = element_line(colour = .resolveCol(col.T,col.L,col.R,0), size = grid.major.size ),
      tern.axis.ticks.major.T       = element_line(colour = .resolveCol(col.T,col.L,col.R,1)),
      tern.axis.ticks.major.L       = element_line(colour = .resolveCol(col.T,col.L,col.R,2)),
      tern.axis.ticks.major.R       = element_line(colour = .resolveCol(col.T,col.L,col.R,3)),
      
      tern.axis.ticks.minor         = element_line(colour = .resolveCol(col.T,col.L,col.R,0), size=grid.minor.size),
      tern.axis.ticks.minor.T       = element_line(colour = .resolveCol(col.T,col.L,col.R,1)),
      tern.axis.ticks.minor.L       = element_line(colour = .resolveCol(col.T,col.L,col.R,2)),
      tern.axis.ticks.minor.R       = element_line(colour = .resolveCol(col.T,col.L,col.R,3)),
      
      tern.panel.expand             = getOption('tern.expand'),
      tern.panel.rotate             = 0,
      tern.panel.grid.ontop         = FALSE,
      tern.axis.line.ontop          = FALSE
    )
}

#' @param col.T colour of top axis, ticks labels and major gridlines
#' @param col.L colour of left axis, ticks, labels and major gridlines
#' @param col.R colour of right axis, ticks, labels and major gridlines
#' @param col.BG colour of the plot background area
#' @param tern.plot.background colour of background colour to plot area
#' @param tern.panel.background colour of panel background of plot area
#' @param col.grid.minor the colour of the minor grid
#' \code{theme_custom} is a convenience function to allow the user to control the basic theme colours very easily.
#' @rdname ggtern_themes
#' @export
theme_custom  <- function(base_size = 12,
                          base_family = "",
                          tern.plot.background  = 'grey92',
                          tern.panel.background = 'white',
                          col.T='gray95',
                          col.L='gray95',
                          col.R='gray95',
                          col.BG="transparent",
                          col.grid.minor="gray90"){
  .theme_tern(base_size=base_size, base_family=base_family, 
              tern.plot.background  = tern.plot.background,
              tern.panel.background = tern.panel.background,
              col.T      = col.T,
              col.L      = col.L,
              col.R      = col.R,
              col.grid.T = col.T,
              col.grid.L = col.L,
              col.grid.R = col.R,
              col.grid.minor = col.grid.minor,
              grid.linetype=6,
              grid.linetype.minor=1,
              grid.major.size=0.25)
}

.resolveCol <- function(t,l,r,ix=0){
  stopifnot( (ix %in% 0:3) )
  if(t == l && l == r) ##ALL THE SAME
    if(ix == 0) return(t) else return(NULL)
  else if(ix == 0) return(NULL) else return( c(t,l,r)[ix])
}


