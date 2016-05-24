#' Radial Coordinates Plots
#'
#' Creates radial coordinates plots, with m variables (rows) and n samples (columns).
#' Each variable is mapped onto a circle, using data points as spring constants.
#' Each column is re-scaled to have numeric values between 0 and 1.
#'
#' @param dat A matrix with \code{m} rows and \code{n} columns, where columns represent dimensions.
#' @param group A vector of length \code{m}, specifying groups (e.g., phenotypes or conditions for \code{m} samples).
#' @param color A vector of hex color codes to represent groups.
#' @param alpha A numeric value for transparency.
#' @param vjust A parameter to vertically adjust axis names around the circle; \code{vjust} arguments for \code{geom_text}
#' @param hjust A parameter to horizontally adjust axis names around the circle; \code{hjust} arguments for \code{geom_text}
#'
#' @return \code{svd.radial} creates and draws a figure, which is a \code{ggplot} object.
#'
#' @export radviz
#' @author Neo Christopher Chung \email{nchchung@gmail.com}
#' @seealso \link{svd.radial}
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @import ggplot2
#'
#' @references Ankerst M., Keim D. A., Kriegel H.-P. Circle Segments: A Technique for Visually Exploring Large Multidimensional Data Sets, IEEE Visualization, 1996.
#' @references K.A. Olsen, R.R. Korfhage, K.M. Sochats, M.B. Spring and J.G. Williams. Visualisation of a Document Collection: The VIBE System, Information Processing and Management, Vol. 29, No. 1, pp. 69-81, Pergamon Press Ltd, 1993.
#'
#' @examples
#' set.seed(1234)
#' dat = matrix(rnorm(9*4), 9, 4, dimnames=list(paste(1:9),letters[1:4]))
#' radviz(dat, group=c(rep("Group 1",3),rep("Group 2",3),rep("Group 3",3)))
radviz <- function(dat,group=NULL,color=NULL,hjust=0,vjust=0,alpha=1) {
  m=dim(dat)[1]
  n=dim(dat)[2]

  if(is.null(group)) group = rep(1, nrow(dat))
  if(is.null(colnames(dat))) colnames(dat) = paste0("V",1:n)
  group=as.factor(group)
  nGroup = length(unique(group))
  nameGroup = levels(group)

  countGroup=1:nGroup
  cnames=colnames(dat)

  # set color
  if(is.null(color)) {
    if(length(unique(group)) > 12) {
      warning("Too many groups can be difficult to distinguish in a plot.")
      color = colorRampPalette(brewer.pal(12, "Set3"))(nGroup)
    } else if(nGroup > 9 & nGroup <= 12) {
      color = brewer.pal(n=nGroup, name="Set3")
    } else if(nGroup > 2 & nGroup <= 9) {
      color = brewer.pal(n=nGroup, name="Set1")
    } else if(nGroup == 2) {
      color = c("#377EB8","#FF7F00")
    } else if(nGroup == 1) {
      color = c("#000000")
    }
  }

    #created projection of each observation
    dat.norm = apply(dat, 2, scales::rescale)
    # if minval, maxval are different from 0,1
    # dat.norm = apply(dat, 2, function(x) scales::rescale(x, to=c(min,max)))

    # projection codes adopted from radviz2 function from dprep
    sumrows=rowSums(dat.norm)
    columns=seq(0,(n-1))
    angles=(2*pi*columns)/n
    cosines=cos(angles)
    sines=sin(angles)
    proj.x=(dat.norm %*% cosines)
    proj.x=proj.x/sumrows
    proj.y=(dat.norm %*% sines)
    proj.y=proj.y/sumrows
    data.rad = data.frame(x=proj.x, y=proj.y, group=group)

  dat <- circleFun(c(0,0),diameter = 2,npoints = 100)
  #geom_path will do open circles, geom_polygon will do filled circles
  g = ggplot(dat, aes_string(x="x",y="y")) + geom_path() + theme_bw() + ylim(-1.2,1.2) + xlim(-1.2,1.2) +
    theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
  ################
  dimensions = data.frame(cosines, sines, label=cnames, hjust=hjust, vjust=vjust)
  g = g + addticks(dimensions)
  g = g + addpoints(data.rad, group, alpha=alpha)
  return(g)
}

# Making circle using joran's function http://stackoverflow.com/questions/6862742/draw-a-circle-with-ggplot2
circleFun <- function(center = c(0,0), diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

# Add ticks around a circle for a radviz
addticks <- function(dataticks) {
  list(
    geom_point(data=dataticks, aes_string(x="cosines", y="sines")),
    geom_text(data=dataticks, aes_string(x="cosines", y="sines", label="label", hjust="hjust", vjust="vjust"))
  )
}

# Add data points to a radviz
addpoints <- function(data.rad, group=1, alpha=1) {
  dat = cbind(data.rad, group)
  if(length(unique(group)) == 1) {
    list(
      geom_point(data=dat, aes_string(x="x", y="y"), alpha=alpha)
    )
  } else {
    list(
      geom_point(data=dat, aes_string(x="x", y="y", colour="group"), alpha=alpha)
    )
  }
}
