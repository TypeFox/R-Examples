#' Spnet basic examples
#' 
#' Create SpatialNetwork object examples for demonstration and testing purpose.
#' 
#' @param map logical; if \code{TRUE} an example of map is provided.
#' @param color logical; if \code{TRUE} an example of map colorization is provided.
#' @param symbol logical; if \code{TRUE} an example of symbol use is provided.
#' @param network1 logical; if \code{TRUE} a first example of network is provided.
#' @param network2 logical; if \code{TRUE} a second example of network is provided.
#' @param barplot logical; if \code{TRUE} a example of barplot rendering of a numeric variable is provided.
#' @param title logical; if \code{TRUE} a example of title is provided.
#' @return a \code{SpatialNetwok} object.
#' @examples
#' data(world.map.simplified, package = "spnet")
#' net1 <- spnet.example.basic()
#' plot(net1)
#' @rdname spnet.example.basic
#' @importFrom utils data
#' @export
spnet.example.basic <- function(
  map = TRUE,
  color = TRUE,
  symbol = TRUE,
  network1 = TRUE,
  network2 = TRUE,
  barplot = TRUE,
  title = TRUE
) {
  example.basic.env <- new.env()
  data("world.map.simplified", package = "spnet", envir = example.basic.env)
  node <- c("France", "United States", "Brazil", "Australia", "Antarctica" )
  position <- match(node, example.basic.env$world.map.simplified@data[,'NAME']) - 1
  net1 <- spnet.create(
    data.frame(
      'NODE' =  node,
      'POSITION' = position
    )
  )
  graph.title.main(net1) <- ""
  if(map) {
    graph.map(net1) <- example.basic.env$world.map.simplified
  }
  if(color) {
    net1$continent <- c("Europa", "America", "America", "Oceania", "Antarctica")
    graph.color.variable(net1) <- "continent"
    graph.color.legend(net1) <- c('Europa' = "#CBB3E466", 'America' = "#D490B366", 'Oceania' = "#CBE4B366")
    graph.color.background(net1) <- "#B3E4E466" # light blue
    graph.color.border(net1) <- "#55555566" # grey
    graph.color.region(net1) <- "#D2A65F66" # light orange
  }
  if(symbol) {
    net1$role <- c('North', 'North', 'South', 'South', 'South')
    graph.symbol.variable(net1) <- 'role'
    graph.symbol.legend(net1) <- c('North' = 'triangle.up', 'South' = 'triangle.down')
    graph.symbol.color(net1) <- '#A52A2A88'
    graph.symbol.cex(net1) <- 1
    graph.symbol.shift.y(net1) <- 6
  }
  if(network1) {
    network1 <- matrix(
      rep(0, length(node)^2),
      nrow = length(node),
      dimnames = list(node, node)
    )
      
    graph.networks.add(net1) <- "network1"
    graph.network.data(net1, "network1") <- network1
    graph.network.data(net1, 'network1')['France', 'United States'] <- 2
    graph.network.data(net1, 'network1')['Australia', 'United States'] <- 1
    graph.network.data(net1, 'network1')['France', 'Brazil'] <- 3
    graph.network.data(net1, 'network1')['Brazil', 'France'] <- 2
    graph.network.label(net1, 'network1') <- 'Holidays'
    graph.network.arrow.shift.y(net1, 'network1') <- 2
    graph.network.arrow.color(net1, 'network1') <- '#33333366'
    graph.network.arrow.thickness(net1, 'network1') <- 0.5
  }
  if(network2) {
    network2 <- matrix(
      rep(0, length(node)^2),
      nrow = length(node),
      dimnames = list(node, node)
    )
      
    graph.networks.add(net1) <- "network2"
    graph.network.data(net1, "network2") <- network2
    graph.network.data(net1, 'network2')['Brazil', 'Australia'] <- 2
    graph.network.label(net1, 'network2') <- 'Studies'
    graph.network.arrow.shift.y(net1, 'network2') <- -2
    graph.network.arrow.opacity(net1, 'network2') <- 0.9
    graph.network.arrow.color(net1, 'network2') <- 'grey'
    graph.network.arrow.thickness(net1, 'network2') <- 0.5
  }
  if(barplot) {
    net1$num.var <- c(0.1,0.3,0.5,0.9,0.0)
    graph.barplot.variable(net1) <- "num.var"
    graph.barplot.bound.upper(net1) <- c(-13,20)
    graph.barplot.bound.lower(net1) <- c(-13,3)
    graph.barplot.fgcolor(net1) <- "#333333DD"
    graph.barplot.bgcolor(net1) <- "#E6E6E6DD"
    graph.barplot.width(net1) <- 10
  }
  if(title) {
    graph.title.main(net1) <- "Places visited by John, Elsa, Brian and Kate"
    graph.title.sub(net1) <- "For holidays and studies"
  }
  return(net1)
}

#' @rdname spnet.example.basic
#' @export
spnet.example.basic.full <- function(){
  net1 <- spnet.example.basic()
  return(net1)
}
#' @rdname spnet.example.basic
#' @export
spnet.example.basic.map <- function(){
  net1 <- spnet.example.basic(
    map = TRUE,
    color = FALSE,
    symbol = FALSE,
    network1 = FALSE,
    network2 = FALSE,
    barplot = FALSE
  )
  return(net1)
}