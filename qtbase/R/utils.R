## Holds some convenience wrappers for various objects

##' Get the dimensions of rectangles and rectangular objects
##'
##' @param x A rectangular object, like \code{QRect} or \code{QGraphicsItem}
##' @method dim QRectF
##' @rdname dim-methods
dim.QRectF <- function(x) c(x$width(), x$height())
##' @method dim QRect
##' @rdname dim-methods
dim.QRect <- dim.QRectF
##' @method dim QGraphicsScene
##' @rdname dim-methods
dim.QGraphicsScene <- function(x) dim(x$sceneRect)
##' @method dim QGraphicsItem
##' @rdname dim-methods
dim.QGraphicsItem <- function(x) dim(x$boundingRect)
##' @method dim QGraphicsView
##' @rdname dim-methods
dim.QGraphicsView <- function(x) dim(x$viewport()$rect)

##' Finds a child \code{QObject} by name. Mirrors a function in Qt
##' that is inaccessible via the ordinary interface because it is a
##' template.
##'
##' This is particularly useful when working with QtDesigner, where
##' objects in the UI file are named.
##' @title Find child by name
##' @param x The parent \code{QObject}
##' @param name The name of the child
##' @return The child \code{QObject}
##' @author Michael Lawrence
##' @examples
##' parent <- Qt$QObject()
##' child <- Qt$QObject(parent)
##' child$objectName <- "foo"
##' qfindChild(parent, "foo")

qfindChild <- function(x, name) {
  ## not exactly sure how to use this one
  Qt$QGlobalSpace$qt_qFindChild_helper(x, name, x$metaObject())
}

mapRTypeToQtType <- function(x) {
  map <- c("raw" = "unsigned char *",
           "logical" = "bool",
           "integer" = "int",
           "numeric" = "double",
           "character" = "QString")
  map[x]
}
