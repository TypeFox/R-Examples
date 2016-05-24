venn.diagram.crispr<-
function (x, filename=NULL, height = 3000, width = 3000, resolution = 500, 
          imagetype = "tiff", units = "px", compression = "lzw", na = "stop", 
          main = NULL, sub = NULL, main.pos = c(0.5, 1.05), main.fontface = "plain", 
          main.fontfamily = "serif", main.col = "black", main.cex = 1, 
          main.just = c(0.5, 1), sub.pos = c(0.5, 1.05), sub.fontface = "plain", 
          sub.fontfamily = "serif", sub.col = "black", sub.cex = 1, 
          sub.just = c(0.5, 1), category.names = names(x), force.unique = TRUE,
          ...) 
{
  VennDiagram::venn.diagram(...)
#   #requireNamespace("VennDiagram")
#   
#   if (force.unique) {
#     for (i in 1:length(x)) {
#       x[[i]] <- unique(x[[i]])
#     }
#   }
#   if ("none" == na) {
#     x <- x
#   }
#   else if ("stop" == na) {
#     for (i in 1:length(x)) {
#       if (any(is.na(x[[i]]))) {
#         stop("NAs in dataset", call. = FALSE)
#       }
#     }
#   }
#   else if ("remove" == na) {
#     for (i in 1:length(x)) {
#       x[[i]] <- x[[i]][!is.na(x[[i]])]
#     }
#   }
#   else {
#     stop("Invalid na option: valid options are \"none\", \"stop\", and \"remove\"")
#   }
#   if (0 == length(x) | length(x) > 5) {
#     stop("Incorrect number of elements.", call. = FALSE)
#   }
#   if (1 == length(x)) {
#     list.names <- category.names
#     if (is.null(list.names)) {
#       list.names <- ""
#     }
#     grob.list <- VennDiagram::draw.single.venn(area = length(x[[1]]), 
#                                                category = list.names, ind = FALSE, ...)
#   }
#   else if (2 == length(x)) {
#     grob.list <- VennDiagram::draw.pairwise.venn(area1 = length(x[[1]]), 
#                                                  area2 = length(x[[2]]), cross.area = length(intersect(x[[1]], 
#                                                                                                        x[[2]])), category = category.names, ind = FALSE, 
#                                                  ...)
#   }
#   else if (3 == length(x)) {
#     A <- x[[1]]
#     B <- x[[2]]
#     C <- x[[3]]
#     list.names <- category.names
#     nab <- intersect(A, B)
#     nbc <- intersect(B, C)
#     nac <- intersect(A, C)
#     nabc <- intersect(nab, C)
#     grob.list <- VennDiagram::draw.triple.venn(area1 = length(A), 
#                                                area2 = length(B), area3 = length(C), n12 = length(nab), 
#                                                n23 = length(nbc), n13 = length(nac), n123 = length(nabc), 
#                                                category = list.names, ind = FALSE, list.order = 1:3, 
#                                                ...)
#   }
#   else if (4 == length(x)) {
#     A <- x[[1]]
#     B <- x[[2]]
#     C <- x[[3]]
#     D <- x[[4]]
#     list.names <- category.names
#     n12 <- intersect(A, B)
#     n13 <- intersect(A, C)
#     n14 <- intersect(A, D)
#     n23 <- intersect(B, C)
#     n24 <- intersect(B, D)
#     n34 <- intersect(C, D)
#     n123 <- intersect(n12, C)
#     n124 <- intersect(n12, D)
#     n134 <- intersect(n13, D)
#     n234 <- intersect(n23, D)
#     n1234 <- intersect(n123, D)
#     grob.list <- VennDiagram::draw.quad.venn(area1 = length(A), 
#                                              area2 = length(B), area3 = length(C), area4 = length(D), 
#                                              n12 = length(n12), n13 = length(n13), n14 = length(n14), 
#                                              n23 = length(n23), n24 = length(n24), n34 = length(n34), 
#                                              n123 = length(n123), n124 = length(n124), n134 = length(n134), 
#                                              n234 = length(n234), n1234 = length(n1234), category = list.names, 
#                                              ind = FALSE, ...)
#   }
#   else if (5 == length(x)) {
#     A <- x[[1]]
#     B <- x[[2]]
#     C <- x[[3]]
#     D <- x[[4]]
#     E <- x[[5]]
#     list.names <- category.names
#     n12 <- intersect(A, B)
#     n13 <- intersect(A, C)
#     n14 <- intersect(A, D)
#     n15 <- intersect(A, E)
#     n23 <- intersect(B, C)
#     n24 <- intersect(B, D)
#     n25 <- intersect(B, E)
#     n34 <- intersect(C, D)
#     n35 <- intersect(C, E)
#     n45 <- intersect(D, E)
#     n123 <- intersect(n12, C)
#     n124 <- intersect(n12, D)
#     n125 <- intersect(n12, E)
#     n134 <- intersect(n13, D)
#     n135 <- intersect(n13, E)
#     n145 <- intersect(n14, E)
#     n234 <- intersect(n23, D)
#     n235 <- intersect(n23, E)
#     n245 <- intersect(n24, E)
#     n345 <- intersect(n34, E)
#     n1234 <- intersect(n123, D)
#     n1235 <- intersect(n123, E)
#     n1245 <- intersect(n124, E)
#     n1345 <- intersect(n134, E)
#     n2345 <- intersect(n234, E)
#     n12345 <- intersect(n1234, E)
#     grob.list <- VennDiagram::draw.quintuple.venn(area1 = length(A), 
#                                                   area2 = length(B), area3 = length(C), area4 = length(D), 
#                                                   area5 = length(E), n12 = length(n12), n13 = length(n13), 
#                                                   n14 = length(n14), n15 = length(n15), n23 = length(n23), 
#                                                   n24 = length(n24), n25 = length(n25), n34 = length(n34), 
#                                                   n35 = length(n35), n45 = length(n45), n123 = length(n123), 
#                                                   n124 = length(n124), n125 = length(n125), n134 = length(n134), 
#                                                   n135 = length(n135), n145 = length(n145), n234 = length(n234), 
#                                                   n235 = length(n235), n245 = length(n245), n345 = length(n345), 
#                                                   n1234 = length(n1234), n1235 = length(n1235), n1245 = length(n1245), 
#                                                   n1345 = length(n1345), n2345 = length(n2345), n12345 = length(n12345), 
#                                                   category = list.names, ind = FALSE, ...)
#   }
#   else {
#     stop("Invalid size of input object")
#   }
#   if (!is.null(sub)) {
#     grob.list <- VennDiagram::add.title(gList = grob.list, x = sub, pos = sub.pos, 
#                            fontface = sub.fontface, fontfamily = sub.fontfamily, 
#                            col = sub.col, cex = sub.cex)
#   }
#   if (!is.null(main)) {
#     grob.list <- VennDiagram::add.title(gList = grob.list, x = main, pos = main.pos, 
#                            fontface = main.fontface, fontfamily = main.fontfamily, 
#                            col = main.col, cex = main.cex)
#   }
#   if (!is.null(filename)) {
#     current.type <- getOption("bitmapType")
#     if (length(grep("Darwin", Sys.info()["sysname"]))) {
#       options(bitmapType = "quartz")
#     }
#     else {
#       options(bitmapType = "cairo")
#     }
#     if ("tiff" == imagetype) {
#       tiff(filename = filename, height = height, width = width, 
#            units = units, res = resolution, compression = compression)
#     }
#     else if ("png" == imagetype) {
#       png(filename = filename, height = height, width = width, 
#           units = units, res = resolution)
#     }
#     else if ("svg" == imagetype) {
#       svg(filename = filename, height = height, width = width)
#     }
#     else {
#       stop("You have misspelled your 'imagetype', please try again")
#     }
#     grid::grid.draw(grob.list)
#     dev.off()
#     options(bitmapType = current.type)
#     return(1)
#   }
#   plot.new()
#   grid::grid.draw(grob.list)
  #return(grob.list)
}
#<environment: namespace:VennDiagram>