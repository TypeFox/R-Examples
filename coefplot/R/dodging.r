# #' collidev
# #' 
# #' Vertical collision checking
# #' 
# #' A hacky adaptation of ggplot's collide function to be used for vertical collisions.  No warranties on this working except that it looks good for now.
# #' 
# #' The adaptation switch x's to y's and y's to x's and width to height
# #' 
# #' @author Jared P. Lander
# #' @aliases collidev
# #' @import ggplot2
# #' @export collidev
# #' @return Not sure
# #' @param data data
# #' @param height How much to dodge the items
# #' @param name Refer to ggplot2
# #' @param strategy Refer to ggplot2
# #' @param check.height Refer to ggplot2
# #' @examples
# #' # none here
# collidev <- function(data, height = NULL, name, strategy, check.height = TRUE) 
# {
#     if (!is.null(height)) {
#         if (!(all(c("ymin", "ymax") %in% names(data)))) {
# #             data <- within(data, {
# #                 ymin <- y - height/2
# #                 ymax <- y + height/2
# #             })
#             data[, c("ymax", "ymin")] <- data$y + matrix(c(1, -1), ncol=2, nrow=NROW(data), byrow=TRUE)*height/2
#         }
#     }
#     else {
#         if (!(all(c("ymin", "ymax") %in% names(data)))) {
#             data$ymin <- data$y
#             data$ymax <- data$y
#         }
#         heights <- unique(with(data, ymax - ymin))
#         heights <- heights[!is.na(heights)]
#         if (!zero_range(range(heights))) {
#             warning(name, " requires constant height: output may be incorrect", 
#                     call. = FALSE)
#         }
#         height <- heights[1]
#     }
#     data <- data[order(data$ymin), ]
#     intervals <- as.numeric(t(unique(data[c("ymin", "ymax")])))
#     intervals <- intervals[!is.na(intervals)]
#     if (length(unique(intervals)) > 1 & any(diff(scale(intervals)) < 
#                                                 -1e-06)) {
#         warning(name, " requires non-overlapping y intervals", 
#                 call. = FALSE)
#     }
#     if (!is.null(data$xmax)) {
#         # this line is commented out and replaced with the one below to deal with CRANs abhorence of non-visible bindings
#         #ddply(data, .(ymin), strategy, height = height)
#         ddply(data, "ymin", strategy, height = height)
#     }
#     else if (!is.null(data$x)) {
#         message("xmax not defined: adjusting position using x instead")
#         # this line is commented out and replaced with the ones below to deal with CRANs abhorence of non-visible bindings
# #        transform(ddply(transform(data, xmax = x), .(ymin), strategy, 
# #                        height = height), x = xmax)
#         data$xmax <- data$x
#         data <- ddply(data, "ymin", strategy, height = height)
#         data$x <- data$xmax
#         data
#     }
#     else {
#         stop("Neither x nor xmax defined")
#     }
# }
# 
# #' pos_dodgev
# #' 
# #' Vertical dodging internal
# #' 
# #' A hacky adaptation of ggplot's pos_dodge function to be used for vertical collisions.  No warranties on this working except that it looks good for now.
# #' 
# #' The adaptation switch x's to y's and y's to x's and width to height
# #' 
# #' @author Jared P. Lander
# #' @aliases pos_dodgev
# #' @export pos_dodgev
# #' @return Not sure
# #' @param df Refer to ggplot2
# #' @param height Refer to ggplot2
# #' @examples
# #' # none here
# pos_dodgev <- function (df, height) 
# {
#     n <- length(unique(df$group))
#     if (n == 1) 
#         return(df)
#     if (!all(c("ymin", "ymax") %in% names(df))) {
#         df$ymin <- df$y
#         df$ymax <- df$y
#     }
#     d_width <- max(df$ymax - df$ymin)
#     diff <- height - d_width
#     groupidx <- match(df$group, sort(unique(df$group)))
#     df$y <- df$y + height * ((groupidx - 0.5)/n - 0.5)
#     df$ymin <- df$y - d_width/n/2
#     df$ymax <- df$y + d_width/n/2
#     df
# }
# 
# #' position_dodgev
# #' 
# #' Vertical dodging internal
# #' 
# #' A hacky adaptation of ggplot's position_dodge function to be used for vertical collisions.  No warranties on this working except that it looks good for now.
# #' 
# #' No changes were necessary but this is a new geom with a v on the end
# #' 
# #' @author Jared P. Lander
# #' @aliases position_dodgev
# #' @return Not sure
# #' @param width Refer to ggplot2
# #' @param height Refer to ggplot2
# #' @export position_dodgev
# #' @examples
# #' # none here
# position_dodgev <- function (width = NULL, height = NULL) { 
#     #PositionDodgev$new(width = width, height = height)
#     ggplot2::ggproto(NULL, PositionDodgev, width=width, height=height)
# }
# 
# PositionDodgev <- ggplot2::ggproto("dodgev", ggplot2::Position, 
#        required_aes = "y",
#        height = NULL,
#        setup_params = function(self, data) {
#            if (is.null(data$ymin) && is.null(data$ymax) && is.null(self$height)) {
#                warning("Height not defined. Set with `position_dodgev(height = ?)`",
#                        call. = FALSE)
#            }
#            list(height = self$height)
#        },
#        
#        compute_panel = function(data, params, scales) {
#            ggplot2:::collide(data, params$height, "position_dodgev", pos_dodgev, check.width = FALSE)
#        }
# )
