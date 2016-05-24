# Detect and prevent collisions.
# Powers dodging, stacking and filling.
collidev <- function(data, height = NULL, name, strategy, check.height = TRUE) {
    # Determine height
    if (!is.null(height)) {
        # height set manually
        if (!(all(c("ymin", "ymax") %in% names(data)))) {
            data$ymin <- data$y - height / 2
            data$ymax <- data$y + height / 2
        }
    } else {
        if (!(all(c("ymin", "ymax") %in% names(data)))) {
            data$ymin <- data$y
            data$ymax <- data$y
        }
        
        # height determined from data, must be floating point constant
        heights <- unique(data$ymax - data$ymin)
        heights <- heights[!is.na(heights)]
        
        #   # Suppress warning message since it's not reliable
        #     if (!zero_range(range(heights))) {
        #       warning(name, " requires constant height: output may be incorrect",
        #         call. = FALSE)
        #     }
        height <- heights[1]
    }
    
    # Reorder by x position, relying on stable sort to preserve existing
    # ordering, which may be by group or order.
    data <- data[order(data$ymin), ]
    
    # Check for overlap
    intervals <- as.numeric(t(unique(data[c("ymin", "ymax")])))
    intervals <- intervals[!is.na(intervals)]
    
    if (length(unique(intervals)) > 1 & any(diff(scale(intervals)) < -1e-6)) {
        warning(name, " requires non-overlapping y intervals", call. = FALSE)
        # This is where the algorithm from [L. Wilkinson. Dot plots.
        # The American Statistician, 1999.] should be used
    }
    
    if (!is.null(data$xmax)) {
        plyr::ddply(data, "ymin", strategy, height = height)
    } else if (!is.null(data$x)) {
        data$xmax <- data$x
        data <- plyr::ddply(data, "ymin", strategy, height = height)
        data$x <- data$xmax
        data
    } else {
        stop("Neither x nor xmax defined")
    }
}

# Stack overlapping intervals.
# Assumes that each set has the same horizontal position
pos_stackv <- function(df, height) {
    if (nrow(df) == 1) return(df)
    
    n <- nrow(df) + 1
    x <- ifelse(is.na(df$x), 0, df$x)
    if (all(is.na(df$y))) {
        heights <- rep(NA, n)
    } else {
        heights <- c(0, cumsum(x))
    }
    
    df$xmin <- heights[-n]
    df$xmax <- heights[-1]
    df$x <- df$xmax
    df
}

# Stack overlapping intervals and set height to 1.
# Assumes that each set has the same horizontal position.
pos_fillv <- function(df, height) {
    stacked <- pos_stackv(df, height)
    stacked$xmin <- stacked$xmin / max(stacked$xmax)
    stacked$xmax <- stacked$xmax / max(stacked$xmax)
    stacked$x <- stacked$xmax
    stacked
}

# Dodge overlapping interval.
# Assumes that each set has the same horizontal position.
pos_dodgev <- function(df, height) {
    n <- length(unique(df$group))
    if (n == 1) return(df)
    
    if (!all(c("ymin", "ymax") %in% names(df))) {
        df$ymin <- df$y
        df$ymax <- df$y
    }
    
    d_height <- max(df$ymax - df$ymin)
    
    # df <- data.frame(n = c(2:5, 10, 26), div = c(4, 3, 2.666666,  2.5, 2.2, 2.1))
    # ggplot(df, aes(n, div)) + geom_point()
    
    # Have a new group index from 1 to number of groups.
    # This might be needed if the group numbers in this set don't include all of 1:n
    groupidy <- match(df$group, sort(unique(df$group)))
    
    # Find the center for each group, then use that to calculate xmin and xmax
    df$y <- df$y + height * ((groupidy - 0.5) / n - .5)
    df$ymin <- df$y - d_height / n / 2
    df$ymax <- df$y + d_height / n / 2
    
    df
}


#' Adjust position by dodging overlaps to the side.
#'
#' @inheritParams ggplot2::position_identity
#' @param height Dodging height, when different to the height of the individual
#'   elements. This is useful when you want to align narrow geoms with wider
#'   geoms. See the examples for a use case.
#' @family position adjustments
#' @export
#' @examples
#' ggplot(mtcars, aes(factor(cyl), fill = factor(vs))) +
#'   geom_bar(position = "dodge")
#' \donttest{
#' ggplot(diamonds, aes(price, fill = cut)) +
#'   geom_histogram(position="dodge")
#' # see ?geom_boxplot and ?geom_bar for more examples
#'
#' # To dodge items with different heights, you need to be explicit
#' df <- data.frame(x=c("a","a","b","b"), y=2:5, g = rep(1:2, 2))
#' p <- ggplot(df, aes(x, y, group = g)) +
#'   geom_bar(
#'     stat = "identity", position = "dodge",
#'     fill = "grey50", colour = "black"
#'   )
#' p
#'
#' # A line range has no height:
#' p + geom_linerange(aes(ymin = y-1, ymax = y+1), position = "dodge")
#' # You need to explicitly specify the height for dodging
#' p + geom_linerange(aes(ymin = y-1, ymax = y+1),
#'   position = position_dodge(height = 0.9))
#'
#' # Similarly with error bars:
#' p + geom_errorbar(aes(ymin = y-1, ymax = y+1), height = 0.2,
#'   position = "dodge")
#' p + geom_errorbar(aes(ymin = y-1, ymax = y+1, height = 0.2),
#'   position = position_dodge(height = 0.90))
#' }
position_dodgev <- function(height = NULL) {
    ggproto(NULL, PositionDodgeV, height = height)
}



PositionDodgeV <- ggproto("PositionDodgeV", Position,
                          required_aes = "y",
                          height = NULL,
                          setup_params = function(self, data) {
                              if (is.null(data$ymin) && is.null(data$ymax) && is.null(self$height)) {
                                  warning("height not defined. Set with `position_dodgev(height = ?)`",
                                          call. = FALSE)
                              }
                              list(height = self$height)
                          },
                          
                          compute_panel = function(data, params, scales) {
                              collidev(data, params$height, "position_dodgev", pos_dodgev, check.height = FALSE)
                          }
)
