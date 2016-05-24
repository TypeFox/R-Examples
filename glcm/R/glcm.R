# Function to calculate edge for glcm when processing block by block
calc_glcm_edge <- function(shift, window) {
    if ((length(shift) == 2) && is.numeric(shift)) shift <- list(shift)
    if ((!(is.vector(shift) && all(lapply(shift, length) == 2)) &&
         !(is.matrix(shift) && ncol(shift) == 2)) ||
        !(all(floor(unlist(shift)) == unlist(shift)))) {
        stop('shift must be a list of length 2 integer vectors, or a 2 column matrix')
    }
    if (!is.matrix(shift)) {
        shift <- matrix(unlist(shift), ncol=2, byrow=TRUE)
    }
    neg_shifts <- shift[, 2][shift[, 2] < 0]
    pos_shifts <- shift[, 2][shift[, 2] > 0]
    if (length(neg_shifts) == 0) neg_shifts <- 0
    if (length(pos_shifts) == 0) pos_shifts <- 0
    return(c(abs(min(neg_shifts)) + ceiling(window[2] / 2) - 1,
             abs(max(pos_shifts)) + ceiling(window[2] / 2) - 1))
}

#' Image texture measures from grey-level co-occurrence matrices (GLCM)
#'
#' This function supports calculating texture statistics derived from 
#' grey-level co-occurrence matrices (GLCMs). The default textures are 
#' calculated using a 45 degree shift. See Details for other options.
#'
#' The \code{statistics} parameter should be a list, and can include any (one 
#' or more) of the following: 'mean', 'mean_ENVI', 'variance', 'variance_ENVI', 
#' 'homogeneity', 'contrast', 'dissimilarity', 'entropy', 'second_moment', 
#' and/or 'correlation'. By default all of the statistics except for 
#' "mean_ENVI" and "variance_ENVI" will be returned .
#'
#' \code{shift} can be one of:
#' \enumerate{
#'   \item a two element integer vector giving the shift (Q in Gonzalez and 
#'   Woods, 2008), as (number of rows, number of columns).
#'
#'  \item a list of integer vectors of length 2 specifying multiple (row, col) 
#'  shifts over which to calculate the GLCM textures. For example:
#'  \code{shift=list(c(1,1), c(-1,-1))}
#'
#'  \item a matrix with two columns specifying, in rows, multiple (row, col) 
#'  shifts over which to calculate the GLCM textures. For example:
#'  \code{shift=matrix(c(1,1,-1,-1), byrow=TRUE, ncol=2)}
#' }
#'
#' If multiple shifts are supplied, \code{glcm} will calculate each texture 
#' statistic using all the specified shifts, and return the mean value of the 
#' texture for each pixel. To calculate GLCM textures over "all directions" (in 
#' the terminology of commonly used remote sensing software), use: 
#' \code{shift=list(c(0,1), c(1,1), c(1,0), c(1,-1))}. This will calculate the 
#' average GLCM texture using shifts of 0 degrees, 45 degrees, 
#' 90 degrees, and 135 degrees.
#' @export
#' @encoding UTF-8
#' @import Rcpp
#' @importFrom utils stack
#' @usage glcm(x, n_grey = 32, window = c(3, 3), shift = c(1, 1), statistics = 
#' c("mean", "variance", "homogeneity", "contrast", "dissimilarity", "entropy", 
#' "second_moment", "correlation"), min_x=NULL, max_x=NULL, na_opt="any", 
#' na_val=NA, scale_factor=1, asinteger=FALSE)
#' @param x a \code{RasterLayer} or \code{matrix}
#' @param n_grey number of grey levels to use in texture calculation
#' @param window the window size to consider for texture calculation as a two 
#' element integer vector (number of rows, number of columns)
#' @param shift a list or matrix specifying the shift to use. See Details.
#' @param statistics A list of GLCM texture measures to calculate (see 
#' Details).
#' @param min_x minimum value of input \code{RasterLayer} (optional, 
#' \code{glcm} will calculate if not supplied). Useful when running \code{glcm} 
#' over blocks of a raster.
#' @param max_x maximum value of input \code{RasterLayer} (optional, 
#' \code{glcm} will calculate if not supplied). Useful when running \code{glcm} 
#' over blocks of a raster.
#' @param na_opt How to handle NA values in \code{x}. Can be set to "ignore", 
#' "any" or "center". If set to "any", all textures statistics for a given 
#' pixel will be set to NA if there are any NA values in the \code{window} 
#' around that pixel. If set to "center" this will only occur if the center 
#' value is an NA. If set to "ignore", NA values in \code{window} will be 
#' ignored.
#' @param na_val the value to use to fill NA values on edges of \code{x} where 
#' textures cannot be calculated due to the window falling outside of the 
#' image, and as necessary depending on the chosen \code{na_opt}.
#' @param scale_factor factor by which to multiply results.  Useful if rounding 
#' results to integers (see \code{asinteger} argument).
#' @param asinteger whether to round results to nearest integer. Can be used to 
#' save space by saving results as, for example, an 'INT2S' \code{raster}.
#' @return A \code{RasterLayer} or \code{RasterStack} with the requested GLCM 
#' texture measures.
#' @references
#' Lu, D., and M. Batistella. 2005. Exploring TM image texture and its 
#' relationships with biomass estimation in Rondônia, Brazilian Amazon.  Acta 
#' Amazonica 35:249--257.
#'
#' Gonzalez, R. C. 2008. Digital image processing. 3rd ed. Prentice Hall, Upper 
#' Saddle River, N.J, pages 830--836.
#'
#' Haralick, R. M., K. Shanmugam, and I. Dinstein. 1973. Textural features for 
#' image classification. IEEE Transactions on Systems, Man and Cybernetics 
#' SMC-3:610--621.
#'
#' Pratt, W. K. 2007. Digital image processing: PIKS Scientific inside. 4th ed.
#' Wiley-Interscience, Hoboken, N.J pages 540--541, 563--566.
#' @examples
#' \dontrun{
#' require(raster)
#' # Calculate GLCM textures using default 90 degree shift
#' textures_shift1 <- glcm(raster(L5TSR_1986, layer=1))
#' plot(textures_shift1)
#'
#' # Calculate GLCM textures over all directions
#' textures_all_dir <- glcm(raster(L5TSR_1986, layer=1),
#'                          shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)))
#' plot(textures_all_dir)
#' }
glcm <- function(x, n_grey=32, window=c(3, 3), shift=c(1, 1),
                 statistics=c('mean', 'variance', 'homogeneity', 'contrast', 
                              'dissimilarity', 'entropy', 'second_moment', 
                              'correlation'), min_x=NULL, max_x=NULL, 
                 na_opt='any', na_val=NA, scale_factor=1, asinteger=FALSE) {
    if (!inherits(x, 'RasterLayer') &&
        !(inherits(x, 'matrix') && (length(dim(x)) == 2))) {
        stop('x must be a RasterLayer or two-dimensional matrix')
    }

    if (length(window) != 2 || !all(floor(window) == window)) {
        stop('window must be integer vector of length 2')
    }
    # Convert a length 2 vector to a 2 element list (for
    # backwards compatibility)
    if ((length(shift) == 2) && is.numeric(shift)) shift <- list(shift)
    if ((!(is.vector(shift) && all(lapply(shift, length) == 2)) &&
         !(is.matrix(shift) && ncol(shift) == 2)) ||
        !(all(floor(unlist(shift)) == unlist(shift)))) {
        stop('shift must be a list of length 2 integer vectors, or a 2 column matrix')
    }
    if (!is.matrix(shift)) {
        shift <- matrix(unlist(shift), ncol=2, byrow=TRUE)
    }
    min_rast_rows <- window[1] + abs(max(shift[, 1]))
    min_rast_cols <- window[2] + abs(max(shift[, 2]))
    if ((window[1] %% 2 == 0) || (window[2] %% 2 == 0)) {
        stop('both elements of window must be odd')
    } else if (min_rast_rows > nrow(x)) {
        stop("window[1] + the maximum x shift value must be less than nrow(x)")
    } else if (min_rast_cols > ncol(x)) {
        stop("window[2] + the maximum y shift value must be less than ncol(x)")
    } else if (!inherits(statistics, 'character')) {
        stop('statistics must be a character vector')
    }
    avail_stats <- c('mean', 'mean_ENVI', 'variance', 'variance_ENVI', 
                     'homogeneity', 'contrast', 'dissimilarity', 'entropy', 
                     'second_moment', 'correlation')
    stat_check <- unlist(lapply(statistics, function(stat) stat %in% avail_stats))
    if (sum(stat_check) != length(stat_check)) {
        stop(paste('invalid statistic(s):',
                   paste(statistics[!stat_check], collapse=', ')))
    }
    if (!(na_opt %in% c('any', 'center', 'ignore'))) {
        stop('na_opt must be "any", "center", or "ignore"')
    }

    if (inherits(x, 'RasterLayer')) {
        if (is.null(min_x)) min_x <- raster::cellStats(x, 'min')
        if (is.null(max_x)) max_x <- raster::cellStats(x, 'max')

        if (raster::canProcessInMemory(x, length(statistics) + 2)) {
            x_cut <- raster::cut(x, breaks=seq(min_x, max_x,
                                               length.out=(n_grey + 1)),
                                 include.lowest=TRUE, right=FALSE)
            x_cut <- raster::as.matrix(x_cut)
            textures <- calc_texture(x_cut, n_grey, window, shift, statistics, na_opt, 
                                     na_val)
        } else {
            edge <- calc_glcm_edge(shift, window)
    
            min_block_rows <- window[1] + edge[1] + edge[2]
            bs <- raster::blockSize(x, minrows=min_block_rows, minblocks=1)
            # Handle case of very small last blocks by combining the two final 
            # blocks when the last block is very small
            if ((bs$n > 1) & bs$nrows[length(bs$nrows)] < min_block_rows) {
                bs$nrows[length(bs$nrows) - 1] <- sum(bs$nrows[(length(bs$nrows) - 1):length(bs$nrows)])
                bs$nrows <- bs$nrows[-length(bs$nrows)]
                bs$row <- bs$row[-length(bs$row)]
                bs$n <- length(bs$row)
            }
            n_blocks <- bs$n

            # bs_mod is the blocksize that will contain blocks that have been expanded 
            # to avoid edge effects
            bs_mod <- bs
            if (n_blocks > 1) {
                # Expand blocks to account for edge effects on the top:
                bs_mod$row[2:n_blocks] <- bs_mod$row[2:n_blocks] - edge[1]
                # Need to read additional rows from these blocks to avoid an offset
                bs_mod$nrows[2:n_blocks] <- bs_mod$nrows[2:n_blocks] + edge[1]
                # Read more bottom rows to account for bottom edge effects
                bs_mod$nrows[1:(n_blocks - 1)] <- bs_mod$nrows[1:(n_blocks - 1)] + edge[2]
            }

            # TODO: Should detect this and fix automatically
            if (any(bs_mod$row < 1)) {
                stop("underflow: cannot read without edge effects - report to package author")
            } else if (any((bs_mod$nrows + bs_mod$row - 1) > nrow(x))) {
                stop("overflow: cannot read without edge effects - report to package author")
            } else if ((n_blocks > 1) & any(bs_mod$nrows < min_block_rows)) {
                stop("cannot read without edge effects - report to package author")
            }
            
            started_writes <- FALSE
            for (block_num in 1:n_blocks) {
                this_block <- raster::getValues(x, row=bs_mod$row[block_num], 
                                        nrows=bs_mod$nrows[block_num],
                                        format='matrix')
                x_cut <- matrix(findInterval(this_block, seq(min_x, max_x, length.out=(n_grey + 1)),
                                             all.inside=TRUE), nrow=nrow(this_block))
                out_block <- calc_texture(x_cut, n_grey, window, shift, 
                                          statistics, na_opt, na_val)
                layer_names <- dimnames(out_block)[[3]]
                # Drop the padding added to top to avoid edge effects, unless we are 
                # really on the top of the image, where top edge effects cannot be 
                # avoided
                if ((block_num != 1) && (edge[1] > 0)) {
                    out_block <- out_block[-(1:edge[1]), , ]
                    # The below line is needed to maintain a 3 dimensional array, 
                    # even when an n x m x 1 array is returned from 
                    # calc_texture_full_image because a single statistic was chosen. 
                    # Without the below line, removing a row will coerce the 3d array 
                    # to a 2d matrix, and the bottom padding removal will fail as it 
                    # references a 3d matrix).
                    if (length(dim(out_block)) < 3) dim(out_block) <- c(dim(out_block), 1)
                }
                # Drop the padding added to bottom to avoid edge effects, unless we are 
                # really on the bottom of the image, where bottom edge effects cannot 
                # be avoided
                if ((block_num != n_blocks) && (edge[2] > 0)) {
                    out_block <- out_block[-((nrow(out_block)-edge[2]+1):nrow(out_block)), , ]
                    if (length(dim(out_block)) < 3) dim(out_block) <- c(dim(out_block), 1)
                }
                if (!started_writes) {
                    # Setup an output raster with number of layers equal to the number 
                    # of layers in out_block, and extent/resolution equal to extent and 
                    # resolution of x
                    if (dim(out_block)[3] == 1) {
                        textures <- raster::raster(x)
                    } else {
                        textures <- raster::brick(x, nl=dim(out_block)[3], values=FALSE)
                    }
                    textures <- raster::writeStart(textures, filename=raster::rasterTmpFile())
                    names(textures) <- layer_names
                    started_writes <- TRUE
                }
                # To write to a RasterBrick the out_block needs to be structured as 
                # a 2-d matrix with bands in columns and columns as row-major vectors
                if (dim(out_block)[3] == 1) {
                    out_block <- aperm(out_block, c(3, 2, 1))
                    out_block <- matrix(out_block, ncol=nrow(out_block))
                } else {
                    out_block <- aperm(out_block, c(3, 2, 1))
                    out_block <- matrix(out_block, ncol=nrow(out_block), byrow=TRUE)
                }
                textures <- raster::writeValues(textures, out_block, bs$row[block_num])
            }
            textures <- raster::writeStop(textures)
        }
        if (!inherits(textures, c('RasterLayer', 'RasterStack', 
                                  'RasterBrick'))) {
            if (dim(textures)[3] > 1) {
                textures <- raster::stack(apply(textures, 3, raster::raster, 
                                                template=x))
            } else {
                textures <- raster::raster(textures[, , 1], template=x)
            }
        }
        names(textures) <- paste('glcm', statistics, sep='_')
    } else if (inherits(x, 'matrix')) {
        if (is.null(min_x)) min_x <- min(x)
        if (is.null(max_x)) max_x <- max(x)
        x_cut <- matrix(findInterval(x, seq(min_x, max_x, length.out=(n_grey + 1)),
                                     all.inside=TRUE), nrow=nrow(x))
        textures <- calc_texture(x_cut, n_grey, window, shift, statistics, na_opt, 
                                 na_val)
        dimnames(textures) <- list(NULL, NULL, paste('glcm', statistics, sep='_'))
    } else {
        stop('x must be a RasterLayer or two-dimensional matrix')
    }

    if (scale_factor != 1) {
        textures <- textures * scale_factor
    }

    if (asinteger) {
        textures <- round(textures)
    }

    return(textures)
}
