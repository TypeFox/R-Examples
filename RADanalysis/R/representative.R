#'  Computes representative normalized RAD of a group of normalized RADs.
#' @param norm_rad A matrix which contains the normalized RADs (samples in rows).
#' @param sample_ids Vector of row numbers of the desired group,
#'     from which a representative RAD is going to be produced.
#' @param plot A logical. If \code{TRUE}, plots the repRAD. The plot would be added to the previous plot.
#' @param min_rank The minimum rank to be considered for making repRADs.
#' @param confidence Confidence interval of plotted repRAD. Default is 0.9.
#' @param with_conf A logical. If \code{TRUE}, plots the confidence interval in addition to repRAD. Only works when \code{plot} is \code{TRUE}.
#' @param ... Other graphical parameters to use for plotting. This function uses internally the
#'     functions \code{lines} and \code{polygon} to plot.

#' @return A list of following parameters:
#' @return $average: Contains a vector of length equal to the columns of \code{norm_rad}. This in the representative normalized RAD which is
#'     the average of normalized RAD of the group.
#' @return $population_sd: A vector of length equal to the columns of \code{norm_rad} which contains the standard deviation
#'     for each rank.
#' @return $standard_error: A vector of length equal to the columns of \code{norm_rad} which contains the standard deviation
#'     of the mean for each rank. This vector is the result of \code{population_sd / sqrt(n)},
#'     when n is the number of members of the group (length of \code{sample_ids}).
#' @return If \code{plot = TRUE}, plot of the repRAD is produced and would be added to the previous plot.
#' @return If \code{with_conf = TRUE}, confidence interval would be added to the repRAD plot.
#' @seealso \code{\link{RADnormalization}} for normalize an abundance vector. This function return more details compared to \code{\link{RADnormalization_matrix}},
#'          \code{\link{RADnormalization_matrix}} for normalize an entire otutable,
#'          \code{\link{representative_point}} for study the representative of groups of samples in a multi-dimensional scaling plot,

#' @examples
#' line_cols <- c("green","red","blue")
#' sample_classes <- c(1,1,1,1,2,2,3,3,1,1,2,3,3,1,1,2,3,3)
#' maxrank <- 400
#' data("gut_nrads")
#' nrads <- gut_nrads
#' nrads <- nrads$norm_matrix
#'
#' #plot nrads
#' plot(1e10,xlim = c(1,maxrank),ylim = c(2e-5,1),log="xy",
#'      xlab = "rank",ylab = "abundance",cex.lab = 1.5,axes = FALSE)
#' sfsmisc::eaxis(side = 1,at = c(1,10,100,1000,10000))
#' sfsmisc::eaxis(side = 2,at = c(1e-4,1e-3,1e-2,1e-1,1),las = 0)
#' for(i in 1:nrow(nrads)){
#'     points(nrads[i,],type = 'l',col = line_cols[sample_classes[i]],lwd = 0.8)
#' }
#' #plot confidence intervals of representative nrads
#' a <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == 1),
#'                       plot = TRUE,confidence = 0.9,with_conf = TRUE,
#'                       col = scales::alpha(line_cols[1],0.5),border = NA)
#' a <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == 2),
#'                       plot = TRUE,confidence = 0.9,with_conf = TRUE,
#'                       col = scales::alpha(line_cols[2],0.5),border = NA)
#' a <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == 3),
#'                       plot = TRUE,confidence = 0.9,with_conf = TRUE,
#'                       col = scales::alpha(line_cols[3],0.5),border = NA)
#' #plot representative nrads
#' a <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == 1),
#'                       plot = TRUE,with_conf = FALSE,
#'                       col = scales::alpha(line_cols[1],0.8),lwd = 4)
#' a <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == 2),
#'                       plot = TRUE,with_conf = FALSE,
#'                       col = scales::alpha(line_cols[2],0.8),lwd = 4)
#' a <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == 3),
#'                       plot = TRUE,with_conf = FALSE,
#'                       col = scales::alpha(line_cols[3],0.8),lwd = 4)
#' legend("bottomleft",bty = "n",legend = c("pre Cp","under Cp","post Cp"),
#' col = line_cols,lwd = 3)
#'
#' @export representative_RAD
representative_RAD <- function(norm_rad,sample_ids = NULL,plot = F,min_rank = 1,confidence = 0.95,with_conf = TRUE,...){
    if (is.null(sample_ids)) sample_ids <- 1:(dim(norm_rad)[1])
    if (length(sample_ids) == 0) return(NULL)
    if (length(sample_ids) == 1) sample_ids <- c(sample_ids,sample_ids)

    s_data <- norm_rad[sample_ids,]
    av <- apply(s_data, 2, mean)
    population_sd <- apply(s_data, 2, stats::sd)
    confi <- confidence
    a <- stats::qnorm(confi + (1-confi)/2)
    se <- a * population_sd / sqrt(length(sample_ids))
    if(plot){
        x <- min_rank:(dim(norm_rad)[2]+min_rank-1)
        graphics::lines(x = x,y = av,...)
        if (with_conf) {
            if (length(sample_ids) == 2) {
                graphics::polygon(x = c(x , rev(x)), y = c(apply(s_data, 2, max),rev(apply(s_data, 2, min))),...)
            }else{
                graphics::polygon(x = c(x , rev(x)), y = c((av + se),rev(av-se)),...)
            }
        }
    }
    results <- list(average = av,population_sd = population_sd,standard_error = se)
    return(results)
}






#'  Computes representative point based on the coordinates of points which are in
#'     the same group.
#' @param input A matrix which contains the coordinates of samples. Usually this is the
#'     result of ordination of normalized RADs using multi-dimensional scaling (\code{cmdscale}).
#'     In the input matrix each row contains vector of coordinates of one sample.
#' @param ids Vector of row numbers of the desired group,
#'     from which a representative point is going to be represented
#' @param plot A logical. If \code{TRUE}, shows the representative points on the previous plot.
#' @param coord_names A vector which contains the coordintes number that should be used to create representative point.
#'     Default is \code{c(1,2)}.
#' @param standard_error_mean A logical. If \code{TRUE}, uses the standard error of the mean
#'     and plot it with representative points. It works only if \code{plot = TRUE}.
#' @param ... other graphical parameters to use for plotting. This function uses
#'     internally the functions \code{points} and \code{arrows} to plot.
#'
#' @return A list of following parameters:
#' @return $mean: Contains the average of points. A vector with the length of coordinates
#'     used for computing the average. These coordinates are preset in \code{coord_names}.
#' @return $sd: A vector with a length similar to \code{mean} which contains the
#'     standard deviation for each coordinate.
#' @return $mean_standard_error: A vector with a length similar to \code{mean} which
#'     contain the standard deviation of the mean for each coordinate. This vector is the result of \code{sd / sqrt(n)},
#'     when n is the number of members of the group (length of \code{sample_ids}).
#' @return If \code{plot = TRUE}, representative points would be added to the previous plot.
#' @return If \code{standard_error_mean = TRUE}, the standard error of the mean would be added to the representative points.
#' @seealso \code{\link{RADnormalization}} for normalize an abundance vector. This function return more details compared to \code{\link{RADnormalization_matrix}},
#'          \code{\link{RADnormalization_matrix}} for normalize an entire otutable,
#'          \code{\link{representative_RAD}} for study the representative of group of norm rads.

#' @examples
#' line_cols <- c("green","red","blue")
#' sample_classes <- c(1,1,1,1,2,2,3,3,1,1,2,3,3,1,1,2,3,3)
#' maxrank <- 400

#' data("gut_nrads")
#' nrads <- gut_nrads
#' nrads <- nrads$norm_matrix
#'
#' #distance matrix using manhattan distance
#' d <- dist(x = nrads,method = "manhattan")
#' #ordination using classical multi-dimensional scaling
#' mds <- cmdscale(d = d,k = 5,eig = TRUE)
#'
#' #plot the points
#' plot(mds$points,xlab = "First coordinate",ylab = "Second coordinate",pch = 19,cex =1,
#'     col = line_cols[sample_classes],
#'     main = "MDS plot with representative points \n of each group and error bars")
#'
#' #add the representative points wit erorr bar to the previous plot
#' a <- representative_point(input = mds$points,ids = which(sample_classes == 1),
#'     col = scales::alpha(line_cols[1],0.5),
#'     plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
#' a <- representative_point(input = mds$points,ids = which(sample_classes == 2),
#'     col = scales::alpha(line_cols[2],0.5),
#'     plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
#' a <- representative_point(input = mds$points,ids = which(sample_classes == 3),
#'     col = scales::alpha(line_cols[3],0.5),
#'     plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
#'
#' legend("bottomleft",bty = "n",legend = c("pre Cp","under Cp","post Cp"),
#'     col = line_cols,pch = 19)
#'
#' @export representative_point
representative_point <- function(input,ids = NULL,coord_names = c(1,2),standard_error_mean = TRUE,plot = FALSE,...) {
    if (is.null(ids)) ids <- 1:(dim(input)[1])
    if (length(ids) == 0) return(NULL)
    if (length(ids) == 1) ids <- c(ids,ids)

    subinput <- input[ids,coord_names]
    m <- apply(subinput, 2, mean)
    names(m) <- paste("coord_",coord_names,sep = "")
    s <- apply(subinput, 2, stats::sd)
    names(s) <- paste("coord_",coord_names,sep = "")
    e <- s / sqrt(length(ids))
    names(e) <- paste("coord_",coord_names,sep = "")

    if(plot) {
        if(length(coord_names) == 2) {
            if(standard_error_mean) {
                graphics::points(x = m[1],y = m[2],...)
                graphics::arrows(x0 = m[1], y0 = m[2], x1 = m[1] + e[1], y1 = m[2], length = 0.1, angle = 90)
                graphics::arrows(x0 = m[1], y0 = m[2], x1 = m[1] - e[1], y1 = m[2], length = 0.1, angle = 90)
                graphics::arrows(x0 = m[1], y0 = m[2], x1 = m[1], y1 = m[2] + e[2], length = 0.1, angle = 90)
                graphics::arrows(x0 = m[1], y0 = m[2], x1 = m[1], y1 = m[2] - e[2], length = 0.1, angle = 90)
            }else{
                graphics::points(x = m[1],y = m[2],...)
                graphics::arrows(x0 = m[1], y0 = m[2], x1 = m[1] + s[1], y1 = m[2], length = 0.1, angle = 90)
                graphics::arrows(x0 = m[1], y0 = m[2], x1 = m[1] - s[1], y1 = m[2], length = 0.1, angle = 90)
                graphics::arrows(x0 = m[1], y0 = m[2], x1 = m[1], y1 = m[2] + s[2], length = 0.1, angle = 90)
                graphics::arrows(x0 = m[1], y0 = m[2], x1 = m[1], y1 = m[2] - s[2], length = 0.1, angle = 90)
            }

        }else{
            print("Only two dimensional plot is accepted")
        }
    }
    result <- list(m,s,e)
    names(result) <- c("mean","sd","mean_standard_error")
    return(result)
}



