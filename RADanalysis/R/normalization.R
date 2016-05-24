#'  Normalizes an abundance vector to the desired number of ranks.
#' @param input A vector which contains the abundance values (an abundance vector).
#' @param max_rank The desired rank to which this method normalizes the input.
#' @param average_over Number of times, a normalized RAD is created and averaged to
#'     produce the result.
#' @param min_rank The minimum rank to which this method normalizes the input.
#' @param labels A logical. If \code{TRUE} the label of each rank (ids in the input vector)
#'     will be returned.
#' @param count_data A logical. \code{TRUE} means that the input vector contains counts
#'     (integer values) otherwise contains the relative abundances. In the current version only counts are accepted.
#' @param method Sets the stop criterion for normalization. This should be one of
#'     "lowerlimit", "middle" or
#'     "upperlimit". Method affects the final result.
#'     lowerlimit: Samples from species pool one by one, until reaches max_rank.
#'     middle: Samples from species pool with random size until the sampled vector has desired ranks (max_rank).
#'     upperlimit: Removes from species pool one by one, until reaches max_rank.
#'
#' @return A list of following items:
#' @return $norm_rad: Normalized RAD sum up to 1. If \code{labels = TRUE}, it would also contain the labels.
#' @return $norm_rad_count:  A matrix of \code{average_over} rows and \code{max_rank} columns.
#'     Each row contains one normalized RAD. These normalized RADs are averaged and sum up to 1 in order to make norm_rad
#' @return $norm_rad_mean_sd: Standard deviation of the mean for all the ranks in \code{norm_rad}.
#'     This vector is created using the values in \code{norm_rad_count}
#' @return $inputs: A list which contains inputs used for creating normalized RADs.
#' @examples
#' data("gut_otu_table")
#' rads <- gut_otu_table
#' original_rad <- sort(rads[1,],decreasing = TRUE)
#' #removing zeros
#' original_rad <- original_rad[original_rad > 0]
#' plot(original_rad,ylim = c(1,max(original_rad)),log = "xy", xlab = "Rank",ylab = "Abundance",
#'      main = "RAD of first sample",pch = 19,type = "b",cex = 0.5)


#' print(paste("number of ranks present in the original rad is:",length(original_rad)))

#' norm_rad <- RADnormalization(input = rads[1,],max_rank = 500,average_over = 50)
#' points(x = norm_rad$norm_rad * sum(norm_rad$norm_rad_count[1,]) ,pch = 19,cex = 1, type = "l",
#'        col = "blue",lwd = 4)
#' points(x = norm_rad$norm_rad_count[1,],pch = 19,cex = 1, type = "l",col = "red",lwd = 3)
#' points(x = norm_rad$norm_rad_count[2,],pch = 19,cex = 1, type = "l",col = "green",lwd = 3)
#' legend("bottomleft",legend = c("Original RAD","possible norm rad","possible norm rad",
#'                                paste("nrad averaged over 50 realizations, times",
#'                                sum(norm_rad$norm_rad_count[1,]))),
#'        col = c("black","red","green","blue"),lwd = 2,bty = "n")
#' @seealso \code{\link{RADnormalization_matrix}} for normalize an entire otutable,
#'          \code{\link{representative_point}} for study the representative of groups of samples in a multi-dimensional scaling plot,
#'          \code{\link{representative_RAD}} for study the representative of group of norm rads.
#' @export RADnormalization

RADnormalization <- function(input,max_rank,average_over = 1,min_rank = 1,labels = FALSE,count_data = TRUE,method = "upperlimit"){
    ### methods: lowerlimit, middle, upperlimit

    ## known bug: for the very special case of something like pool = c(1,2,5,6,3,2,1,5,4,2,1,5,125) (last one a new one) when we want to norm to original richness it does not work
    if(!(method %in% c("upperlimit","lowerlimit","middle"))) stop("method must one of: lowerlimit, middle or upperlimit")
    input <- as.vector(input)
    if(labels == T){if(average_over>1){average_over <- 1;print(paste("average_over changed to 1 because labels set to TRUE."));}}
    #input must be a vector of integers

    if(!count_data) {
        stop("This version only accepts count data (integer values)")
#        temp <- RADnormalization_relative_numbers(input = input,max_rank = max_rank,average_over = average_over,min_rank = min_rank,labels = labels,count_data = count_data)
#        return(temp)
    }

    #  if(count_data){
    if (!isTRUE(all(input == floor(input)))) stop("'input' must only contain integer values")
    non_zero_ids <- which(input>0)
    input_without_zeros <- input[non_zero_ids]
    individual.pool <- non_zero_ids[rep(1:length(non_zero_ids),input_without_zeros)]

    if((length(non_zero_ids) == max_rank) && (method == "upperlimit")){
        if(labels == T){
            s_pool <- individual.pool
            new.otutable <- sort(table(s_pool),decreasing = T)
            norm_rad_count <- cbind(as.integer(rownames(new.otutable)),as.vector(new.otutable))
            colnames(norm_rad_count) <- c("label_ids","Abundance")
            Abundance <- norm_rad_count[min_rank:max_rank,2] / sum(norm_rad_count[min_rank:max_rank,2])
            norm_rad <- cbind(norm_rad_count[,1],Abundance)
            colnames(norm_rad) <- c("label_ids","Abundance")
            inputs <- list(input,min_rank,max_rank,average_over,labels,method)
            names(inputs) <- c("input","min_rank","max_rank","average_over","labels","method")
            out <- list(norm_rad,norm_rad_count,inputs)
            names(out) <- c("norm_rad","norm_rad_count","inputs")
            return(out)

        }else{
            norm_rad_count <- sort(input,decreasing = T)[1:max_rank]
            norm_rad <- norm_rad_count[min_rank:max_rank] / sum(norm_rad_count[min_rank:max_rank])
            norm_rad_mean_sd <- NA
            inputs <- list(input,min_rank,max_rank,average_over,labels,method)
            names(inputs) <- c("input","min_rank","max_rank","average_over","labels","method")
            out <- list(norm_rad,norm_rad_count,norm_rad_mean_sd,inputs)
            names(out) <- c("norm_rad","norm_rad_count","norm_rad_mean_sd","inputs")
            return(out)
        }
    }

    #Normalize to max rank
    norm_rad_count <- c()
    s_pool <- c()


    for(i in 1:average_over){
        s_pool <- sample(individual.pool)
        l_pool <- length(individual.pool)
        mark0 <- 1
        mark1 <- l_pool
        mark <- (mark0+mark1) %/% 2
        current_richness <- length(unique(s_pool[1:mark]))

        ### MIDDLE
        if(method == "middle"){
            while(current_richness != max_rank){
                if(current_richness < max_rank){
                    mark0 <- mark
                    mark  <- ceiling((mark+mark1) / 2)
                }else{
                    mark1 <- mark
                    mark  <- (mark+mark0) %/% 2
                }
                current_richness <- length(unique(s_pool[1:mark]))
            }
        }
        ### UPPERLIMIT
        if(method == "upperlimit"){
            continue <- T
            while(continue){
                if(current_richness <= max_rank){
                    mark0 <- mark
                    mark  <- (mark+mark1) %/% 2
                    if(mark0 == mark) continue <- F
                }else{
                    mark1 <- mark
                    mark  <- (mark+mark0) %/% 2
                }
                current_richness <- length(unique(s_pool[1:mark]))
            }
        }
        ### LOWERLIMIT
        if(method == "lowerlimit"){
            continue <- T
            while(continue){
                if(current_richness < max_rank){
                    mark0 <- mark
                    mark  <- ceiling((mark+mark1) / 2)
                }else{
                    mark1 <- mark
                    mark  <- ceiling((mark+mark0) / 2)
                    if(mark == mark1) continue <- F
                }
                current_richness <- length(unique(s_pool[1:mark]))
            }
        }
        ### COMPUTE NORM RAD
        norm_pool <- s_pool[1:mark]
        norm_table <- table(norm_pool)
        abundance <- sort(as.vector(norm_table),decreasing = T)
        norm_rad_count <- rbind(norm_rad_count,abundance)
    }


    #prepare the result and return them
    if(labels == T){
        new.otutable <- sort(table(s_pool),decreasing = T)
        norm_rad_count <- cbind(as.integer(rownames(new.otutable)),as.vector(new.otutable))
        colnames(norm_rad_count) <- c("label_ids","Abundance")
        Abundance <- norm_rad_count[min_rank:max_rank,2] / sum(norm_rad_count[min_rank:max_rank,2])
        norm_rad <- cbind(norm_rad_count[,1],Abundance)
        colnames(norm_rad) <- c("label_ids","Abundance")
        inputs <- list(input,min_rank,max_rank,average_over,labels,method)
        names(inputs) <- c("input","min_rank","max_rank","average_over","labels","method")
        out <- list(norm_rad,norm_rad_count,inputs)
        names(out) <- c("norm_rad","norm_rad_count","inputs")
        return(out)
    }else{
        if(average_over != 1){
            norm_rad <- colMeans(norm_rad_count[,min_rank:max_rank] / rowSums(norm_rad_count[,min_rank:max_rank]) )
        }else{
            norm_rad <- norm_rad_count[,min_rank:max_rank] / sum(norm_rad_count[,min_rank:max_rank])
        }
        #    norm_rad <- colSums(norm_rad_count[,min_rank:max_rank]) / sum(norm_rad_count[,min_rank:max_rank])
        if(average_over != 1) {
            norm_rad_mean_sd <- apply((norm_rad_count[,min_rank:max_rank] / rowSums(norm_rad_count[,min_rank:max_rank]) ), 2, stats::sd) / sqrt(average_over)
        }else{
            norm_rad_mean_sd <- NA
        }
        inputs <- list(input,min_rank,max_rank,average_over,labels,method)
        names(inputs) <- c("input","min_rank","max_rank","average_over","labels","method")
        out <- list(norm_rad,norm_rad_count,norm_rad_mean_sd,inputs)
        names(out) <- c("norm_rad","norm_rad_count","norm_rad_mean_sd","inputs")
        return(out)
    }
}







#'  Normalizes an abundance table to the desired number of ranks
#' @param input A vector or matrix which contains the abundance values (an abundance table).
#' @param max_rank The desired rank to which this method normalizes the input.
#' @param average_over Number of times, a normalized RAD is created and averaged to produce the result.
#' @param min_rank The minimum rank to which this method normalizes the input.
#' @param labels A logical. If \code{TRUE} the label of each rank (ids in the input vector) will be returned.
#' @param count_data A logical. \code{TRUE} means that the input vector contains counts (integer values) otherwise contains the relative abundances. In the current version only counts are accepted.
#' @param sample_in_row A logical. \code{TRUE} means that the abundance vector of samples are represented in rows otherwise in columns.
#' @param method Sets the stop criterion for normalization. This should be one of
#'     "lowerlimit", "middle" or
#'     "upperlimit". Method affects the final result.
#'     lowerlimit: Samples from species pool one by one, until reaches max_rank.
#'     middle: Samples from species pool with random size until the sampled vector has desired ranks (max_rank).
#'     upperlimit: Removes from species pool one by one, until reaches max_rank.
#' @param verbose A logical. If \code{TRUE}, prints the progress in percent in console.
#'
#' @return A list of following items:
#' @return $norm_matrix A matrix which contains normalized RADs sum up to 1. If \code{labels = TRUE}, it would also contain the labels.
#' @return $inputs A list which contains inputs used for creating normalized RADs. It does not contain \code{input} because it could be very big.
#' @examples
#' data("gut_otu_table")
#' rads <- gut_otu_table
#' #plot original rads
#' line_cols <- c("green","red","blue")
#' sample_classes <- c(1,1,1,1,2,2,3,3,1,1,2,3,3,1,1,2,3,3)
#' plot(1,xlim = c(1,2000),ylim = c(1,20000),col = "white",log  = "xy",
#'      xlab = "Rank",ylab = "Abundance",main = "Original RADs from antibiotic data set")
#' for(i in 1:nrow(rads)){
#'     temp <- sort(rads[i,],decreasing = TRUE)
#'     temp <- temp[temp>0]
#'     lines(x = temp,lwd = 2,col = line_cols[sample_classes[i]])
#' }
#' legend("bottomleft",bty = "n",legend = c("pre Cp","under Cp","post Cp"),col = line_cols,lwd = 3)
#'
#'
#' nrads <- RADnormalization_matrix(input = rads,max_rank = 400,average_over = 20,sample_in_row = TRUE)
#' nrads <- nrads$norm_matrix
#'
#' plot(1,xlim = c(1,400),ylim = c(4e-5,1),col = "white",log  = "xy",
#'      xlab = "Rank",ylab = "Abundance",
#'      main = "NRADs from antibiotic data set with R = 400 \n with average_over = 20")
#' for(i in 1:nrow(nrads)){
#'     lines(x = nrads[i,],lwd = 2,col = line_cols[sample_classes[i]])
#' }
#' legend("bottomleft",bty = "n",legend = c("pre Cp","under Cp","post Cp"),col = line_cols,lwd = 3)
#' @export RADnormalization_matrix
#' @seealso \code{\link{RADnormalization}} for normalize an abundance vector. This function return more details compared to \code{\link{RADnormalization_matrix}},
#'          \code{\link{representative_point}} for study the representative of groups of samples in a multi-dimensional scaling plot,
#'          \code{\link{representative_RAD}} for study the representative of group of norm rads.
#' @import scales sfsmisc

RADnormalization_matrix <- function(input,max_rank,average_over = 1,min_rank = 1,labels = FALSE,count_data = TRUE,sample_in_row = TRUE, method = "upperlimit",verbose = T) {
    if(!sample_in_row) input <- t(input)
    if(is.null(dim(input))) return(RADnormalization(input = input,max_rank = max_rank,average_over = average_over,min_rank = min_rank,labels = labels,count_data = count_data,method = method))

    norm <- c()
    norm_matrix <- c()
    for (i in 1:dim(input)[1]) {
        nrad <- RADnormalization(input = input[i,],max_rank = max_rank,average_over = average_over,min_rank = min_rank,labels = labels,count_data = count_data,method = method)
        norm <- c(norm,nrad)
        norm_matrix <- rbind(norm_matrix,nrad$norm_rad)
    if(verbose)  cat(i,"(",round(100*i/nrow(input),digits = 2),"%)","|")
    }
    inputs <- list(min_rank,max_rank,average_over,labels,sample_in_row,method,verbose)
    names(inputs) <- c("min_rank","max_rank","average_over","labels","sample_in_row","method","verbose")
    out <- list(norm_matrix,inputs)
    names(out) <- c("norm_matrix","inputs")
    return(out)
}
