#' Identify site specific kinase substrate relationships using dynamic data.
#'
#' @description Using this package you can combine known site specific 
#' kinase substrate relationships with dynamic experimental data and determine active 
#' kinases and their substrates.
#' 
#' @author Westa Domanova
#' @docType package
#' @name ksrlive
NULL

#' Create a kinase substrate relationship list from a data frame
#'
#' \code{KSR.list} returns a list of kinase substrate relationships
#'
#' The function KSR.list creates a list of kinase substrate relationships from
#' a data frame and can combine kinase families into one list. Substrates occuring 
#' in multiple lists can be excluded. 
#'
#' @param df data frame of kinase substrate relationships with substrate 
#' identifier in the first column and kinase identifier in the second column.
#' @param kinasefamilies named list of kinase identifiers that have to be combined, 
#' one list per kinase family, list will be named after first family member
#' @param exclusive logical, if TRUE only substrates exclusive to the kinase
#' will be included in the list (substrates with multiple kinases will be excluded)
#' @return named list of substrate identifiers, with the corresponding kinase
#' identifiers as the list names
#'
#' @examples
#' data(phosphonetworkdf)
#' data(datakin)
#' 
#' # first column has to be substrate id, second kinase id
#' kin_data <- KSR.list(phosphonetwork_df[, c("SUB_IDENT", "KIN_ACC_ID")]) 
#' # Akt1 and Akt2 belong to the same kinase family, combine their substrates 
#' # into one list and name the list after the first family member
#' fam <- list(akt = c("P31749", "P31751"))
#' kin_data_fam <- KSR.list(phosphonetwork_df[, c("SUB_IDENT", "KIN_ACC_ID")], 
#' kinasefamilies = fam)
#' 
#' # only include phosphosites appearing once
#' kin_data_fam_exc <- KSR.list(phosphonetwork_df[, c("SUB_IDENT", "KIN_ACC_ID")], 
#'                              kinasefamilies = fam,
#'                              exclusive = TRUE)
#'@export
#'@importFrom stats na.omit

KSR.list <- function(df, kinasefamilies = NULL, exclusive = FALSE){
  temp <- split(df[,1], f = df[,2])
  # out<-temp
  out_cl <- lapply(temp, unique) ## delete duplicates
  out_cl <- lapply(out_cl, function(x){as.character(na.omit(x))}) ## delete NAs
  # delete empty lists
  full <- which(sapply(out_cl, function(x){length(x) > 0}))
  if (length(full) > 0) {
    out_cl <- out_cl[full]
    # combine kinasefamilies together
    if (is.null(kinasefamilies)) {
      out_fam <- out_cl
    }else{
      out_fam <- lapply(kinasefamilies, 
                        function(x){unique(unlist(out_cl[unlist(x)]))})
      names(out_fam) <- sapply(kinasefamilies, "[[", 1)
      
      out_cl <- out_cl[-which(names(out_cl) %in% unlist(kinasefamilies))]
      out_fam <- append(out_cl, out_fam)
    }
    if (!exclusive) {
      out_final <- out_fam
    }else{
      fam_df <- data.frame(sub = unlist(out_fam), 
                           kin = names(out_fam)[rep(seq_along(out_fam), 
                                                    lapply(out_fam, length))],
                           stringsAsFactors = FALSE)
      substr <- split(fam_df[ , 2], f = fam_df[ , 1])
      substr_cl <- lapply(substr, unique) ## delete duplicates
      ### find exclusive substrates
      sub_kinases <- sapply(substr_cl, length)
      ## only one kinase
      onekin <- which(sub_kinases == 1)
      # twokin<-which(sub.kinases==2)
      fam_df<-fam_df[fam_df[ , 1] %in% names(substr_cl)[onekin], ]
      temp <- split(fam_df[ , 1], f = fam_df[ , 2])
      out_final <- lapply(temp, unique) ## delete duplicates
    }
  }else{
    print("No lists available")
    out_final <- NULL
  }
  return(out_final)
}

#' Create random data
#'
#' \code{random.data} returns a data frame of random numeric values 
#'
#' The function random.data returns a data frame of random numeric values with the same 
#' number of columns as the input data and with n-nrow(data) rows. By default the values are drawn
#' from a uniform distribution of values between the minimum and the maximum of the input data. Values
#' can be drawn from background data instead if included.
#'
#' @param data data frame of time course of substrates, each substrate is a row
#' @param back_data data frame of numeric values that can to be used as background data, 
#' if not provided a values are drawn from a uniform distribution between minimum and maximum
#' of input data
#' @param n numeric specifying how many rows should be contained in the resulting data frame
#' @param random.seed numeric used as seed
#' @return data frame of random numeric values with n-nrow(data) rows and same number of 
#' columns as input data
#'
#' @examples
#' data(phosphonetworkdf)
#' data(datakin)
#' # only need what is present in data
#' phosphonetwork_data <- phosphonetwork_df[
#' phosphonetwork_df[,"SUB_IDENT"] %in% data_kin[,"SUB_IDENT"]
#' ,]
#' fam <- list(akt = c("P31749", "P31751"))
#' kin_data_fam_exc <- KSR.list(phosphonetwork_data[, c("SUB_IDENT", "KIN_ACC_ID")], 
#'                              kinasefamilies = fam,
#'                              exclusive = TRUE)
#' # only do for Akt and Mtor (P31749, P42345)
#' substrate_profiles <- lapply(kin_data_fam_exc[c("P31749", "P42345")], 
#' function(x){data_kin[match(x, data_kin[,"SUB_IDENT"]),1:9]})
#'
#' substrate_profiles_random <- lapply(substrate_profiles, 
#' function(x){rbind(x, random.data(x, random.seed = 123))})
#' @export
#' @importFrom stats runif

### add random data
random.data <- function(data, back_data = NULL, n = 50, random.seed = NULL){
      if(n <= nrow(data)){
            n <- 0
            message("Data has more observations than given maximum.")
            return(NULL)
      }else{
            if(!is.null(back_data)){
                  indata <- match(rownames(data), rownames(back_data))
                  back_data <- back_data[-indata, ]
                  if(!is.null(random.seed)){
                        set.seed(random.seed)
                  }
                  random_data <- back_data[sample(rownames(back_data), n - nrow(data)), ]
                  rownames(random_data) <- paste(rownames(random_data), rep("_random", nrow(random_data)), sep = "")
            }else{
                  if(!is.null(random.seed)){
                        set.seed(random.seed)
                  }
                  random_data <- t(replicate(n-nrow(data), runif(ncol(data), min = min(data), max = max(data))))
                  rownames(random_data) <- paste(c(1:nrow(random_data)), rep("_random", nrow(random_data)), sep = "")
                  colnames(random_data) <- colnames(data)
            }
            return(random_data)
      }
}

#' Return clustering assignments produced by tight.clust
#'
#' \code{clustering} returns vectors of clustering assignments
#'
#' The function clustering creates a named list of cluster assignments for substrates. 
#'
#' @param tightclust list of objects returned by the tight.clust function
#' @param data data frame of time course of substrates, each substrate is a row
#' @return named list containing named vectors of cluster assignments, names correspond to rownames in data
#' and names of list are kinase identifiers
#'
#' @importFrom tightClust tight.clust
#' @examples
#' data(phosphonetworkdf)
#' data(datakin)
#' # only need what is present in data
#' phosphonetwork_data <- phosphonetwork_df[
#' phosphonetwork_df[,"SUB_IDENT"] %in% data_kin[,"SUB_IDENT"]
#' ,]
#' fam <- list(akt = c("P31749", "P31751"))
#' kin_data_fam_exc <- KSR.list(phosphonetwork_data[, c("SUB_IDENT", "KIN_ACC_ID")], 
#'                              kinasefamilies = fam,
#'                              exclusive = TRUE)
#' # only do for Akt and Mtor (P31749, P42345)
#' substrate_profiles <- lapply(kin_data_fam_exc[c("P31749", "P42345")], 
#' function(x){data_kin[match(x, data_kin[,"SUB_IDENT"]),1:9]})
#'
#' substrate_profiles_random <- lapply(substrate_profiles, 
#' function(x){rbind(x, random.data(x, random.seed = 123))})
#' 
#' target <- 3
#' substrate_profiles_tight <- lapply(substrate_profiles_random, function(x){
#' tightClust::tight.clust(x, target = target, k.min = 7, resamp.num = 100, random.seed = 12345)
#' })
#' 
#' kin_clust<- mapply(function(x,y){clustering(x, y)}, 
#'                         substrate_profiles_tight, substrate_profiles, SIMPLIFY = FALSE)
#'@export
clustering <- function(tightclust, data){
      clustnum <- unique(tightclust$cluster)
      clustnum <- clustnum[-which(clustnum == -1)]
      clust <- lapply(c(1:length(clustnum)), function(x){which(tightclust$cluster == x)})
      
      ### delete all clusters that are only random data
      cluster <- lapply(clust, function(x){
            as.character(na.omit(rownames(data)[x]))
      })
      big <- which(sapply(cluster, length) > 1)
      ### create cluster assignment vector
      clust_assg <- numeric(length = nrow(data))
      names(clust_assg) <- rownames(data)
      if(length(big) > 0){
            for(i in 1:length(big)){
                  clust_assg[cluster[[big[i]]]] <- i
            }
      }else{
            message("Does not have a cluster.")
      }
      return(clust_assg)
}

#' Find clusters containing core substrates
#'
#' \code{clust.expand} returns a list of kinase substrate relationships
#'
#' The function clust.expand takes the resulting core substrates from the exclusive clustering and finds
#' the corresponding substrate clusters in the clustering using all substrates. 
#'
#' @param clust named list containing named vectors of cluster assignments, names correspond to rownames in data
#' and names of list are kinase identifiers (result of clustering performed using exclusive substrates)
#' @param clust_all named list containing named vectors of cluster assignments, names correspond to rownames in data
#' and names of list are kinase identifiers (result of clustering performed using all substrates)
#' @param diff character vector of substrate identifiers that are differentially regulated 
#' @return named list containing named vectors of cluster assignments, names correspond to rownames in data
#' and names of list are kinase identifiers
#' @importFrom tightClust tight.clust
#' @examples
#' data(phosphonetworkdf)
#' data(datakin)
#' # only need what is present in data
#' phosphonetwork_data <- phosphonetwork_df[
#' phosphonetwork_df[,"SUB_IDENT"] %in% data_kin[,"SUB_IDENT"]
#' ,]
#' fam <- list(akt = c("P31749", "P31751"))
#' kin_data_fam_exc <- KSR.list(phosphonetwork_data[, c("SUB_IDENT", "KIN_ACC_ID")], 
#'                              kinasefamilies = fam,
#'                              exclusive = TRUE)
#'                              
#' # only do for Akt and Mtor (P31749, P42345)
#' substrate_profiles <- lapply(kin_data_fam_exc[c("P31749", "P42345")], 
#' function(x){data_kin[match(x, data_kin[,"SUB_IDENT"]),1:9]})
#'
#' substrate_profiles_random <- lapply(substrate_profiles, 
#' function(x){rbind(x, random.data(x, random.seed = 123))})
#' 
#' target <- 3
#' substrate_profiles_tight <- lapply(substrate_profiles_random, function(x){
#' tightClust::tight.clust(x, target = target, k.min = 7, resamp.num = 100, random.seed = 12345)
#' })
#' 
#' kin_clust<- mapply(function(x,y){clustering(x, y)}, 
#'                         substrate_profiles_tight, substrate_profiles, SIMPLIFY = FALSE)
#'                         
#' # do clustering using all available substrates
#' kin_data_fam <- KSR.list(phosphonetwork_data[, c("SUB_IDENT", "KIN_ACC_ID")], 
#'                          kinasefamilies = fam)
#' 
#' substrate_profiles_all <- lapply(kin_data_fam[c("P31749", "P42345")], 
#' function(x){data_kin[match(x, data_kin[,"SUB_IDENT"]),1:9]})
#'
#' substrate_profiles_random_all <- lapply(substrate_profiles_all, 
#'                        function(x){rbind(x, random.data(x, random.seed = 123))})
#' 
#' target <- 3
#' substrate_profiles_tight_all <- lapply(substrate_profiles_random_all, function(x){
#' tightClust::tight.clust(x, target = target, k.min = 7, resamp.num = 100, random.seed = 12345)
#' })
#' 
#' kin_clust_all <- mapply(function(x,y){clustering(x, y)}, 
#'                         substrate_profiles_tight_all, substrate_profiles_all, 
#'                         SIMPLIFY = FALSE)
#'                         
#' expand_all <- mapply(function(x,y){clust.expand(x, y)}, 
#'                         kin_clust, kin_clust_all, SIMPLIFY = FALSE)
#'@export
clust.expand <- function(clust, clust_all, diff = NULL){
      # find all possible cluster cores
      imp_clu <- names(table(clust))
      if(any(imp_clu == 0)){
            imp_clu <- imp_clu[-which(imp_clu == 0)]
      }
      clust_mem <- lapply(imp_clu, function(x){
            names(clust)[which(clust == as.numeric(x))]
      })
      # if any cluster member is not in a cluster when clustering all
      # then remove it from the member list
      if (any(clust_all[unlist(clust_mem)] == 0)) {
            ind <- lapply(clust_mem, function(x){
                  which(clust_all[x] == 0)
            })
            clust_mem <- lapply(c(1:length(ind)), function(x){
                  if (length(ind[[x]]) != 0){
                        clust_mem[[x]][-ind[[x]]]
                  }else{
                        clust_mem[[x]]
                  }
            })
      }
      
      # test whether the core is differentially regulated
      if (!is.null(diff)) {
            if (any(!(unlist(clust_mem) %in% diff))) {
                  ind <- lapply(clust_mem, function(x){
                        which(!(unlist(x) %in% diff))
                  })
                  clust_mem <- lapply(c(1:length(ind)), function(x){
                        if (length(ind[[x]]) != 0) {
                              clust_mem[[x]][-ind[[x]]]
                        }else{
                              clust_mem[[x]]
                        }
                  })
            }
      }
      # if no cluster found
      if (all(sapply(clust_mem, length) == 0)) {
            message("Does not have a cluster")
            expand_clust <- NA
      }else{
            expand <- lapply(clust_mem, function(x){
                  clust_all[unlist(x)]
            })
            expand_t <- lapply(expand, table)
            cluster <- lapply(expand_t, function(x){as.numeric(names(unlist(x)))})
            expand_clust <- clust_all
            expand_clust[! (expand_clust %in% unlist(cluster))] <- 0
            hascl <- which(sapply(cluster, length) > 0)
            clustnum_vec <- c(1:length(hascl))
            clust_ind <- list()
            for(i in 1:length(hascl)){
                  clust_ind[[i]] <- which(expand_clust %in% cluster[[hascl[i]]])
            }
            for(i in 1:length(clust_ind)){
                  expand_clust[clust_ind[[i]]] <- clustnum_vec[i]
            }
      }
return(expand_clust)
}
