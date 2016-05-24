ConfigShowSimilarEntries <- function(configuration, dataset_name = NULL, var_name = NULL, main_path = NULL, file_path = NULL, nc_var_name = NULL, suffix = NULL, varmin = NULL, varmax = NULL, n_results = 10) {
  ## Simon White: http://www.catalysoft.com/articles/StrikeAMatch.html
  getBigrams <- function(str) {
    bigramlst <- list()
    for (i in 1:(nchar(str) - 1)) {
      bigramlst[[i]] <- substr(str, i, i + 1)
    }
    return(bigramlst)
  }

  strSimilarity <- function(str1, str2) {
    str1 <- tolower(str1)
    str2 <- tolower(str2)    

    if (is.null(str1)) {
      str1 <- ""
    } else if (is.na(str1)) {
      str1 <- ""
    } 
    if (is.null(str2)) {
      str2 <- ""
    } else if (is.na(str2)) {
      str2 <- ""
    }
    if (nchar(str1) <= 1 && nchar(str2) <= 1) {
      return (ifelse(str1 == str2, 1, 0))
    } else if (nchar(str1) == 1) {
      return (ifelse(grepl(str1, str2, fixed = TRUE), 1, 0))
    } else if (nchar(str2) == 1) {
      return (ifelse(grepl(str2, str1, fixed = TRUE), 1, 0))
    } else if (nchar(str1) == 0 || nchar(str2) == 0) {
      return (0)
    } else {
      pairs1 <- getBigrams(str1)
      pairs2 <- getBigrams(str2)
      unionlen <- length(pairs1) + length(pairs2)
      hit_count <- 0
      for (x in 1:length(pairs1)) {
        for(y in 1:length(pairs2)) {
          if (pairs1[[x]] == pairs2[[y]]) {
            hit_count <- hit_count + 1
          }
        }
      }
      return ((2.0 * hit_count) / unionlen)
    }
  }
  
  strSimilarityVec <- Vectorize(strSimilarity, c('str1', 'str2'), USE.NAMES = FALSE)

  all_tables <- c('experiments', 'observations')
  all_fields <- c('dataset_name', 'var_name', 'main_path', 'file_path', 'nc_var_name', 'suffix', 'varmin', 'varmax')
  selected_fields <- which(unlist(lapply(as.list(match.call())[all_fields], function (x) !is.null(x))))
  values <- unlist(as.list(match.call())[all_fields[selected_fields]], use.names = FALSE)

  if (length(selected_fields) < 1) {
    stop("There must be at least one selected field ('dataset_name', 'var_name', 'main_path', 'file_path', 'nc_var_name', 'suffix', 'varmin' or 'varmax').")
  }

  similarities <- list()
  for (table in all_tables) {
    similarities[[table]] <- vector("list", 4)
    for (level in 1:4) {
      if (length(configuration[[table]][[level]]) > 0) {
        similarities[[table]][[level]] <- unlist(lapply(configuration[[table]][[level]], function(x) mean(strSimilarityVec(x[selected_fields], values))))
      }
    }
  }
 
  n_results <- min(n_results, length(unlist(similarities)))
  threshold <- sort(unlist(similarities, use.names = FALSE), decreasing = TRUE)[n_results]
  n_minimums <- sum(sort(unlist(similarities, use.names = FALSE), decreasing = TRUE)[1:n_results] == threshold)
  
  matches <- list()
  n_picked_minimums <- 0
  for (table in all_tables) {
    matches[[table]] <- list()
    line_numbers <- c()
    offset <- 0
    for (level in 1:4) {
      matches_to_pick <- which(similarities[[table]][[level]] > threshold)
      if (n_picked_minimums < n_minimums) {
        minimums <- which(similarities[[table]][[level]] == threshold)
        if (length(minimums) + n_picked_minimums > n_minimums) {
          minimums <- minimums[1:(n_minimums - n_picked_minimums)]
        }
        matches_to_pick <- c(matches_to_pick, minimums)
        n_picked_minimums <- n_picked_minimums + length(minimums)
      }
      line_numbers <- c(line_numbers, matches_to_pick + offset)
      offset <- offset + length(similarities[[table]][[level]])
      matches[[table]][[level]] <- configuration[[table]][[level]][matches_to_pick]
    }
    dataset_type <- ifelse(grepl('experiments', table), 'experiments', 'observations')
    ConfigShowTable(matches, dataset_type, line_numbers)
  }
  
  invisible(matches)
}
