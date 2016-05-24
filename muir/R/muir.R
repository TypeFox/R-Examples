#' Explore Datasets with Trees
#'
#' This function allows users to easily and dynamically explore or document a
#' data.frame using a tree data structure. Columns of interest in the data.frame can
#' be provided to the function, as well as critieria for how they should be represented
#' in discrete nodes, to generate a data tree representing those columns and filters.
#' @param data A data.frame to be explored using trees
#' @param node.levels A character vector of columns from \code{data} that will be used to
#' construct the tree that are provided in the order that they should appear in the tree levels.
#'
#'  For each column, the user can add a suffix to the columnn name to indicate whether to generate
#'  nodes for all distinct values of the column in the date.frame, a specific number of values
#'  (i.e., the "Top (n)" values), and whether or not to aggregate remaining values into a separate
#'  "Other" node, or to use user-provided filter criteria for the column as provided in
#'  the \code{level.criteria} parameter. This does mean that the column names cannot have a ":"
#'  and must be replaced in the data.frame before being passed in to \code{muir} as
#'  the \code{data} param.
#'
#'  Values can be provided as "colname", "colname:*", "colname:3", "colname:+",
#'  or "colname:*+". The separator character ":" and the special characters in the suffix that
#'  follow (as outlined below) indicate which approach to take for each column.
#'  \itemize{
#'  \item Providing just the column name itself (e.g, "hp") will return results
#'  based on the operators and values provided in the \code{level.criteria} parameter
#'  for that column name. See \code{level.criteria} for more details.
#'  \item Providing the column name with an ":*"  suffix (e.g., "hp:*") will return a node for
#'  all distinct values for that column up to the limit imposed by the \code{node.limit} value.
#'  If the number of distinct values is greater than the \code{node.limit}, only the top "n"
#'  values (based on number of occurences) will be returned.
#'  \item Providing the column name with an ":\code{n}" suffix (e.g., "hp:3"), where
#'   \code{n} = a positive integer, will return a node for all distinct values for
#'   that column up to the limit imposed by the integer provided in \code{n}.
#'   If the number of distinct values is greater than the value provided in \code{n},
#'   only the top "n" values (based on number of occurences) will be returned.
#'  \item Providing the column name ending with an ":+" suffix (e.g., "hp:+") will return all the
#'  values provided in the \code{level.criteria} parameter for that column plus an extra node
#'  titled "Other" for that column that aggregates all the remaining values not included
#'  in the filter criteria provided in \code{level.criteria} for that column.
#'  \item Providing a column name ending with both symbols (e.g., "hp:*+", "hp:3+") in the suffix
#'   will return a node for all distinct values for that column up to the limit imposed by either
#'   the \code{node.limit} or the \code{n} value plus an additional "Other" node aggregating
#'   any remaining values beyond the \code{node.limit} or \code{n}, if applicable.
#'   If the number of distinct values is <= the \code{node.limit} or \code{n} then the "Other"
#'  node will not be created.
#'  }
#' @param node.limit Numeric value. When providing a column in \code{node.levels} with an ":*" suffix,
#' the \code{node.limit} will limit how many distinct values to actually process to prevent
#' run-away queries and unreadable trees. The limit defaults to 3 (not including an additional
#' 4th if requesting to provide an "Other" node as well with a ":*+" suffix). If the
#' number of distinct values for the column is greater than the \code{node.limit}, the tree
#' will include the Top "X" values based on count, where "X" = \code{node.limit}. If the
#' \code{node.limit} is greater than the number of distinct values for the column, it will
#' be ignored.
#' @param level.criteria A data.frame consisting of 4 character columns containing
#' column names (matching -- without suffixes -- the columns in \code{node.levels} that will
#' use the criteria in \code{level.criteria} to determine the filters used for each node),
#' an operator or boolean function (e.g., "==",">", "is.na", "is.null"), a value,
#' and a corresponding node title for the node displaying that criteria.
#'
#' E.g.,"wt, ">=", "4000", "Heavy Cars"
#'
#' @param label.vals Character vector of additional values to include in the node provided as a
#' character vector. The values must take the form of dplyr \code{summarise} functions
#' (as characters) and include the columns the functions should be run against (e.g.,
#' "min(hp)", "mean(hp)", etc.). If no custom suffix is added, the summary function itself
#' will be used as the label. Similar to \code{node.levels} a custom suffix can be added
#' using ":" to print a more meaningful label (e.g., "mean(hp):Avg HP"). In this example,
#' the label printed in the node will be "Avg HP:", otherwise it would be mean_hp (note
#' that the parens "(" and ")" are removed to be rendered in HTML without error). As with
#' \code{node.levels}, the column name itself cannot have a ":" and must be replaced in
#' the data.frame before being passed in to \code{muir} as the \code{data} param.
#' @param tree.dir Character. The direction the tree graph should be rendered. Defaults to "LR"
#' \enumerate{
#' \item Use "LR" for left-to-right
#' \item Use "RL" for right-to left
#' \item Use "TB" for top-to-bottom
#' \item User "BT" for bottom-to-top
#' }
#' @param show.percent Logical. Should nodes show the percent of records represented by
#' that node compared to the total number of records in \code{data.} Defaults to TRUE
#' @param num.precision Number of digits to print numeric label values out to
#' @param show.empty.child Logical. Show a balanced tree with children nodes that are all
#' empty or stop expanding the tree once there is a parent node that is empty.
#' Defaults to FALSE -- don't show empty children nodes
#' @param tree.height Numeric. Control tree height to zoom in/out on nodes. Passed to DiagrammeR
#' as \code{height} param. Defaults to -1, which appears to optimize the tree size
#' for viewing (still researching why exactly that works! :-))
#' @param tree.width Numberic. Control tree width to zoom in/out on nodes. Passed to DiagrammeR
#' as \code{width} param. Defaults to -1, which appears to best optimize the tree size
#' for viewing (still researching why exactly that works! :-))
#' @return An object of class \code{htmlwidget} (via DiagrammeR) that will
#' intelligently print itself into HTML in a variety of contexts
#' including the R console, within R Markdown documents,
#' and within Shiny output bindings.
#' @examples
#' \dontrun{
#' # Load in the 'mtcars' dataset
#' data(mtcars)
#'
#' # Basic exploration - show all values
#' mtTree <- muir(data = mtcars, node.levels = c("cyl:*", "carb:*"))
#' mtTree
#'
#' # Basic exploration - show all values overriding default node.limit
#' mtTree <- muir(data = mtcars, node.levels = c("cyl:*", "carb:*"), node.limit = 5)
#' mtTree
#'
#' # Show all values overriding default node.limit differently for each column
#' mtTree <- muir(data = mtcars, node.levels = c("cyl:2", "carb:5"))
#' mtTree
#'
#' # Show all values overriding default node.limit for each column
#' # and aggregating all distinct values above the node.limit into a
#' # separate "Other" column to collect remaining values
#'
#' # Top 2 occurring 'carb' values will be returned in their own nodes,
#' # remaining values/counts will be aggregated into a separate "Other" node
#' mtTree <- muir(data = mtcars, node.levels = c("cyl:2", "carb:2+"))
#' mtTree
#'
#' # Add additional calculations to each node output (dplyr::summarise functions)
#' mtTree <- muir(data = mtcars, node.levels = c("cyl:2", "carb:2+"),
#' label.vals = c("min(wt)", "max(wt)"))
#' mtTree
#'
#' # Make new label values more reader-friendly
#' mtTree <- muir(data = mtcars, node.levels = c("cyl:2", "carb:2+"),
#' label.vals = c("min(wt):Min Weight", "max(wt):Max Weight"))
#' mtTree
#'
#' # Instead of just returning top counts for columns provided in \code{node.levels},
#' # provide custom filter criteria and custom node titles in \code{label.vals}
#' # (criteria could also be read in from a csv file as a data.frame)
#' criteria <- data.frame(col = c("cyl", "cyl", "carb"),
#' oper = c("<", ">=", "=="),
#' val = c(4, 4, 2),
#' title = c("Less Than 4 Cylinders", "4 or More Cylinders", "2 Carburetors"))
#'
#' mtTree <- muir(data = mtcars, node.levels = c("cyl", "carb"),
#' level.criteria = criteria,
#' label.vals = c("min(wt):Min Weight", "max(wt):Max Weight"))
#' mtTree
#'
#' # Use same criteria but show all other values for the column where NOT
#' # EQUAL to the combination of the filters provided for that column (e.g., for cyl
#' # where !(cyl < 4 | cyl >= 4) in an "Other" node
#' mtTree <- muir(data = mtcars, node.levels = c("cyl:+", "carb:+"),
#' level.criteria = criteria,
#' label.vals = c("min(wt):Min Weight", "max(wt):Max Weight"))
#' mtTree
#'
#' # Show empty child nodes (balanced tree)
#' mtTree <- muir(data = mtcars, node.levels = c("cyl:+", "carb:+"),
#' level.criteria = criteria,
#' label.vals = c("min(wt):Min Weight", "max(wt):Max Weight"),
#' show.empty.child = TRUE)
#' mtTree
#'
#' # Save tree to HTML file with \code{htmlwidgets} package to working directory
#' mtTree <- muir(data = mtcars, node.levels = c("cyl:2", "carb:2+"))
#' htmlwidgets::saveWidget(mtTree, "mtTree.html")
#' }
#'
#' @import dplyr
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_trim
#' @importFrom stringr str_split
#' @importFrom stringr str_extract
#' @export
#'
#' @rdname muir

muir <- function(data, node.levels, node.limit = 3, level.criteria = NULL, label.vals = NULL,
                 tree.dir = "LR", show.percent = TRUE, num.precision = 2,
                 show.empty.child = FALSE, tree.height = -1, tree.width = -1) {

  # Hack for CRAN R CMD check
  parent <- NULL; rm("parent")
  index <- NULL; rm("index")
  node <- NULL; rm("node")
  cnt <- NULL; rm("cnt")
  colindex <- NULL; rm("colindex")

  # validate function parameters
  if(inherits(data, "data.table")) {
    data <- as.data.frame(data) #data.tables are not currently supported
  }

  if(!inherits(data,"data.frame")) {
    stop("data param must be a data.frame")
  }

  if(sum(grepl(":", colnames(data))) > 0) {
    # Using colnames that have an ":" embedded in the name will cause issues with muir processing
    # the custom node.level and label.vals suffixes.
    warning("There are columns names in the 'data' data.frame that include a ':'. ",
            "These columns should not be used in 'node.levels' or 'label.vals' or there may be errors.")
  }

  if(class(node.levels) != "character") {
    stop("node.levels must be a vector of type character. Provide the column names surrounded by quotes.")
  }

  if(!(tree.dir %in% c("LR", "RL", "TB", "BT"))) {
    stop("tree.dir must be either 'LR', 'RL', 'TB', or 'BT'.")
  }

  ## Remove factors from data so there are not filter or summarise errors later
  ## And makse sure all NULL and "" values are coerced to NA so they can be handled consitently
  i <- sapply(data, is.factor)
  data[i] <- lapply(data[i], as.character)

  ## Parse node.level columns to separate colnames from the requested criteria
  node.criteria <- data.frame(index = 1:length(node.levels),
                              do.call(rbind, stringr::str_split(node.levels, ":")),
                              stringsAsFactors = FALSE)

  # check whether any ":" suffixes were provided (and split off) for any node.criteria cols
  if (length(colnames(node.criteria)) == 3) {

    colnames(node.criteria) <- c("index", "col", "criteria")

    #make sure cols with no node.level criteria are marked as NA (they will be repeated in df)
    node.criteria[,3][node.criteria[,2] == node.criteria[,3]] <- NA

  } else { # no ":" found so no split occurred. Force 3rd column as NA

    node.criteria <- dplyr::mutate(node.criteria, criteria = NA)

  }

  colnames(node.criteria) <- c("index", "col", "criteria")
  node.levels <- as.vector(unlist(node.criteria$col))

  ## Ensure all node.levels provided exist as columns in the data df
  if(sum(unique(node.levels) %in% colnames(data)) != length(unique(node.levels))) {
    stop("All of the columns provided in node.levels do not appear in the data data.frame.")
  }

  ## Check if level.criteria will be used based on node.levels and if so that
  ## the value provided for level.criteria is in the correct format ad supports
  ## the columns provided in node.levels

  ## Get columns from node.levels that did not request all or 'n' number of values
  crit.cols <- unique(node.criteria$col[is.na(node.criteria$criteria)|
                                          node.criteria$criteria == "+"])

  ## Check level.criteria and stop if invalid
  if (length(crit.cols) > 0) {

    # validate the level.critera param
    validate.level.criteria(level.criteria) # will stop if error is found

    ## if the level.criteria param passed all the checks above, continue processing
    colnames(level.criteria) <- c("col", "oper", "val", "title") # set names
    level.criteria[] <- lapply(level.criteria, as.character) # and remove factors

    if(sum(unique(crit.cols) %in% unique(level.criteria$col)) !=
         length(unique(crit.cols))) {
      # all column names without ":*" or a ":{digit}" are not present in the level.criteria df
      stop("The level.criteria data.frame is missing values for the following columns: ",
           paste(setdiff(unique(crit.cols), unique(level.criteria$col)), collapse = ", "))
    }
  }

  ## Get data values for columns where the user wants nodes for each value and update criteria
  ## to use those values instead of any values passed in via level.criteria. Constrain number
  ## of values to the node.limit in descending order by total count of occurrences in the data df
  for(i in 1:nrow(node.criteria)) {
    if (grepl("\\*|\\d", node.criteria$criteria[i])) {

      ## Coerce NULLs and "" into NAs for columns being used without level.criteria
      data[,node.criteria$col[i]] <- sapply(data[,node.criteria$col[i]],
                                            function(x) ifelse(x == "", NA, x))

      # Getncounts for each distinct value for col[i] and arrange most to least
      col.values <- s_group_by(data, node.criteria$col[i])
      col.values <- dplyr::summarize(col.values, cnt = n())
      col.values <- dplyr::arrange(col.values, desc(cnt))

      if (grepl("\\*", node.criteria$criteria[i])) {

        # validate node.limit is valid
        if(!(is.numeric(node.limit)) | node.limit == 0) {
          stop("node.limit must be a positive integer.")
        }
        # get the top n values based on the node.limit value
        col.values <- dplyr::slice(col.values, 1:node.limit)

      } else {

        # get the top n values based on the param passed with the invidivual column
        col.values <- dplyr::slice(col.values, 1:stringr::str_extract(node.criteria$criteria[i], "\\d+"))
      }

      col.values <- unlist(s_select(col.values, node.criteria$col[i]))

      ##### TBD - come back to this and make the NAs work as operators
      new.criteria <- data.frame(col = as.character(node.criteria$col[i]),
                                 oper = "==",
                                 val = as.character(col.values),
                                 title = paste0(as.character(node.criteria$col[i]), " = ",
                                                as.character(col.values)),
                                 stringsAsFactors = FALSE)

      # If crtieria is NA, change operator from "==" to "is.na"
      new.criteria$oper[is.na(new.criteria$val)] <- "is.na"

      ## update/instantiate level.criteria with the node cols and values based on the top n values
      if(!is.null(level.criteria)) {
        # only filter if the current level.criteria isn't NULL and passes validation
        validate.level.criteria(level.criteria) # will stop if error is found
        level.criteria <- dplyr::filter(level.criteria, col != node.criteria$col[i])
      }

      level.criteria <- dplyr::bind_rows(level.criteria, new.criteria)
    }
  }

  # how many nodes will be written for each tree level
  colcounts <- dplyr::group_by(level.criteria, col)
  colcounts <- dplyr::summarise(colcounts, cnt = n())
  colcounts <- dplyr::inner_join(colcounts, node.criteria, by = "col")
  colcounts <- dplyr::arrange(colcounts, index)

  # number of df columns used based on node.levels
  numcols <- length(node.levels)

  #Add 1 to the count to account for the additonal "other" nodes where requested
  for(i in 1:nrow(colcounts)) {
    if (grepl("\\+", colcounts$criteria[i])) {
      colcounts$cnt[i] <- colcounts$cnt[i] + 1
    }
  }

  numnodes <- 1 + max(cumsum(cumprod(colcounts$cnt))) # plus 1 for level 0 node

  # establish initital column names
  nodedf_cols <- c("node", "colindex", "parent", "filter", "leaf_filter", "title", "nl_n")

  # If label vals are provided, add them as columns in the node df
  add.labels = NULL
  if (!is.null(label.vals)) {
    # check if custom label was provided, if not create label from dplyr function itself
    # Parse label.vals to separate out custom labels if provided
    label.vals <- data.frame(index = 1:length(label.vals),
                             do.call(rbind, stringr::str_split(label.vals, ":")),
                             stringsAsFactors = FALSE)

    # check whether any ":" suffixes were provided (and split off) for any label.vals
    if (length(colnames(label.vals)) == 2) {
      # no ":" found so no split occurred. Force 3rd column as NA
      label.vals <- dplyr::mutate(label.vals, label = NA)
    }

    colnames(label.vals) <- c("index", "fun", "label")

    # make sure label.vals with no customer label are marked as NA (they will be repeated in df)
    label.vals$label[label.vals$fun == label.vals$label] <- NA

    # remove special chars that may not render in HTML correctly
    label.vals$label <- stringr::str_replace_all(label.vals$label, "[[:punct:]]", " ")
    label.vals$label <- stringr::str_trim(label.vals$label)

    for (l in 1:nrow(label.vals)) {

      if(is.na(label.vals$label[l])) {

        init <- label.vals$fun[l]
        x1 <- str_replace(init, "\\(", "_")
        x2 <- str_replace(x1, "\\)","")
        add.labels <- c(add.labels, paste0("nl_", x2))

      } else {

        add.labels <- c(add.labels, paste0("nl_", label.vals$label[l]))
      }

    }
    # add columns matching the names that should be printed for each additional label
    nodedf_cols <- c(nodedf_cols,add.labels)
  }

  # Instantiate df based on expected number of nodes
  nodedf <- data.frame(matrix(ncol = length(nodedf_cols), nrow = numnodes))
  colnames(nodedf) <- nodedf_cols

  ## Set parent node (level 0) values
  nodedf$node[1] <- 1
  nodedf$colindex[1] <- 0
  nodedf$parent[1] <- "None"
  nodedf$filter[1] <- "None"
  nodedf$leaf_filter[1] <- "None"
  nodedf$title[1] <- "All"
  nodedf$nl_n[1] <- as.integer(dplyr::summarise(data, n()))

  ## Add values for additional lables if provided
  if (!is.null(add.labels)) {
    for (l in 1:length(add.labels)) {
      if (nodedf$nl_n[1] > 0) {

        nodedf[, add.labels[l]][1] <- s_summarise(data, label.vals$fun[l])[[1]]

        ## if value is numeric/double, format to decimals places = num.precision
        ## character columns/values will be coerced to 'NA'
        nodedf[, add.labels[l]][1] <- format(round(as.numeric(nodedf[, add.labels[l]][1]),
                                                   digits = num.precision),
                                             nsmall = num.precision)
      }
    }
  }

  ### Build out tree table nodedf
  cur_node <- 2 # cur_node = 1 was the parent (level 0) node

  for (r in 1:numcols) { # loop through the columns provided from the dataset

    parnodes <- dplyr::filter(nodedf, colindex == (r - 1))
    leaves <- dplyr::filter(level.criteria, col == node.levels[r])

    for (pn in 1:nrow(parnodes)) { # Loop through the parent nodes at the current level

      root <- parnodes$node[pn]
      leaf_filter <- NULL

      if (!show.empty.child & parnodes$nl_n[pn] == 0) {

        ## if parent node has no values, set the child leaves as NA/0
        nodedf$node[cur_node] <- cur_node
        nodedf$colindex[cur_node] <- r
        nodedf$parent[cur_node] <- root
        nodedf$filter[cur_node] <- NA
        nodedf$leaf_filter[cur_node] <- NA
        nodedf$title[cur_node] <- NA
        nodedf$nl_n[cur_node] <- 0

        if (!is.null(add.labels)) {
          for (l in 1:length(add.labels)) {

            nodedf[, add.labels[l]][cur_node] <- NA
          }
        }

        cur_node <- cur_node + 1

      } else {

        if (parnodes$filter[pn] == "None") {
          pfltr <- NULL
        } else {
          pfltr <- parnodes$filter[pn]
        }

        for (n in 1:nrow(leaves)) { # Add leaves for each parent node at this level


          nodedf$node[cur_node] <- cur_node
          nodedf$colindex[cur_node] <- r
          nodedf$title[cur_node] <- leaves[[4]][n]
          nodedf$parent[cur_node] <- root
          if (leaves[[3]][n] == "" | is.na(leaves[[3]][n])) {
            nodedf$leaf_filter[cur_node] <- ff(leaves[[2]][n], leaves[[1]][n])
          } else {
            nodedf$leaf_filter[cur_node] <- paste0(leaves[[1]][n], leaves[[2]][n], "\"", leaves[[3]][n], "\"")
          }

          if (!is.null(pfltr)) {
            cur_filter <- paste0(pfltr, ",", nodedf$leaf_filter[cur_node])
          } else {
            cur_filter <- nodedf$leaf_filter[cur_node]
          }

          nodedf$filter[cur_node] <- cur_filter
          cur_filter_df <- s_filter(data, cur_filter)
          nodedf$nl_n[cur_node] <- as.integer(dplyr::summarise(cur_filter_df, n()))

          ## Add values for additional lables if provided
          if (!is.null(add.labels)) { #
            for (l in 1:length(add.labels)) {
              if (nodedf$nl_n[cur_node] > 0) {

                nodedf[, add.labels[l]][cur_node] <- s_summarise(cur_filter_df, label.vals$fun[l])[[1]]

                ## if value is numeric/double, format to decimals places = num.precision
                ## character columns/values will be coerced to 'NA'
                nodedf[, add.labels[l]][cur_node] <- format(round(as.numeric(nodedf[, add.labels[l]][cur_node]),
                                                                  digits = num.precision),
                                                            nsmall = num.precision)
              }
            }
          }

          cur_node <- cur_node + 1

        }

        ## Add "Other" Node to capture remaining data not included in the child filters
        if(grepl("\\+",node.criteria$criteria[node.criteria$col == node.levels[r]])) {
          nodedf$node[cur_node] <- cur_node
          nodedf$colindex[cur_node] <- r
          nodedf$title[cur_node] <- "Other"

          leafdf<- filter(nodedf, parent == root)
          leafdf <- select(leafdf, leaf_filter)
          prev_leaf_filters <- leafdf
          nodedf$leaf_filter[cur_node] <- paste0("!", paste0(leafdf[[1]], collapse = ",!"))

          if (!is.null(pfltr)) {
            cur_filter <- paste0(pfltr, ",", nodedf$leaf_filter[cur_node])
          } else {
            cur_filter <- nodedf$leaf_filter[cur_node]
          }

          nodedf$filter[cur_node] <- cur_filter
          cur_filter_df <- s_filter(data, cur_filter)
          nodedf$nl_n[cur_node] <- as.integer(dplyr::summarise(cur_filter_df, n()))
          nodedf$parent[cur_node] <- root

          ## Add values for additional lables if provided
          if (!is.null(add.labels)) {
            for (l in 1:length(add.labels)) {
              if (nodedf$nl_n[cur_node] > 0) {

                nodedf[, add.labels[l]][cur_node] <- s_summarise(cur_filter_df, label.vals$fun[l])[[1]]

                ## if value is numeric/double, format to decimals places = num.precision
                ## character values will be coerced to 'NA'
                nodedf[, add.labels[l]][cur_node] <- format(round(as.numeric(nodedf[, add.labels[l]][cur_node]),
                                                                  digits = num.precision),
                                                            nsmall = num.precision)
              }
            }
          }
          cur_node <- cur_node + 1
        }

      }

    }
  }

  # add percent of total n() value to the nodedf for inclusion as a node label later
  if (show.percent) {
    nl_pct <- format(round((as.numeric(nodedf$nl_n)/as.numeric(nodedf$nl_n[1])) * 100,
                          digits = num.precision), nsmall = num.precision)
    nodedf <- dplyr::mutate(nodedf, "nl_%" = nl_pct)
  }

  # remove any left-over NA nodes due to parent nodes with no values and show.empty.child = TRUE
  # and format number column to include commas
  nodedf <- filter(nodedf, !is.na(node))
  nodedf$nl_n <- formatC(nodedf$nl_n, format = "d", big.mark = ",")

  # create htmlwidget/DiagrammeR object
  tree <- build_tree(nodedf, tree.dir, tree.height, tree.width)

  # return/render generated tree
  tree
}
