#' Preview the levels in a parsed Eprime file
#'
#' @details \code{preview_levels} prints out the unique combinations of
#' Eprime.Level number, Procedure, and Running in the frame list.
#' \code{preview_frames} prints out example frame from each of the unique
#' levels. \code{preview_eprime} does both.
#'
#' @param frame_list a FrameList (a list of EprimeFrames)
#' @return Nothing. Preview text is printed to the console.
#' @export
preview_eprime <- function(frame_list) {
  preview_levels(frame_list)
  preview_frames(frame_list)
  invisible(NULL)
}

#' @rdname preview_eprime
#' @export
preview_levels <- function(frame_list) {
  prep <- preview_prep(frame_list)
  cat("Level Counts: \n")
  print(prep$row_counts, row.names = FALSE)
  invisible(NULL)
}

#' @rdname preview_eprime
#' @export
preview_frames <- function(frame_list) {
  prep <- preview_prep(frame_list)

  for(chunk_num in seq_along(prep$unique_frames)) {
    curr_row <- prep$unique_rows[chunk_num, ]
    curr_chunk <- prep$unique_frames[[chunk_num]]

    cat("\n")
    print(curr_row, row.names = FALSE)
    str(curr_chunk)
  }

  invisible(NULL)
}

preview_prep <- function(frame_list) {
  keys <- c("Eprime.Level", "Running", "Procedure")
  main_cols <- pick_apply(keys, frame_list)
  full_table <- to_data_frame(main_cols)[keys]

  unique_rows <- unique(full_table)
  unique_frames <- as.FrameList(frame_list[as.numeric(row.names(unique_rows))])

  # Include frequency count
  row_counts <- plyr::count(full_table)
  row_counts <- plyr::join(unique_rows, row_counts, by = keys)

  list(row_counts = row_counts, unique_rows = unique_rows,
       unique_frames = unique_frames)
}




