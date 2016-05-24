#' Extract Drop Out from a Data.Frame
#' 
#' Find drop in Data.frame that contains multiple 
#' questions that had been asked sequentially.
#' 
#' @param df a data.frame
#' @param q_pos columns that contain questions
#' @examples
#' data(dropout)
#' dropout$pos <- extract_drop_out_from_df(dropout,2:10)
#' dropout$pos
#' @export
extract_drop_out_from_df <- function(df,q_pos){
  nms <- names(df[,q_pos])
  tf <- is.na(df[,q_pos])
  tpos <- apply(tf,1,which)
  tpos[lapply(tpos,length) == 0] <- NA
  # out vector contains drop out position
  # dpos <- 
  sapply(tpos,find_drop_out,clnms = nms)
  
}

