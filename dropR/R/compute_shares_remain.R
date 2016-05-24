#' Compute the Share of Remaining Participants
#' 
#' Compute the share of remaining participant per question.
#' 
#' @param df a data.frame
#' @param drop_out_pos drop out position
#' @param number_of_questions integer number of questions
#' @export
compute_shares_remain <- function(df,drop_out_pos,
                                  number_of_questions){
  grd <- expand.grid(1:number_of_questions,0)
  # count drop out positions, like table() 
  pos_count <- plyr::count(drop_out_pos)
  names(grd) <- names(pos_count)
  grd[pos_count$x,'freq'] <- pos_count$freq
  shares <- 1-(cumsum(grd$freq)/nrow(df))
  shares[-number_of_questions]
}

