tableStyles <- function(x, useRowNames = TRUE, header = NULL,
                        cgroup = NULL, rgroup = NULL)
{
   styles <- getStyles()
   tableDim <- dim(x)
   if(useRowNames) tableDim[2] <- tableDim[2] + 1

   has <- function(x) !is.null(x) && x != ""

   textForm <- if(has(styles$cellText)) matrix(rep(styles$cellText, tableDim[1] * tableDim[2]), nrow = tableDim[1]) else NULL
   cellForm <- if(has(styles$cell)) matrix(rep(styles$cell, tableDim[1] * tableDim[2]), nrow = tableDim[1]) else NULL

   if(!is.null(header))
   {
      headerText <- if(has(styles$headerText)) matrix(rep(styles$headerText, tableDim[2]), nrow = 1) else NULL
      headerCell <- if(has(styles$header)) matrix(rep(styles$header, tableDim[2]), nrow = 1) else NULL
   } else {
      headerText <- NULL
      headerCell <- NULL
   }

   if(!is.null(cgroup) && is.data.frame(cgroup))
   {
      cgroupText <- if(has(styles$headerText)) matrix(rep(styles$headerText, nrow(cgroup)), nrow = 1) else NULL
      cgroupCell <- if(has(styles$header)) matrix(rep(styles$header, nrow(cgroup)), nrow = 1) else NULL
   } else {
      cgroupText <- NULL
      cgroupCell <- NULL
   }

   if(!is.null(rgroup) && is.data.frame(rgroup))
   {
      rgroupText <- if(has(styles$cellText)) matrix(rep(styles$cellText, nrow(rgroup)), nrow = 1) else NULL
      rgroupCell <- if(has(styles$cell)) matrix(rep(styles$cell, nrow(rgroup)), nrow = 1) else NULL
   } else {
      rgroupText <- NULL
      rgroupCell <- NULL
   }

   list(table = styles$table, text = textForm, cell = cellForm, header = headerText, headerCell = headerCell,
        cgroupText = cgroupText, cgroupCell = cgroupCell, rgroupText = rgroupText, rgroupCell = rgroupCell)
}
