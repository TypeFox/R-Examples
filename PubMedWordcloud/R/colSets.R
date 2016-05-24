#' @title plot colors
#' 
#' @param type palette names from the lists: Accent, Dark2, Pastel1, Pastel2, Paired, Set1, Set2, Set3
#' @export
#' @examples
#' # colors= colSets(type="Accent")
#' # colors= colSets(type="Paired")
#' # colors= colSets(type="Set3")
colSets <- function(type){
  n=switch(type,
          Accent = 8,
          Dark2 = 12,
          Pastel1 = 9,
          Pastel2 = 8,
          Paired = 12,
          Set1 = 9,
          Set2 = 8,
          Set3 = 12
      )
  brewer.pal(n,type)
}
