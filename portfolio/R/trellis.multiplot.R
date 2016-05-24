################################################################################
##
## $Id: trellis.multiplot.R 346 2006-10-01 05:08:55Z enos $
##
## Plot some trellis objects on the same page.
##
################################################################################

.trellis.multiplot <- function(x){

  if(!class(x) == "list"){
    stop("x must be of class 'list'")
  }
  if(!all(sapply(x, class) == "trellis")){
    stop("All elements of x must be of class 'trellis'")
  }

  num <- length(x)
  n <- ceiling(sqrt(num))
  m <- n
  
  if(num == 2){
    m <- 1
  }
  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(n, m)))

  elt.num <- 1

  for(i in 1:n){
    for(j in 1:m){
      if(elt.num > num) break()
      row <- i
      col <- j
      if(num == 3 && row == 2) col <- c(1,2)
        
      pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
      
      print(x[[elt.num]], newpage = FALSE)
      popViewport(1)
      elt.num <- elt.num + 1
    }
  }
}
