"IndexToSet" <- 
function(set.index) {
   E <- if (any(set.index == 0)) which(set.index == 0)
        else NULL 
   L <- if (any(set.index == -1)) which(set.index == -1)
        else NULL
   R <- if (any(set.index == 1)) which(set.index == 1)
        else NULL
   obj <- list(Elbow = E, Left = L, Right = R)
obj
}
