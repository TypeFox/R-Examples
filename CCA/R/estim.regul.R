"estim.regul" <-
function (X, Y, grid1 = NULL, grid2 = NULL, plt = TRUE) 
{
    if (is.null(grid1)) {grid1 = seq(0.001, 1, length = 5)} 
    if (is.null(grid2)) {grid2 = seq(0.001, 1, length = 5)}
    grid = expand.grid(grid1, grid2)
    res = apply(grid, 1, function(x) {loo(X, Y, x[1], x[2])})
    res.grid = cbind(grid, res)
    mat = matrix(res, nc = length(grid2))

    if (isTRUE(plt)) {
        img.estim.regul(list(grid1 = grid1, grid2 = grid2, mat = mat))
    }

    opt = res.grid[res.grid[, 3] == max(res.grid[, 3]), ]
    cat("  lambda1 = ", opt[[1]], "\n", " lambda2 = ", opt[[2]], "\n",
    "CV-score = ", opt[[3]], "\n")
    return(invisible(list(lambda1 = opt[[1]], lambda2 = opt[[2]], 
    CVscore = opt[[3]], grid1 = grid1, grid2 = grid2, mat = mat)))
}

