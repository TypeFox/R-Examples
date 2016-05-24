"img.estim.regul" <-
function (estim) 
{
    grid1 = estim$grid1
    grid2 = estim$grid2
    mat = estim$mat
    nlevel = min(512, length(c(unique(mat))))

    par(mar = c(4.5, 4.5, 4.5, 6))
    image(grid1, grid2, mat, col = heat.colors(nlevel), 
    xlab = expression(lambda[1]), ylab=expression(lambda[2]), 
    main = expression(CV(lambda[1], lambda[2])), axes = FALSE,
    zlim = c(min(mat), max(mat)), oldstyle = TRUE)

    if (length(grid1) > 10) {
        grid1 = seq(min(grid1), max(grid1), length = 11)
    }
    if (length(grid2) > 10) {
        grid2 = seq(min(grid2), max(grid2), length = 11)
    }

    axis(1, at = grid1, labels = as.character(round(grid1, 4)))
    axis(2, at = grid2, labels = as.character(round(grid2, 4)))
    box()

    image.plot(legend.only = TRUE, nlevel = nlevel, 
            zlim = c(min(mat), max(mat)), col = heat.colors(nlevel),
            legend.width = 1.8)   
}

