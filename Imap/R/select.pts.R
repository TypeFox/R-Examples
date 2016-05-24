select.pts <-
function (pts, list.of.lists.obj = 1, outside.poly = FALSE, col = "blue", alpha = 0.5, lty = 1, ...) 
{

    if(is.list(pts) == TRUE & !is.null(names(pts[[list.of.lists.obj]]))) {
       if(names(pts[[list.of.lists.obj]])[1] == 'll') 
         pts <- pts[[list.of.lists.obj]]$ll
    }

    pts <- pts # Needed if e.g. imap() is called inside of select.pts: select.pts(imap())

    col <- col.alpha(col, alpha)

    Poly <- draw.polygon(col = col, lty = lty)
    tf <- inside.polygon(pts, Poly)
    points(pts[tf, 1], pts[tf, 2], ...)

    if(!any(tf, na.rm=TRUE))
      stop('No points selected. A smoother polygon may be needed.')

    if (outside.poly) 
        pts[!tf, ]
    else pts[tf, ]
} 
