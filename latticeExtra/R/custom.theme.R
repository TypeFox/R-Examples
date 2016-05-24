

## Construct a custom theme based on supplied colors.  Defaults to
## colors from RColorBrewer

custom.theme <-
    function(symbol = brewer.pal(n = 8, name = "Dark2"),
             fill = brewer.pal(n = 12, name = "Set3"),
             region = brewer.pal(n = 11, name = "Spectral"),
             reference = "#e8e8e8",
             bg = "transparent",
             fg = "black",
             ...)
{
    theme <-
        list(plot.polygon      = list(col = fill[1], border = fg[1]),
             box.rectangle     = list(col= symbol[1]),
             box.umbrella      = list(col= symbol[1]),
             dot.line          = list(col = reference),
             dot.symbol        = list(col = symbol[1]),
             plot.line         = list(col = symbol[1]),
             plot.symbol       = list(col= symbol[1]),
             regions           = list(col = colorRampPalette(region)(100)),
             reference.line    = list(col = reference),
             superpose.line    = list(col = symbol),
             superpose.symbol  = list(col = symbol),
             superpose.polygon = list(col = fill, border = fg),

             background        = list(col = bg),
             add.line          = list(col = fg),
             add.text          = list(col = fg),
             box.dot           = list(col = fg),
             axis.line         = list(col = fg),
             axis.text         = list(col = fg),
             strip.border      = list(col = fg),
             box.3d            = list(col = fg),
             par.xlab.text     = list(col = fg),
             par.ylab.text     = list(col = fg),
             par.zlab.text     = list(col = fg),
             par.main.text     = list(col = fg),
             par.sub.text      = list(col = fg))
    modifyList(modifyList(standard.theme("pdf"), theme), simpleTheme(...))
}

custom.theme.2 <- function(...)
{
    doit <-
        function(symbol = brewer.pal(n = 9, name = "Set1")[c(2:1, 3:5, 7:9)], ## blue first
                 fill = brewer.pal(n = 8, name = "Accent"),
                 region = brewer.pal(n = 11, name = "RdBu"),
                 ...)
        {
            custom.theme(symbol = symbol, fill = fill, region = region, ...)
        }
    doit(...)
}
