library(grid)

##################
# Testing
##################

# Test of a combination of shape and plotmath
library(gridSVG)
grid.newpage()
grid.rect(x=.4, width=.4, height=.8, gp=gpar(fill="grey", col=NA))
grid.text(expression((a + b)^2), x=.4, y=.6)
grid.export("test2.svg")

test <- function(e, file, cex=1) {
    require(grid)
    require(gridSVG)
    grid.newpage()
    grid.rect()
    nrowcol <- n2mfrow(length(e))
    pushViewport(viewport(layout=grid.layout(nrowcol[1], nrowcol[2])))
    for (i in 1:length(e)) {
        col <- (i - 1) %% nrowcol[2] + 1
        row <- (i - 1) %/% nrowcol[2] + 1
        pushViewport(viewport(layout.pos.row=row, layout.pos.col=col))
        grid.segments(c(0, .5), c(.5, 0), c(1, .5), c(.5, 1),
                      gp=gpar(col="grey"))
        just <- c("left", "bottom")
        tg <- textGrob(e[i], gp=gpar(cex=cex[(i - 1) %% length(cex) + 1]),
                       just=just)
        grid.rect(width=grobWidth(tg), height=grobHeight(tg),
                  just=just, gp=gpar(col="red"))
        grid.draw(tg)
        popViewport()
    }
    grid.export(file)
    grid.refresh()
}

test(expression(
    x + y,
    x + y,
    x + y,
    x - y,
    x / y,
    x * y,
    x %+-% y,
    x %/% y,
    x %*% y,
    x %.% y,
    x[i],
    x^2,
    paste(x, y, z),
    sqrt(x),
    sqrt(x, y),
    X + Y,
    a * b,
    a == b,
    x != y,
    x < y,
    x <= y,
    x > y,
    x >= y,
    x %~~% y,
    x %=~% y,
    x %==% y,
    x %prop% y,
    bolditalic(x + plain(y + italic(z + bold(k)))),
    , # space
    symbol("m") + symbol("\042"),
    list(x, y, z),
    1*...*n,
    1*cdots*n,
    1*ldots*n,
    x %subset% y,
    x %subseteq% y,
    x %notsubset% y,
    x %supset% y,
    x %supseteq% y,
    x %in% y,
    x %notin% y),
     "testml-1.svg",
     cex <- c(1, 1.5, .5, rep(1, 100)))

test(expression(
    hat(x),
    hat(xyz),
    tilde(x),
    dot(x),
    ring(x),
    bar(x),
    widehat(xyz),
    widetilde(xyz),
    x %<->% y,
    x %->% y,
    x %<-% y,
    x %up% y,
    x %down% y,
    x %<=>% y,
    x %=>% y,
    x %<=% y,
    x %dblup% y,
    x %dbldown% y,
    alpha - omega,
    Alpha - Omega,
    theta1*phi1*sigma1*omega1,
    Upsilon1,
    aleph,
    infinity,
    partialdiff,
    nabla,
    32*degree,
    60*minute,
    30*second,
    displaystyle(sum(x[i], i==1, n)),
    textstyle(sum(x[i], i==1, n)),
    textstyle(x[i]),
    scriptstyle(x[i]),
    scriptscriptstyle(x[i]),
    underline(xyz),
    x ~~ y,
    x + phantom(0) + y,
    x + over(1, phantom(0)),
    frac(x, y),
    over(x, y),
    atop(x, y),
    ),
     "testml-2.svg")

test(expression(
    sum(x[i], i==1, n),
    prod(plain(P)(X==x), x),
    integral(f(x)*dx, a, b),
    union(A[i], i==1, n),
    intersect(A[i], i==1, n),
    lim(f(x), x %->% 0),
    min(g(x), x > 0),
    inf(S),
    sup(S),
    x^y + z,
    x^(y + z),
    x^{y + z},
    group("(", list(a, b), "]"),
    bgroup("(", atop(x, y), ")"),
    group(lceil, x, rceil),
    ),
     "testml-3.svg")

