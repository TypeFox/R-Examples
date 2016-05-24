w <- gwindow("container", visible = FALSE)
size(w) <- c(500,400)

pg <- gpanedgroup(cont = w, horizontal=TRUE)

## left part of paned group
lg <- ggroup(horizontal = FALSE, cont = pg)

# gframe
gf <- gframe("frame label", cont = lg, expand=TRUE)
gbutton("button", cont = gf)

# svalue
names(gf)[1] <- "new frame label"

## expand
ge <- gexpandgroup("expand group", cont = lg, expand=TRUE)
gbutton("button", cont = ge)

# visible
visible(ge) <- TRUE


## right part of paned group

# notebook
nb <- gnotebook(cont = pg)

gbutton("button", label = "tab label", cont = nb)
gbutton("button", label = "tab label", cont = nb)
tab <- glayout(label = "tab layout", cont = nb)
tab[1,1] <- "label test"
tab[2,2] <- gedit("edit widget", cont = tab)
tab[3,1:2] <- gedit("expand two", cont = tab)
tab[4,1, anchor=c(-1,1)] <- "anchor"


# svalue
svalue(nb) <- 1

# names<-
names(nb)[2] <- "TAB LABEL"

size(w) <- c(400,400)
svalue(pg) <- 0.5

visible(w) <- TRUE
