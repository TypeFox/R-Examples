## ----echo = FALSE--------------------------------------------------------
library(knitr)
opts_chunk$set(fig.pos = "")

library(circlize)
chordDiagram = function(...) {
    circos.par(unit.circle.segments = 200)
    circlize::chordDiagram(...)
}

## ------------------------------------------------------------------------
mat = matrix(1:9, 3)
rownames(mat) = letters[1:3]
colnames(mat) = LETTERS[1:3]
mat

## ------------------------------------------------------------------------
df = data.frame(from = letters[1:3], to = LETTERS[1:3], value = 1:3)
df

## ------------------------------------------------------------------------
set.seed(999)
mat = matrix(sample(18, 18), 3, 6) 
rownames(mat) = paste0("S", 1:3)
colnames(mat) = paste0("E", 1:6)
mat

df = data.frame(from = rep(rownames(mat), times = ncol(mat)),
    to = rep(colnames(mat), each = nrow(mat)),
    value = as.vector(mat),
    stringsAsFactors = FALSE)
df = df[sample(18, 10), ]
df

## ----chord_diagram_basic_simple, eval = FALSE----------------------------
#  chordDiagram(mat)
#  circos.clear()

## ----eval = FALSE--------------------------------------------------------
#  # code is not run when building the vignette
#  chordDiagram(df)
#  circos.clear()

## ----chord_diagram_basic_gap_degree, eval = FALSE------------------------
#  circos.par(gap.degree = c(rep(2, nrow(mat)-1), 10, rep(2, ncol(mat)-1), 10))
#  chordDiagram(mat)
#  circos.clear()

## ----eval = FALSE--------------------------------------------------------
#  # code is not run when building the vignette
#  circos.par(gap.degree = c(rep(2, length(unique(df[[1]]))-1), 10,
#                            rep(2, length(unique(df[[2]]))-1), 10))
#  chordDiagram(df)
#  circos.clear()

## ----chord_diagram_basic_start_degree, eval = FALSE----------------------
#  circos.par(start.degree = 90)
#  chordDiagram(mat)
#  circos.clear()

## ----chord_diagram_basic_order, eval = FALSE-----------------------------
#  chordDiagram(mat, order = c("S1", "E1", "E2", "S2", "E3", "E4", "S3", "E5", "E6"))

## ----chord_diagram_basic, echo = FALSE, fig.align = "center", out.width = "\\textwidth", fig.cap = "Basic usages of {\\tt chordDiagram}. A) default style; B) set {\\tt gap.degree}; C) set {\\tt start.degree}; D) set orders of sectors."----
par(mfrow = c(2, 2))
chordDiagram(mat)
circos.clear()
text(-0.9, 0.9, "A", cex = 1.5)
circos.par(gap.degree = c(rep(2, nrow(mat)-1), 10, rep(2, ncol(mat)-1), 10))
chordDiagram(mat)
circos.clear()
text(-0.9, 0.9, "B", cex = 1.5)
circos.par(start.degree = 90)
chordDiagram(mat)
circos.clear()
text(-0.9, 0.9, "C", cex = 1.5)
chordDiagram(mat, order = c("S1", "E1", "E2", "S2", "E3", "E4", "S3", "E5", "E6"))
text(-0.9, 0.9, "D", cex = 1.5)
par(mfrow = c(1, 1))

## ----chord_diagram_color_grid, eval = FALSE------------------------------
#  grid.col = c(S1 = "red", S2 = "green", S3 = "blue",
#      E1 = "grey", E2 = "grey", E3 = "grey", E4 = "grey", E5 = "grey", E6 = "grey")
#  chordDiagram(mat, grid.col = grid.col)

## ----chord_diagram_color_transparency, eval = FALSE----------------------
#  chordDiagram(mat, grid.col = grid.col, transparency = 0)

## ----chord_diagram_color_mat, eval = FALSE-------------------------------
#  col_mat = rand_color(length(mat), transparency = 0.5)
#  dim(col_mat) = dim(mat)  # to make sure it is a matrix
#  chordDiagram(mat, grid.col = grid.col, col = col_mat)

## ----eval = FALSE--------------------------------------------------------
#  # code is not run when building the vignette
#  col = rand_color(nrow(df))
#  chordDiagram(df, grid.col = grid.col, col = col)

## ----chord_diagram_color_fun, eval = FALSE-------------------------------
#  col_fun = colorRamp2(range(mat), c("#FFEEEEEE", "#FF0000"))
#  chordDiagram(mat, grid.col = grid.col, col = col_fun)

## ----chord_diagram_color_row_col, eval = FALSE, echo = -2----------------
#  chordDiagram(mat, grid.col = grid.col, row.col = 1:3)
#  text(-0.9, 0.9, "E", cex = 1.5)
#  chordDiagram(mat, grid.col = grid.col, column.col = 1:6)

## ----chord_diagram_color, echo = FALSE, fig.align = "center", out.width = "0.6\\textheight", out.height = "0.9\\textheight", fig.width = 7, fig.height = 10.5, fig.cap = "Color settings in {\\tt chordDiagram}. A) set {\\tt grid.col}; B) set {\\tt transparency}; C) set {\\tt col} as a matrix; D) set {\\tt col} as a function; E) set {\\tt row.col}; F) set {\\tt column.col}."----
par(mfrow = c(3, 2))
grid.col = c(S1 = "red", S2 = "green", S3 = "blue",
    E1 = "grey", E2 = "grey", E3 = "grey", E4 = "grey", E5 = "grey", E6 = "grey")
chordDiagram(mat, grid.col = grid.col)
text(-0.9, 0.9, "A", cex = 1.5)
chordDiagram(mat, grid.col = grid.col, transparency = 0)
text(-0.9, 0.9, "B", cex = 1.5)
col_mat = rand_color(length(mat), transparency = 0.5)
dim(col_mat) = dim(mat)  # to make sure it is a matrix
chordDiagram(mat, grid.col = grid.col, col = col_mat)
text(-0.9, 0.9, "C", cex = 1.5)
col_fun = colorRamp2(range(mat), c("#FFEEEEEE", "#FF0000"))
chordDiagram(mat, grid.col = grid.col, col = col_fun)
text(-0.9, 0.9, "D", cex = 1.5)
chordDiagram(mat, grid.col = grid.col, row.col = 1:3)
text(-0.9, 0.9, "E", cex = 1.5)
chordDiagram(mat, grid.col = grid.col, column.col = 1:6)
text(-0.9, 0.9, "F", cex = 1.5)
par(mfrow = c(1, 1))

## ----chord_diagram_style_scalar, eval = FALSE----------------------------
#  chordDiagram(mat, grid.col = grid.col, link.lwd = 2, link.lty = 2, link.border = "black")

## ----chord_diagram_style_fullmat, eval = FALSE---------------------------
#  lwd_mat = matrix(1, nrow = nrow(mat), ncol = ncol(mat))
#  rownames(lwd_mat) = rownames(mat)
#  colnames(lwd_mat) = colnames(mat)
#  lwd_mat[mat > 12] = 2
#  
#  border_mat = matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
#  rownames(border_mat) = rownames(mat)
#  colnames(border_mat) = colnames(mat)
#  border_mat[mat > 12] = "black"
#  
#  chordDiagram(mat, grid.col = grid.col, link.lwd = lwd_mat, link.border = border_mat)

## ----chord_diagram_style_submatrix, eval = FALSE-------------------------
#  border_mat2 = matrix("black", nrow = 1, ncol = ncol(mat))
#  rownames(border_mat2) = rownames(mat)[2]
#  colnames(border_mat2) = colnames(mat)
#  
#  chordDiagram(mat, grid.col = grid.col, link.lwd = 2, link.border = border_mat2)

## ----chord_diagram_style_dataframe, eval = FALSE-------------------------
#  lty_df = data.frame(c("S1", "S2", "S3"), c("E5", "E6", "E4"), c(1, 2, 3))
#  lwd_df = data.frame(c("S1", "S2", "S3"), c("E5", "E6", "E4"), c(2, 2, 2))
#  border_df = data.frame(c("S1", "S2", "S3"), c("E5", "E6", "E4"), c(1, 1, 1))
#  chordDiagram(mat, grid.col = grid.col, link.lty = lty_df, link.lwd = lwd_df,
#      link.border = border_df)

## ----eval = FALSE--------------------------------------------------------
#  # code is not run when building the vignette
#  chordDiagram(df, grid.col = grid.col, link.lty = sample(1:3, nrow(df), replace = TRUE),
#      link.lwd = runif(nrow(df))*2, link.border = sample(0:1, nrow(df), replace = TRUE))

## ----chord_diagram_style, echo = FALSE, fig.align = "center", out.width = "\\textwidth", fig.cap = "Link style settings in {\\tt chordDiagram}. A) graphic parameters set as scalar; B) graphic parameters set as matrix; C) graphic parameters set as sub matrix. D) graphic parameters set as a three-column data frame."----
par(mfrow = c(2, 2))
chordDiagram(mat, grid.col = grid.col, link.lwd = 2, link.lty = 2, link.border = "black")
text(-0.9, 0.9, "A", cex = 1.5)
lwd_mat = matrix(1, nrow = nrow(mat), ncol = ncol(mat))
rownames(lwd_mat) = rownames(mat)
colnames(lwd_mat) = colnames(mat)
lwd_mat[mat > 12] = 2

border_mat = matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
rownames(border_mat) = rownames(mat)
colnames(border_mat) = colnames(mat)
border_mat[mat > 12] = "black"

chordDiagram(mat, grid.col = grid.col, link.lwd = lwd_mat, link.border = border_mat)
text(-0.9, 0.9, "B", cex = 1.5)
border_mat2 = matrix("black", nrow = 1, ncol = ncol(mat))
rownames(border_mat2) = rownames(mat)[2]
colnames(border_mat2) = colnames(mat)

chordDiagram(mat, grid.col = grid.col, link.lwd = 2, link.border = border_mat2)
text(-0.9, 0.9, "C", cex = 1.5)
lty_df = data.frame(c("S1", "S2", "S3"), c("E5", "E6", "E4"), c(1, 2, 3))
lwd_df = data.frame(c("S1", "S2", "S3"), c("E5", "E6", "E4"), c(2, 2, 2))
border_df = data.frame(c("S1", "S2", "S3"), c("E5", "E6", "E4"), c(1, 1, 1))
chordDiagram(mat, grid.col = grid.col, link.lty = lty_df, link.lwd = lwd_df,
    link.border = border_df)
text(-0.9, 0.9, "D", cex = 1.5)
par(mfrow = c(1, 1))

## ----chord_diagram_highlight_row, eval = FALSE---------------------------
#  chordDiagram(mat, grid.col = grid.col, row.col = c("#FF000080", "#00FF0010", "#0000FF10"))

## ----chord_diagram_highlight_mat, eval = FALSE---------------------------
#  col_mat[mat < 12] = "#00000000"
#  chordDiagram(mat, grid.col = grid.col, col = col_mat)

## ----chord_diagram_highlight_fun, eval = FALSE---------------------------
#  col_fun = function(x) ifelse(x < 12, "#00000000", "#FF000080")
#  chordDiagram(mat, grid.col = grid.col, col = col_fun)

## ----chord_diagram_highlight_df, eval = FALSE----------------------------
#  col_df = data.frame(c("S1", "S2", "S3"), c("E5", "E6", "E4"),
#      c("#FF000080", "#00FF0080", "#0000FF80"))
#  chordDiagram(mat, grid.col = grid.col, col = col_df)

## ----eval = FALSE--------------------------------------------------------
#  # code is not run when building the vignette
#  col = rand_color(nrow(df))
#  col[df[[3]] < 10] = "#00000000"
#  chordDiagram(df, grid.col = grid.col, col = col)

## ----chord_diagram_highlight, echo = FALSE, fig.align = "center", out.width = "\\textwidth", fig.cap = "Highlight links by colors. A) set {\\tt row.col}; B) set by matrix; C) set by color function; D) set by a three-column data frame."----
par(mfrow = c(2, 2))
chordDiagram(mat, grid.col = grid.col, row.col = c("#FF000080", "#00FF0010", "#0000FF10"))
text(-0.9, 0.9, "A", cex = 1.5)
col_mat[mat < 12] = "#00000000"
chordDiagram(mat, grid.col = grid.col, col = col_mat) 
text(-0.9, 0.9, "B", cex = 1.5)
col_fun = function(x) ifelse(x < 12, "#00000000", "#FF000080")
chordDiagram(mat, grid.col = grid.col, col = col_fun)
text(-0.9, 0.9, "C", cex = 1.5)
col_df = data.frame(c("S1", "S2", "S3"), c("E5", "E6", "E4"), 
    c("#FF000080", "#00FF0080", "#0000FF80"))
chordDiagram(mat, grid.col = grid.col, col = col_df)
text(-0.9, 0.9, "D", cex = 1.5)
par(mfrow = c(1, 1))

## ----chord_diagram_link_order1, eval = FALSE, echo = c(1, 3)-------------
#  chordDiagram(mat, grid.col = grid.col, link.border = 1,
#      link.sort = TRUE, link.decreasing = TRUE)
#  text(-0.9, 0.9, "A", cex = 1.5)
#  chordDiagram(mat, grid.col = grid.col, link.border = 1,
#      link.sort = TRUE, link.decreasing = FALSE)
#  text(-0.9, 0.9, "B", cex = 1.5)

## ----chord_diagram_link_order, echo = FALSE, fig.align = "center", out.width = "\\textwidth", out.height = "0.5\\textwidth", fig.width = 7, fig.height = 3.5, fig.cap = "Orders of links. A) set {\\tt link.decreasing} to {\\tt TRUE}; B) set {\\tt link.decreasing} to {\\tt FALSE}."----
par(mfrow = c(1, 2))
chordDiagram(mat, grid.col = grid.col, link.border = 1, 
    link.sort = TRUE, link.decreasing = TRUE)
text(-0.9, 0.9, "A", cex = 1.5)
chordDiagram(mat, grid.col = grid.col, link.border = 1, 
    link.sort = TRUE, link.decreasing = FALSE)
text(-0.9, 0.9, "B", cex = 1.5)
par(mfrow = c(1, 1))

## ----chord_diagram_directional_simple, eval = FALSE, echo = c(1, 3, 5)----
#  chordDiagram(mat, grid.col = grid.col, directional = 1)
#  text(-0.9, 0.9, "A", cex = 1.5)
#  chordDiagram(mat, grid.col = grid.col, directional = 1, diffHeight = 0.08)
#  text(-0.9, 0.9, "B", cex = 1.5)
#  chordDiagram(mat, grid.col = grid.col, directional = -1)
#  text(-0.9, 0.9, "C", cex = 1.5)

## ------------------------------------------------------------------------
mat2 = matrix(sample(100, 35), nrow = 5)
rownames(mat2) = letters[1:5]
colnames(mat2) = letters[1:7]
mat2

## ----chord_diagram_directional_overlap, eval = FALSE---------------------
#  chordDiagram(mat2, grid.col = 1:7, directional = 1, row.col = 1:5)

## ------------------------------------------------------------------------
mat3 = mat2
for(cn in intersect(rownames(mat3), colnames(mat3))) {
    mat3[cn, cn] = 0
}
mat3

## ----chord_diagram_directional_non_selfloop, eval = FALSE----------------
#  chordDiagram(mat3, grid.col = 1:7, directional = 1, row.col = 1:5)

## ----chord_diagram_directional_arrow, eval = FALSE-----------------------
#  arr.col = data.frame(c("S1", "S2", "S3"), c("E5", "E6", "E4"),
#      c("black", "black", "black"))
#  chordDiagram(mat, grid.col = grid.col, directional = 1, direction.type = "arrows",
#      link.arr.col = arr.col, link.arr.length = 0.2)

## ----chord_diagram_directional_arrow2, eval = FALSE----------------------
#  arr.col = data.frame(c("S1", "S2", "S3"), c("E5", "E6", "E4"),
#      c("black", "black", "black"))
#  chordDiagram(mat, grid.col = grid.col, directional = 1,
#      direction.type = c("diffHeight", "arrows"),
#      link.arr.col = arr.col, link.arr.length = 0.2)

## ----chord_diagram_directional_arrow3, eval = FALSE----------------------
#  matx = matrix(rnorm(64), 8)
#  chordDiagram(matx, directional = 1, direction.type = c("diffHeight", "arrows"),
#      link.arr.type = "big.arrow")

## ----chord_diagram_directional_arrow4, eval = FALSE----------------------
#  chordDiagram(matx, directional = 1, direction.type = c("diffHeight", "arrows"),
#      link.arr.type = "big.arrow", diffHeight = -0.04)

## ----eval = FALSE--------------------------------------------------------
#  # code is not run when building the vignette
#  chordDiagram(df, directional = 1)

## ----chord_diagram_directional, echo = FALSE, fig.align = "center", out.width = "\\textwidth", out.height = "\\textwidth", fig.width = 10.5, fig.height = 10.5, fig.cap = "Visualization of directional matrix. A) with default settings; B) set difference of two feet of links; C) set the starting feet; D, E) row names and column names have overlaps; F, G, H, I) directions are represented by arrows."----
par(mfrow = c(3, 3))
chordDiagram(mat, grid.col = grid.col, directional = 1)
text(-0.9, 0.9, "A", cex = 1.5)
chordDiagram(mat, grid.col = grid.col, directional = 1, diffHeight = 0.08)
text(-0.9, 0.9, "B", cex = 1.5)
chordDiagram(mat, grid.col = grid.col, directional = -1)
text(-0.9, 0.9, "C", cex = 1.5)
chordDiagram(mat2, grid.col = 1:7, directional = 1, row.col = 1:5)
text(-0.9, 0.9, "D", cex = 1.5)
chordDiagram(mat3, grid.col = 1:7, directional = 1, row.col = 1:5)
text(-0.9, 0.9, "E", cex = 1.5)
arr.col = data.frame(c("S1", "S2", "S3"), c("E5", "E6", "E4"), 
    c("black", "black", "black"))
chordDiagram(mat, grid.col = grid.col, directional = 1, direction.type = "arrows",
    link.arr.col = arr.col, link.arr.length = 0.2)
text(-0.9, 0.9, "F", cex = 1.5)
arr.col = data.frame(c("S1", "S2", "S3"), c("E5", "E6", "E4"), 
    c("black", "black", "black"))
chordDiagram(mat, grid.col = grid.col, directional = 1, 
    direction.type = c("diffHeight", "arrows"),
    link.arr.col = arr.col, link.arr.length = 0.2)
text(-0.9, 0.9, "G", cex = 1.5)
matx = matrix(rnorm(64), 8)
chordDiagram(matx, directional = 1, direction.type = c("diffHeight", "arrows"),
    link.arr.type = "big.arrow")
text(-0.9, 0.9, "H", cex = 1.5)
chordDiagram(matx, directional = 1, direction.type = c("diffHeight", "arrows"),
    link.arr.type = "big.arrow", diffHeight = -0.04)
text(-0.9, 0.9, "I", cex = 1.5)
par(mfrow = c(1, 1))

## ----chord_diagram_self_link1, eval = FALSE, echo = c(1,2,4)-------------
#  df2 = data.frame(start = c("a", "b", "c", "a"), end = c("a", "a", "b", "c"))
#  chordDiagram(df2, grid.col = 1:3, self.link = 1, link.border = 1)
#  text(-0.9, 0.9, "A", cex = 1.5)
#  chordDiagram(df2, grid.col = 1:3, self.link = 2, link.border = 1)
#  text(-0.9, 0.9, "B", cex = 1.5)

## ----chord_diagram_self_link, echo = FALSE, fig.align = "center", out.width = "\\textwidth", out.height = "0.5\\textwidth", fig.width = 7, fig.height = 3.5, fig.cap = "Deal with self links. A) set {\\tt self.link} to 1; B) set {\\tt self.link} to 2."----
par(mfrow = c(1, 2))
df2 = data.frame(start = c("a", "b", "c", "a"), end = c("a", "a", "b", "c"))
chordDiagram(df2, grid.col = 1:3, self.link = 1, link.border = 1)
text(-0.9, 0.9, "A", cex = 1.5)
chordDiagram(df2, grid.col = 1:3, self.link = 2, link.border = 1)
text(-0.9, 0.9, "B", cex = 1.5)
par(mfrow = c(1, 1))

## ----chord_diagram_symmetric_show, eval = FALSE--------------------------
#  mat3 = matrix(rnorm(25), 5)
#  colnames(mat3) = letters[1:5]
#  chordDiagram(cor(mat3), grid.col = 1:5, symmetric = TRUE,
#      col = colorRamp2(c(-1, 0, 1), c("green", "white", "red")))

## ----chord_diagram_symmetric_hidden, eval = FALSE, echo = FALSE----------
#  chordDiagram(cor(mat3), grid.col = 1:5, col = colorRamp2(c(-1, 0, 1), c("green", "white", "red")))

## ----chord_diagram_symmetric, echo = FALSE, fig.align = "center", out.width = "\\textwidth", out.height = "0.5\\textwidth", fig.width = 7, fig.height = 3.5, fig.cap = "Visualization of symmetric matrix. A) set {\\tt symmetric} to {\\tt TRUE}; B) set {\\tt symmetric} to {\\tt FALSE}."----
par(mfrow = c(1, 2))
mat3 = matrix(rnorm(25), 5)
colnames(mat3) = letters[1:5]
chordDiagram(cor(mat3), grid.col = 1:5, symmetric = TRUE,
    col = colorRamp2(c(-1, 0, 1), c("green", "white", "red")))
text(-0.9, 0.9, "A", cex = 1.5)
chordDiagram(cor(mat3), grid.col = 1:5, col = colorRamp2(c(-1, 0, 1), c("green", "white", "red")))
text(-0.9, 0.9, "B", cex = 1.5)
par(mfrow = c(1, 1))

## ----echo = 2:3----------------------------------------------------------
pdf(NULL)
chordDiagram(mat)
circos.info()
invisible(dev.off())

## ----chord_diagram_default_track_simple, eval = FALSE, echo = c(1, 3, 4, 6)----
#  chordDiagram(mat, grid.col = grid.col, annotationTrack = "grid")
#  text(-0.9, 0.9, "A", cex = 1.5)
#  chordDiagram(mat, grid.col = grid.col, annotationTrack = c("name", "grid"),
#      annotationTrackHeight = c(0.03, 0.01))
#  text(-0.9, 0.9, "B", cex = 1.5)
#  chordDiagram(mat, grid.col = grid.col, annotationTrack = NULL)
#  text(-0.9, 0.9, "C", cex = 1.5)

## ----chord_diagram_default_track, echo = FALSE, fig.align = "center", out.width = "\\textwidth", fig.cap = "Track organization in {\\tt chordDiagram}. A) only show the grid track; B) set label track and grid track with heights; C) do not add label track or grid track."----
par(mfrow = c(2, 2))
chordDiagram(mat, grid.col = grid.col, annotationTrack = "grid")
text(-0.9, 0.9, "A", cex = 1.5)
chordDiagram(mat, grid.col = grid.col, annotationTrack = c("name", "grid"),
    annotationTrackHeight = c(0.03, 0.01))
text(-0.9, 0.9, "B", cex = 1.5)
chordDiagram(mat, grid.col = grid.col, annotationTrack = NULL)
text(-0.9, 0.9, "C", cex = 1.5)
par(mfrow = c(1, 1))

## ----echo = 2:3----------------------------------------------------------
pdf(NULL)
chordDiagram(mat, preAllocateTracks = 2)
circos.info()
invisible(dev.off())

## ----eval = FALSE--------------------------------------------------------
#  list(ylim = c(0, 1),
#       track.height = circos.par("track.height"),
#       bg.col = NA,
#       bg.border = NA,
#       bg.lty = par("lty"),
#       bg.lwd = par("lwd"))

## ----eval = FALSE--------------------------------------------------------
#  chordDiagram(mat, annotationTrack = NULL,
#      preAllocateTracks = list(track.height = 0.3))
#  circos.info(sector.index = "S1", track.index = 1)

## ----eval = FALSE--------------------------------------------------------
#  chordDiagram(mat, annotationTrack = NULL,
#      preAllocateTracks = list(list(track.height = 0.1),
#                               list(bg.border = "black")))

## ----chord_diagram_labels_show, eval = FALSE-----------------------------
#  chordDiagram(mat, grid.col = grid.col, annotationTrack = "grid",
#      preAllocateTracks = list(track.height = 0.3))
#  # we go back to the first track and customize sector labels
#  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
#      xlim = get.cell.meta.data("xlim")
#      ylim = get.cell.meta.data("ylim")
#      sector.name = get.cell.meta.data("sector.index")
#      circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
#          niceFacing = TRUE, adj = c(0, 0.5))
#  }, bg.border = NA) # here set bg.border to NA is important

## ----chord_diagram_labels_inside, eval = FALSE---------------------------
#  chordDiagram(mat, grid.col = grid.col,
#      annotationTrack = "grid", annotationTrackHeight = 0.15)
#  for(si in get.all.sector.index()) {
#      xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
#      ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
#      circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1,
#          facing = "bending.inside", col = "white")
#  }

## ----chord_diagram_labels_multile_style, eval = FALSE--------------------
#  mat2 = matrix(rnorm(100), 10)
#  chordDiagram(mat2, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1))
#  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
#      xlim = get.cell.meta.data("xlim")
#      xplot = get.cell.meta.data("xplot")
#      ylim = get.cell.meta.data("ylim")
#      sector.name = get.cell.meta.data("sector.index")
#  
#      if(abs(xplot[2] - xplot[1]) < 20) {
#          circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
#              niceFacing = TRUE, adj = c(0, 0.5))
#      } else {
#          circos.text(mean(xlim), ylim[1], sector.name, facing = "inside",
#              niceFacing = TRUE, adj = c(0.5, 0))
#      }
#  }, bg.border = NA)

## ----chord_diagram_labels, echo = FALSE, fig.align = "center", out.width = "\\textwidth", fig.cap = "Customize sector labels. A) put sector labels in radical direction; B) sector labels are put inside grids; C) sector labels are put in different direction according the width of sectors."----
par(mfrow = c(2, 2))
chordDiagram(mat, grid.col = grid.col, annotationTrack = "grid", 
    preAllocateTracks = list(track.height = 0.3))
# we go back to the first track and customize sector labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", 
        niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
text(-0.9, 0.9, "A", cex = 1.5)
chordDiagram(mat, grid.col = grid.col, 
    annotationTrack = "grid", annotationTrackHeight = 0.15)
for(si in get.all.sector.index()) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, 
        facing = "bending.inside", col = "white")
}
text(-0.9, 0.9, "B", cex = 1.5)
circos.par(points.overflow.warning = FALSE)
mat2 = matrix(rnorm(100), 10)
chordDiagram(mat2, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")

    if(abs(xplot[2] - xplot[1]) < 20) {
        circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
            niceFacing = TRUE, adj = c(0, 0.5))
    } else {
        circos.text(mean(xlim), ylim[1], sector.name, facing = "inside", 
            niceFacing = TRUE, adj = c(0.5, 0))
    }
}, bg.border = NA)
circos.clear()
text(-0.9, 0.9, "C", cex = 1.5)
par(mfrow = c(1, 1))

## ----chord_diagram_axes_two, eval = FALSE--------------------------------
#  # similar as the previous example, but we only plot the grid track
#  chordDiagram(mat, grid.col = grid.col, annotationTrack = "grid",
#      preAllocateTracks = list(track.height = 0.1))
#  for(si in get.all.sector.index()) {
#      circos.axis(h = "top", labels.cex = 0.3, major.tick.percentage = 0.2,
#          sector.index = si, track.index = 2)
#  }
#  
#  # the second axis as well as the sector labels are added in this track
#  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
#      xlim = get.cell.meta.data("xlim")
#      xplot = get.cell.meta.data("xplot")
#      ylim = get.cell.meta.data("ylim")
#      sector.name = get.cell.meta.data("sector.index")
#  
#      if(abs(xplot[2] - xplot[1]) > 20) {
#          circos.lines(xlim, c(mean(ylim), mean(ylim)), lty = 3) # dotted line
#          for(p in seq(0.2, 1, by = 0.2)) {
#              circos.text(p*(xlim[2] - xlim[1]) + xlim[1], mean(ylim) + 0.1,
#                  p, cex = 0.3, adj = c(0.5, 0), niceFacing = TRUE)
#          }
#      }
#      circos.text(mean(xlim), 1, sector.name, niceFacing = TRUE, adj = c(0.5, 0))
#  }, bg.border = NA)
#  circos.clear()

## ----chord_diagram_axes, echo = FALSE, fig.align = "center", out.width = "0.5\\textwidth", out.height = "0.5\\textwidth", fig.width = 3.5, fig.height = 3.5, fig.cap = "Customize sector axes"----
# similar as the previous example, but we only plot the grid track
chordDiagram(mat, grid.col = grid.col, annotationTrack = "grid", 
    preAllocateTracks = list(track.height = 0.1))
for(si in get.all.sector.index()) {
    circos.axis(h = "top", labels.cex = 0.3, major.tick.percentage = 0.2,
        sector.index = si, track.index = 2)
}

# the second axis as well as the sector labels are added in this track
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    
    if(abs(xplot[2] - xplot[1]) > 20) {
        circos.lines(xlim, c(mean(ylim), mean(ylim)), lty = 3) # dotted line
        for(p in seq(0.2, 1, by = 0.2)) {
            circos.text(p*(xlim[2] - xlim[1]) + xlim[1], mean(ylim) + 0.1, 
                p, cex = 0.3, adj = c(0.5, 0), niceFacing = TRUE)
        }
    }
    circos.text(mean(xlim), 1, sector.name, niceFacing = TRUE, adj = c(0.5, 0))
}, bg.border = NA)
circos.clear()

## ----chord_diagram_compare_1, eval = FALSE-------------------------------
#  mat1 = matrix(sample(20, 25, replace = TRUE), 5)
#  
#  gap.degree = c(rep(2, 4), 10, rep(2, 4), 10)
#  circos.clear()
#  circos.par(gap.degree = gap.degree, start.degree = -10/2)
#  chordDiagram(mat1, directional = 1, grid.col = rep(1:5, 2))
#  circos.clear()

## ----chord_diagram_compare_2, eval = FALSE-------------------------------
#  mat2 = mat1 / 2

## ----chord_diagram_compare_3, eval = FALSE-------------------------------
#  percent = sum(abs(mat2)) / sum(abs(mat1))
#  blank.degree = (360 - sum(gap.degree)) * (1 - percent)

## ----chord_diagram_compare_4, eval = FALSE-------------------------------
#  big.gap = (blank.degree - sum(rep(2, 8)))/2
#  gap.degree = c(rep(2, 4), big.gap, rep(2, 4), big.gap)
#  circos.par(gap.degree = gap.degree, start.degree = -big.gap/2)
#  chordDiagram(mat2, directional = 1, grid.col = rep(1:5, 2), transparency = 0.5)
#  circos.clear()

## ----chord_diagram_compare, echo = FALSE, fig.align = "center", out.width = "\\textwidth", out.height = "0.5\\textwidth", fig.width = 7, fig.height = 3.5, fig.cap = "Compare two Chord Diagrams and make them in same scale. bottom matrix has half the values as in the upper matrix."----
par(mfrow = c(1, 2))
mat1 = matrix(sample(20, 25, replace = TRUE), 5)

gap.degree = c(rep(2, 4), 10, rep(2, 4), 10)
circos.clear()
circos.par(gap.degree = gap.degree, start.degree = -10/2)
chordDiagram(mat1, directional = 1, grid.col = rep(1:5, 2))
circos.clear()
text(-0.9, 0.9, "A", cex = 1.5)
mat2 = mat1 / 2
percent = sum(abs(mat2)) / sum(abs(mat1))
blank.degree = (360 - sum(gap.degree)) * (1 - percent)
big.gap = (blank.degree - sum(rep(2, 8)))/2
gap.degree = c(rep(2, 4), big.gap, rep(2, 4), big.gap)
circos.par(gap.degree = gap.degree, start.degree = -big.gap/2)
chordDiagram(mat2, directional = 1, grid.col = rep(1:5, 2), transparency = 0.5)
circos.clear()
text(-0.9, 0.9, "B", cex = 1.5)
par(mfrow = c(1, 1))

## ------------------------------------------------------------------------
mat = matrix(rnorm(36), 6, 6)
rownames(mat) = paste0("R", 1:6)
colnames(mat) = paste0("C", 1:6)
mat[2, ] = 1e-10
mat[, 3] = 1e-10

## ----chord_diagram_reduce_1, eval=FALSE----------------------------------
#  chordDiagram(mat)

## ----chord_diagram_reduce_2, eval=FALSE----------------------------------
#  chordDiagram(mat, row.col = rep(c("red", "blue"), 3))

## ----chord_diagram_reduce_3, eval=FALSE----------------------------------
#  chordDiagram(mat, grid.col = rep(c("red", "blue"), 6))
#  circos.clear()

## ----chord_diagram_reduce_4, eval=FALSE----------------------------------
#  circos.par("gap.degree" = rep(c(2, 10), 6))
#  chordDiagram(mat)
#  circos.clear()

## ----chord_diagram_reduce, echo = FALSE, fig.align = "center", out.width = "\\textwidth", fig.cap = "Reduced Chord Diagram with removing tiny sectors. A) notice how sector labels are reduced; B) notice how link colors are reduced; C) notice how grid colors are reduced; D) notice how gap degrees are reduced."----
par(mfrow = c(2, 2))
chordDiagram(mat)
text(-0.9, 0.9, "A", cex = 1.5)
chordDiagram(mat, row.col = rep(c("red", "blue"), 3))
text(-0.9, 0.9, "B", cex = 1.5)
chordDiagram(mat, grid.col = rep(c("red", "blue"), 6))
circos.clear()
text(-0.9, 0.9, "C", cex = 1.5)
circos.par("gap.degree" = rep(c(2, 10), 6))
chordDiagram(mat)
circos.clear()
text(-0.9, 0.9, "D", cex = 1.5)

