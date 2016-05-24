## ----echo = FALSE--------------------------------------------------------
library(knitr)
opts_chunk$set(fig.pos = "", fig.align = "center")

library(circlize)
circos.genomicInitialize = function(...) {
    circos.par(unit.circle.segments = 200)
    circlize::circos.genomicInitialize(...)
}

circos.initializeWithIdeogram = function(...) {
    circos.par(unit.circle.segments = 200)
    circlize::circos.initializeWithIdeogram(...)
}

## ------------------------------------------------------------------------
library(circlize)
set.seed(999)

bed = generateRandomBed()
head(bed)
bed = generateRandomBed(nr = 200, nc = 4)
nrow(bed)
bed = generateRandomBed(nc = 2, fun = function(k) runif(k))
head(bed)

## ----genomic_initialize_ideogram_1, eval = FALSE, echo = 1:2-------------
#  par(mar = c(1, 1, 1, 1))
#  circos.initializeWithIdeogram()
#  text(0, 0, "default", cex = 0.7)
#  text(-0.9, 0.9, "A", cex = 1.5)

## ----echo = 2:6----------------------------------------------------------
pdf(NULL)
cytoband.file = paste0(system.file(package = "circlize"), "/extdata/cytoBand.txt")
circos.initializeWithIdeogram(cytoband.file)
cytoband.df = read.table(cytoband.file, colClasses = c("character", "numeric",
    "numeric", "character", "character"), sep = "\t")
circos.initializeWithIdeogram(cytoband.df)
circos.clear()
invisible(dev.off())

## ----eval = FALSE--------------------------------------------------------
#  circos.initializeWithIdeogram(species = "hg18")
#  circos.initializeWithIdeogram(species = "mm10")

## ----genomic_initialize_ideogram_2, eval = FALSE, echo = 1---------------
#  circos.initializeWithIdeogram(chromosome.index = paste0("chr", 10:1))
#  text(0, 0, "subset of chromosomes", cex = 0.7)
#  text(-0.9, 0.9, "B", cex = 1.5)

## ----genomic_initialize_ideogram_3, eval = FALSE, echo = c(1:2, 6:7, 11:12)----
#  cytoband = cytoband.df
#  circos.initializeWithIdeogram(cytoband, sort.chr = TRUE)
#  text(0, 0, "read from cytoband df\nsort.chr = TRUE", cex = 0.7)
#  text(-0.9, 0.9, "C", cex = 1.5)
#  
#  cytoband = cytoband.df
#  circos.initializeWithIdeogram(cytoband, sort.chr = FALSE)
#  text(0, 0, "read from a data frame\nunique(cytoband[[1]])", cex = 0.7)
#  text(-0.9, 0.9, "D", cex = 1.5)
#  
#  cytoband[[1]] = factor(cytoband[[1]], levels = paste0("chr", c(22:1, "X", "Y")))
#  circos.initializeWithIdeogram(cytoband, sort.chr = FALSE)
#  text(0, 0, "read from cytoband file\nfirst column converted to factor\nlevels = paste0('chr', c(22:1, 'X', 'Y'))", cex = 0.7)
#  text(-0.9, 0.9, "E", cex = 1.5)

## ----eval=FALSE----------------------------------------------------------
#  cytoband = read.cytoband()
#  cytoband = read.cytoband(file)
#  cytoband = read.cytoband(df)
#  cytoband = read.cytoband(species)

## ----genomic_initialize_ideogram_4, eval = FALSE, echo = c(1, 4, 8)------
#  circos.initializeWithIdeogram(plotType = c("axis", "labels"))
#  text(0, 0, "plotType = c('axis', 'labels')", cex = 0.7)
#  text(-0.9, 0.9, "F", cex = 1.5)
#  circos.initializeWithIdeogram(plotType = NULL)
#  text(0, 0, "plotType = NULL", cex = 0.7)
#  text(-0.9, 0.9, "G", cex = 1.5)
#  
#  circos.clear()

## ----genomic_initialize_ideogram_5, eval = FALSE, echo = c(1:3, 7:9)-----
#  circos.par("start.degree" = 90)
#  circos.initializeWithIdeogram()
#  circos.clear()
#  text(0, 0, "'start.degree' = 90", cex = 0.7)
#  text(-0.9, 0.9, "H", cex = 1.5)
#  
#  circos.par("gap.degree" = rep(c(2, 4), 12))
#  circos.initializeWithIdeogram()
#  circos.clear()
#  text(0, 0, "'gap.degree' = rep(c(2, 4), 12)", cex = 0.7)
#  text(-0.9, 0.9, "I", cex = 1.5)

## ----genomic_customize_ideogram, out.width = "0.8\\textwidth", fig.cap = "Customize ideogram."----
par(mar = c(1, 1, 1, 1))
circos.initializeWithIdeogram(plotType = NULL)
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
    chr = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.rect(xlim[1], 0, xlim[2], 0.5, col = rand_color(1))
    circos.text(mean(xlim), 0.9, chr, cex = 0.5, facing = "clockwise", 
        niceFacing = TRUE)
}, bg.border = NA)
circos.clear()

## ----eval = FALSE--------------------------------------------------------
#  df = data.frame(
#      name  = c("TP53",  "TP63",    "TP73"),
#      start = c(7565097, 189349205, 3569084),
#      end   = c(7590856, 189615068, 3652765))
#  circos.genomicInitialize(df)

## ----eval=FALSE----------------------------------------------------------
#  circos.genomicInitialize(df)
#  circos.genomicInitialize(df, sector.names = c("tp53", "tp63", "tp73"))
#  circos.genomicInitialize(df, plotType)
#  
#  circos.par(gap.degree = 2)
#  circos.genomicInitialize(df)

## ------------------------------------------------------------------------
load(paste0(system.file(package = "circlize"), "/extdata/tp_family.RData"))
names(tp_family)
names(tp_family[["TP53"]])
head(tp_family[["TP53"]][[1]])

df = data.frame(gene = names(tp_family),
                start = sapply(tp_family, function(x) min(unlist(x))),
                end = sapply(tp_family, function(x) max(unlist(x))))
df

## ----genomic_gene_model_1, eval = FALSE----------------------------------
#  circos.genomicInitialize(df)
#  circos.genomicTrackPlotRegion(ylim = c(0, 1),
#      bg.col = c("#FF000040", "#00FF0040", "#0000FF40"),
#      bg.border = NA, track.height = 0.05)

## ----genomic_gene_model_2, eval = FALSE----------------------------------
#  n = max(sapply(tp_family, length))
#  circos.genomicTrackPlotRegion(ylim = c(0.5, n + 0.5),
#      panel.fun = function(region, value, ...) {
#          gn = get.cell.meta.data("sector.index")
#          tr = tp_family[[gn]]  # all transcripts for this gene
#          for(i in seq_along(tr)) {
#              # for each transcript
#              current_tr_start = min(tr[[i]]$start)
#              current_tr_end = max(tr[[i]]$end)
#              circos.lines(c(current_tr_start, current_tr_end),
#                  c(n - i, n - i), col = "#CCCCCC")
#              circos.genomicRect(tr[[i]], ytop = n - i + 0.4,
#                  ybottom = n - i - 0.4, col = "orange", border = NA)
#          }
#  }, bg.border = NA, track.height = 0.3)
#  circos.clear()

## ----genomic_gene_model, echo = FALSE, out.width = '\\textwidth', fig.cap = "Alternative transcripts for genes."----
par(mar = c(1, 1, 1, 1))
circos.genomicInitialize(df)
circos.genomicTrackPlotRegion(ylim = c(0, 1), 
    bg.col = c("#FF000040", "#00FF0040", "#0000FF40"), 
    bg.border = NA, track.height = 0.05)
n = max(sapply(tp_family, length))
circos.genomicTrackPlotRegion(ylim = c(0.5, n + 0.5), 
    panel.fun = function(region, value, ...) {
        gn = get.cell.meta.data("sector.index")
        tr = tp_family[[gn]]  # all transcripts for this gene
        for(i in seq_along(tr)) {
            # for each transcript
            current_tr_start = min(tr[[i]]$start)
            current_tr_end = max(tr[[i]]$end)
            circos.lines(c(current_tr_start, current_tr_end), 
                c(n - i, n - i), col = "#CCCCCC")
            circos.genomicRect(tr[[i]], ytop = n - i + 0.4, 
                ybottom = n - i - 0.4, col = "orange", border = NA)
        }
}, bg.border = NA, track.height = 0.3)
circos.clear()

## ------------------------------------------------------------------------
cytoband = read.cytoband()
df = cytoband$df
chromosome = cytoband$chromosome

# copy regions for the two zoomed chromosomes
zoom_df = df[df[[1]] %in% chromosome[1:2], ]
zoom_df[[1]] = paste0("zoom_", zoom_df[[1]])
df2 = rbind(df, zoom_df)

# attach ranges for two zoomed chromosomes
xrange = c(cytoband$chr.len, cytoband$chr.len[1:2])
normal_sector_index = seq_along(chromosome)
zoomed_sector_index = length(chromosome) + 1:2

# normalize in normal chromsomes and zoomed chromosomes separately
sector.width = c(xrange[normal_sector_index] / sum(xrange[normal_sector_index]), 
                 xrange[zoomed_sector_index] / sum(xrange[zoomed_sector_index])) 

## ----genomic_zoom_1, eval = FALSE----------------------------------------
#  par(mar = c(1, 1, 1, 1))
#  circos.par(start.degree = 90)
#  circos.genomicInitialize(df2, sector.width = sector.width)

## ------------------------------------------------------------------------
extend_zoomed_chromosome_in_bed = function(bed, chromosome, prefix = "zoom_") {
    zoom_bed = bed[bed[[1]] %in% chromosome, , drop = FALSE]
    zoom_bed[[1]] = paste0(prefix, zoom_bed[[1]])
    rbind(bed, zoom_bed)
}

## ----genomic_zoom_2, eval = FALSE----------------------------------------
#  bed = generateRandomBed(100)
#  circos.genomicTrackPlotRegion(extend_zoomed_chromosome_in_bed(bed, chromosome[1:2]),
#      panel.fun = function(region, value, ...) {
#          circos.genomicPoints(region, value, pch = 16, cex = 0.5)
#  })

## ----genomic_zoom_3, eval = FALSE----------------------------------------
#  circos.link("chr1", get.cell.meta.data("cell.xlim", sector.index = "chr1"),
#      "zoom_chr1", get.cell.meta.data("cell.xlim", sector.index = "zoom_chr1"),
#      col = "#00000020", border = NA)
#  circos.clear()

## ----genomic_zoom, echo = FALSE, out.width = "0.8\\textwidth", fig.cap = "Zoom chromosomes."----
par(mar = c(1, 1, 1, 1))
circos.par(start.degree = 90)
circos.genomicInitialize(df2, sector.width = sector.width)
bed = generateRandomBed(100)
circos.genomicTrackPlotRegion(extend_zoomed_chromosome_in_bed(bed, chromosome[1:2]),
    panel.fun = function(region, value, ...) {
        circos.genomicPoints(region, value, pch = 16, cex = 0.5)
})
circos.link("chr1", get.cell.meta.data("cell.xlim", sector.index = "chr1"),
    "zoom_chr1", get.cell.meta.data("cell.xlim", sector.index = "zoom_chr1"),
    col = "#00000020", border = NA)
circos.clear()

## ----eval=FALSE----------------------------------------------------------
#  circos.genomicTrackPlotRegion(data, panel.fun = function(region, value, ...) {
#      circos.genomicPoints(region, value, ...)
#  })

## ----echo = 2:6----------------------------------------------------------
pdf(NULL)
bed = generateRandomBed(nc = 1)
head(bed)
circos.initializeWithIdeogram(plotType = NULL)
circos.genomicTrackPlotRegion(bed, panel.fun = function(region, value, ...) {
    if(get.cell.meta.data("sector.index") == "chr1") {
        print(head(region))
        print(head(value))
    }
})
circos.clear()
invisible(dev.off())

## ----eval=FALSE----------------------------------------------------------
#  circos.genomicTrackPlotRegion(data, ylim = c(0, 1),
#      panel.fun = function(region, value, ...) {
#          circos.genomicPoints(region, value, ...)
#  })
#  circos.genomicTrackPlotRegion(data, numeric.column,
#      panel.fun = function(region, value, ...) {
#          circos.genomicPoints(region, value, ...)
#  })

## ----eval=FALSE----------------------------------------------------------
#  circos.genomicPoints(region, value, ...)
#  circos.genomicPoints(region, value, numeric.column = c(1, 2))
#  circos.genomicPoints(region, value, cex, pch)
#  circos.genomicPoints(region, value, sector.index, track.index)

## ----eval = FALSE--------------------------------------------------------
#  circos.genomicPoints = function(region, value, numeric.column = 1, ...) {
#      x = (region[[2]] + region[[1]])/2
#      y = value[[numeric.column]]
#      circos.points(x, y, ...)
#  }

## ----eval=FALSE----------------------------------------------------------
#  circos.genomicLines(region, value, ...)
#  circos.genomicLines(region, value, numeric.column = c(1, 2))
#  circos.genomicLines(region, value, lwd, lty = "segment")
#  circos.genomicLines(region, value, area, baseline, border)
#  circos.genomicLines(region, value, sector.index, track.index)

## ----eval=FALSE----------------------------------------------------------
#  circos.genomicText(region, value, ...)
#  circos.genomicText(region, value, y, labels)
#  circos.genomicText(region, value, numeric.column, labels.column)
#  circos.genomicText(region, value, facing, niceFacing, adj)
#  circos.genomicText(region, value, sector.index, track.index)

## ----eval=FALSE----------------------------------------------------------
#  circos.genomicRect(region, value, ytop = 1, ybottom = 0)
#  circos.genomicRect(region, value, ytop.column = 2, ybottom = 0)
#  circos.genomicRect(region, value, col, border)

## ------------------------------------------------------------------------
col_fun = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
col_fun(c(-2, -1, -0.5, 0, 0.5, 1, 2))
col_fun = colorRamp2(breaks = -log10(c(1, 0.05, 1e-4)), 
    colors = c("green", "black", "red"))
p_value = c(0.8, 0.5, 0.001)
col_fun(-log10(p_value))

## ----eval=FALSE----------------------------------------------------------
#  bed = generateRandomBed(nc = 2)
#  # just note `numeric.column` is measured in `bed`
#  circos.genomicTrackPlotRegion(bed, numeric.column = 4,
#      panel.fun = function(region, value, ...) {
#          circos.genomicPoints(region, value, ...)
#          circos.genomicPoints(region, value)
#          # here `numeric.column` is measured in `value`
#          circos.genomicPoints(region, value, numeric.column = 1)
#  })

## ----echo=TRUE, eval=FALSE-----------------------------------------------
#  bedlist = list(generateRandomBed(), generateRandomBed())
#  circos.genomicTrackPlotRegion(bedlist,
#      panel.fun = function(region, value, ...) {
#          i = getI(...)
#          circos.genomicPoints(region, value, col = i, ...)
#  })
#  
#  # column 4 in the first bed and column 5 in the second bed
#  circos.genomicTrackPlotRegion(bedlist, numeric.column = c(4, 5),
#      panel.fun = function(region, value, ...) {
#          i = getI(...)
#          circos.genomicPoints(region, value, col = i, ...)
#  })

## ----eval=FALSE----------------------------------------------------------
#  bed = generateRandomBed(nc = 2)
#  circos.genomicTrackPlotRegion(bed, stack = TRUE,
#      panel.fun = function(region, value, ...) {
#          i = getI(...)
#          circos.genomicPoints(region, value, col = i, ...)
#  })

## ----echo=TRUE, eval=FALSE-----------------------------------------------
#  bedlist = list(generateRandomBed(), generateRandomBed())
#  circos.genomicTrackPlotRegion(bedlist, stack = TRUE,
#      panel.fun = function(region, value, ...) {
#          i = getI(...)
#          circos.genomicPoints(region, value, ...)
#  })

## ----genomic_application_points_0, eval = FALSE--------------------------
#  par(mar = c(1, 1, 1, 1))
#  set.seed(999)
#  
#  circos.par("track.height" = 0.1, start.degree = 90,
#      canvas.xlim = c(0, 1), canvas.ylim = c(0, 1), gap.degree = 270)
#  circos.initializeWithIdeogram(chromosome.index = "chr1", plotType = NULL)

## ----genomic_application_points_A, eval = FALSE, echo = 1:5--------------
#  ### track A
#  bed = generateRandomBed(nr = 300)
#  circos.genomicTrackPlotRegion(bed, panel.fun = function(region, value, ...) {
#      circos.genomicPoints(region, value, pch = 16, cex = 0.5, ...)
#  })
#  pos = get.cell.meta.data("yplot")
#  text(0, mean(pos), "A", adj = c(1.1, 0.5))

## ----genomic_application_points_B, eval = FALSE, echo = 1:10-------------
#  ### track B
#  bed = generateRandomBed(nr = 300)
#  circos.genomicTrackPlotRegion(bed, stack = TRUE,
#      panel.fun = function(region, value, ...) {
#          circos.genomicPoints(region, value, pch = 16, cex = 0.5, ...)
#  
#          i = getI(...)
#          cell.xlim = get.cell.meta.data("cell.xlim")
#          circos.lines(cell.xlim, c(i, i), lty = 2, col = "#00000040")
#  })
#  pos = get.cell.meta.data("yplot")
#  text(0, mean(pos), "B", adj = c(1.1, 0.5))

## ----genomic_application_points_C, eval = FALSE, echo = 1:10-------------
#  ### track C
#  bed1 = generateRandomBed(nr = 300)
#  bed2 = generateRandomBed(nr = 300)
#  bed_list = list(bed1, bed2)
#  circos.genomicTrackPlotRegion(bed_list,
#      panel.fun = function(region, value, ...) {
#          cex = (value[[1]] - min(value[[1]]))/(max(value[[1]]) - min(value[[1]]))
#          i = getI(...)
#          circos.genomicPoints(region, value, cex = cex, pch = 16, col = i, ...)
#  })
#  pos = get.cell.meta.data("yplot")
#  text(0, mean(pos), "C", adj = c(1.1, 0.5))

## ----genomic_application_points_D, eval = FALSE, echo = 1:9--------------
#  ### track D
#  circos.genomicTrackPlotRegion(bed_list, stack = TRUE,
#      panel.fun = function(region, value, ...) {
#          cex = (value[[1]] - min(value[[1]]))/(max(value[[1]]) - min(value[[1]]))
#          i = getI(...)
#          circos.genomicPoints(region, value, cex = cex, pch = 16, col = i, ...)
#          cell.xlim = get.cell.meta.data("cell.xlim")
#          circos.lines(cell.xlim, c(i, i), lty = 2, col = "#00000040")
#  })
#  pos = get.cell.meta.data("yplot")
#  text(0, mean(pos), "D", adj = c(1.1, 0.5))

## ----genomic_application_points_E, eval = FALSE, echo = 1:6--------------
#  ### track E
#  bed = generateRandomBed(nr = 300, nc = 4)
#  circos.genomicTrackPlotRegion(bed,
#      panel.fun = function(region, value, ...) {
#          circos.genomicPoints(region, value, cex = 0.5, pch = 16, col = 1:4, ...)
#  })
#  pos = get.cell.meta.data("yplot")
#  text(0, mean(pos), "E", adj = c(1.1, 0.5))

## ----genomic_application_points_F, eval = FALSE, echo = c(1:11, 14)------
#  ### track F
#  bed = generateRandomBed(nr = 300, nc = 4)
#  circos.genomicTrackPlotRegion(bed, stack = TRUE,
#      panel.fun = function(region, value, ...) {
#          cex = (value[[1]] - min(value[[1]]))/(max(value[[1]]) - min(value[[1]]))
#          i = getI(...)
#          circos.genomicPoints(region, value, cex = cex, pch = 16, col = i, ...)
#  
#          cell.xlim = get.cell.meta.data("cell.xlim")
#          circos.lines(cell.xlim, c(i, i), lty = 2, col = "#00000040")
#  })
#  pos = get.cell.meta.data("yplot")
#  text(0, mean(pos), "F", adj = c(1.1, 0.5))
#  circos.clear()

## ----genomic_application_points, echo = FALSE, out.width = "0.8\\textwidth", fig.cap = "Add points under different modes."----
par(mar = c(1, 1, 1, 1))
set.seed(999)

circos.par("track.height" = 0.1, start.degree = 90,
    canvas.xlim = c(0, 1), canvas.ylim = c(0, 1), gap.degree = 270)
circos.initializeWithIdeogram(chromosome.index = "chr1", plotType = NULL)
### track A
bed = generateRandomBed(nr = 300)
circos.genomicTrackPlotRegion(bed, panel.fun = function(region, value, ...) {
    circos.genomicPoints(region, value, pch = 16, cex = 0.5, ...)
})
pos = get.cell.meta.data("yplot")
text(0, mean(pos), "A", adj = c(1.1, 0.5))
### track B
bed = generateRandomBed(nr = 300)
circos.genomicTrackPlotRegion(bed, stack = TRUE, 
    panel.fun = function(region, value, ...) {
        circos.genomicPoints(region, value, pch = 16, cex = 0.5, ...)
        
        i = getI(...)
        cell.xlim = get.cell.meta.data("cell.xlim")
        circos.lines(cell.xlim, c(i, i), lty = 2, col = "#00000040")
})
pos = get.cell.meta.data("yplot")
text(0, mean(pos), "B", adj = c(1.1, 0.5))
### track C
bed1 = generateRandomBed(nr = 300)
bed2 = generateRandomBed(nr = 300)
bed_list = list(bed1, bed2)
circos.genomicTrackPlotRegion(bed_list, 
    panel.fun = function(region, value, ...) {
        cex = (value[[1]] - min(value[[1]]))/(max(value[[1]]) - min(value[[1]]))
        i = getI(...)
        circos.genomicPoints(region, value, cex = cex, pch = 16, col = i, ...)
})
pos = get.cell.meta.data("yplot")
text(0, mean(pos), "C", adj = c(1.1, 0.5))
### track D
circos.genomicTrackPlotRegion(bed_list, stack = TRUE, 
    panel.fun = function(region, value, ...) {
        cex = (value[[1]] - min(value[[1]]))/(max(value[[1]]) - min(value[[1]]))
        i = getI(...)
        circos.genomicPoints(region, value, cex = cex, pch = 16, col = i, ...)
        cell.xlim = get.cell.meta.data("cell.xlim")
        circos.lines(cell.xlim, c(i, i), lty = 2, col = "#00000040")
})
pos = get.cell.meta.data("yplot")
text(0, mean(pos), "D", adj = c(1.1, 0.5))
### track E
bed = generateRandomBed(nr = 300, nc = 4)
circos.genomicTrackPlotRegion(bed, 
    panel.fun = function(region, value, ...) {
        circos.genomicPoints(region, value, cex = 0.5, pch = 16, col = 1:4, ...)
})
pos = get.cell.meta.data("yplot")
text(0, mean(pos), "E", adj = c(1.1, 0.5))
### track F
bed = generateRandomBed(nr = 300, nc = 4)
circos.genomicTrackPlotRegion(bed, stack = TRUE, 
    panel.fun = function(region, value, ...) {
        cex = (value[[1]] - min(value[[1]]))/(max(value[[1]]) - min(value[[1]]))
        i = getI(...)
        circos.genomicPoints(region, value, cex = cex, pch = 16, col = i, ...)
        
        cell.xlim = get.cell.meta.data("cell.xlim")
        circos.lines(cell.xlim, c(i, i), lty = 2, col = "#00000040")
})
pos = get.cell.meta.data("yplot")
text(0, mean(pos), "F", adj = c(1.1, 0.5))
circos.clear()

## ----genomic_application_lines_0, eval = FALSE---------------------------
#  par(mar = c(1, 1, 1, 1))
#  circos.par("track.height" = 0.1, start.degree = 90,
#      canvas.xlim = c(0, 1), canvas.ylim = c(0, 1), gap.degree = 270)
#  circos.initializeWithIdeogram(chromosome.index = "chr1", plotType = NULL)

## ----genomic_application_lines_A, eval = FALSE, echo = 1:6---------------
#  ### track A
#  bed = generateRandomBed(nr = 500)
#  circos.genomicTrackPlotRegion(bed,
#      panel.fun = function(region, value, ...) {
#          circos.genomicLines(region, value, type = "l", ...)
#  })
#  pos = get.cell.meta.data("yplot")
#  text(0, mean(pos), "A", adj = c(1.1, 0.5))

## ----genomic_application_lines_B, eval = FALSE, echo = 1:9---------------
#  ### track B
#  bed1 = generateRandomBed(nr = 500)
#  bed2 = generateRandomBed(nr = 500)
#  bed_list = list(bed1, bed2)
#  circos.genomicTrackPlotRegion(bed_list,
#      panel.fun = function(region, value, ...) {
#          i = getI(...)
#          circos.genomicLines(region, value, col = i, ...)
#  })
#  pos = get.cell.meta.data("yplot")
#  text(0, mean(pos), "B", adj = c(1.1, 0.5))

## ----genomic_application_lines_C, eval = FALSE, echo = 1:6---------------
#  ### track C
#  circos.genomicTrackPlotRegion(bed_list, stack = TRUE,
#      panel.fun = function(region, value, ...) {
#          i = getI(...)
#          circos.genomicLines(region, value, col = i, ...)
#  })
#  pos = get.cell.meta.data("yplot")
#  text(0, mean(pos), "C", adj = c(1.1, 0.5))

## ----genomic_application_lines_D, eval = FALSE, echo = 1:6---------------
#  ### track D
#  bed = generateRandomBed(nr = 500, nc = 4)
#  circos.genomicTrackPlotRegion(bed,
#      panel.fun = function(region, value, ...) {
#          circos.genomicLines(region, value, col = 1:4, ...)
#  })
#  pos = get.cell.meta.data("yplot")
#  text(0, mean(pos), "D", adj = c(1.1, 0.5))

## ----genomic_application_lines_E, eval = FALSE, echo = 1:7---------------
#  ### track E
#  bed = generateRandomBed(nr = 500, nc = 4)
#  circos.genomicTrackPlotRegion(bed, stack = TRUE,
#      panel.fun = function(region, value, ...) {
#          i = getI(...)
#          circos.genomicLines(region, value, col = i, ...)
#  })
#  pos = get.cell.meta.data("yplot")
#  text(0, mean(pos), "E", adj = c(1.1, 0.5))

## ----genomic_application_lines_F, eval = FALSE, echo = c(1:7, 10)--------
#  ### track F
#  bed = generateRandomBed(nr = 200)
#  circos.genomicTrackPlotRegion(bed,
#      panel.fun = function(region, value, ...) {
#          circos.genomicLines(region, value, type = "segment", lwd = 2,
#              col = rand_color(nrow(region)), ...)
#  })
#  pos = get.cell.meta.data("yplot")
#  text(0, mean(pos), "F", adj = c(1.1, 0.5))
#  circos.clear()

## ----genomic_application_lines, echo = FALSE, out.width = "0.8\\textwidth", fig.cap = "Add lines under different modes."----
par(mar = c(1, 1, 1, 1))
circos.par("track.height" = 0.1, start.degree = 90,
    canvas.xlim = c(0, 1), canvas.ylim = c(0, 1), gap.degree = 270)
circos.initializeWithIdeogram(chromosome.index = "chr1", plotType = NULL)
### track A
bed = generateRandomBed(nr = 500)
circos.genomicTrackPlotRegion(bed, 
    panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, type = "l", ...)
})
pos = get.cell.meta.data("yplot")
text(0, mean(pos), "A", adj = c(1.1, 0.5))
### track B
bed1 = generateRandomBed(nr = 500)
bed2 = generateRandomBed(nr = 500)
bed_list = list(bed1, bed2)
circos.genomicTrackPlotRegion(bed_list, 
    panel.fun = function(region, value, ...) {
        i = getI(...)
        circos.genomicLines(region, value, col = i, ...)
})
pos = get.cell.meta.data("yplot")
text(0, mean(pos), "B", adj = c(1.1, 0.5))
### track C
circos.genomicTrackPlotRegion(bed_list, stack = TRUE, 
    panel.fun = function(region, value, ...) {
        i = getI(...)
        circos.genomicLines(region, value, col = i, ...)
})
pos = get.cell.meta.data("yplot")
text(0, mean(pos), "C", adj = c(1.1, 0.5))
### track D
bed = generateRandomBed(nr = 500, nc = 4)
circos.genomicTrackPlotRegion(bed, 
    panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, col = 1:4, ...)
})
pos = get.cell.meta.data("yplot")
text(0, mean(pos), "D", adj = c(1.1, 0.5))
### track E
bed = generateRandomBed(nr = 500, nc = 4)
circos.genomicTrackPlotRegion(bed, stack = TRUE, 
    panel.fun = function(region, value, ...) {
        i = getI(...)
        circos.genomicLines(region, value, col = i, ...)
})
pos = get.cell.meta.data("yplot")
text(0, mean(pos), "E", adj = c(1.1, 0.5))
### track F
bed = generateRandomBed(nr = 200)
circos.genomicTrackPlotRegion(bed, 
    panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, type = "segment", lwd = 2, 
            col = rand_color(nrow(region)), ...)
})
pos = get.cell.meta.data("yplot")
text(0, mean(pos), "F", adj = c(1.1, 0.5))
circos.clear()

## ----genomic_application_rect_0, eval = FALSE----------------------------
#  par(mar = c(1, 1, 1, 1))
#  circos.par("track.height" = 0.1, start.degree = 90,
#      canvas.xlim = c(0, 1), canvas.ylim = c(0, 1), gap.degree = 270)
#  circos.initializeWithIdeogram(chromosome.index = "chr1", plotType = NULL)
#  f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))

## ----genomic_application_rect_A, eval = FALSE, echo = 1:6----------------
#  ### track A
#  bed = generateRandomBed(nr = 100, nc = 4)
#  circos.genomicTrackPlotRegion(bed, stack = TRUE,
#      panel.fun = function(region, value, ...) {
#          circos.genomicRect(region, value, col = f(value[[1]]), border = NA, ...)
#  })
#  pos = get.cell.meta.data("yplot")
#  text(0, mean(pos), "A", adj = c(1.1, 0.5))

## ----genomic_application_rect_B, eval = FALSE, echo = 1:10---------------
#  ### track B
#  bed1 = generateRandomBed(nr = 100)
#  bed2 = generateRandomBed(nr = 100)
#  bed_list = list(bed1, bed2)
#  circos.genomicTrackPlotRegion(bed_list, stack = TRUE,
#      panel.fun = function(region, value, ...) {
#          i = getI(...)
#          circos.genomicRect(region, value, ytop = i + 0.4, ybottom = i - 0.4,
#              col = f(value[[1]]), ...)
#  })
#  pos = get.cell.meta.data("yplot")
#  text(0, mean(pos), "B", adj = c(1.1, 0.5))

## ----genomic_application_rect_C, eval = FALSE, echo = 1:7----------------
#  ### track C
#  circos.genomicTrackPlotRegion(bed_list, ylim = c(0, 3),
#      panel.fun = function(region, value, ...) {
#          i = getI(...)
#          circos.genomicRect(region, value, ytop = i + 0.4, ybottom = i - 0.4,
#              col = f(value[[1]]), ...)
#  })
#  pos = get.cell.meta.data("yplot")
#  text(0, mean(pos), "C", adj = c(1.1, 0.5))

## ----genomic_application_rect_D, eval = FALSE, echo = c(1:10, 13)--------
#  ### track D
#  bed = generateRandomBed(nr = 200)
#  circos.genomicTrackPlotRegion(bed,
#      panel.fun = function(region, value, ...) {
#          circos.genomicRect(region, value, ytop.column = 1, ybottom = 0,
#              col = ifelse(value[[1]] > 0, "red", "green"), ...)
#  
#          cell.xlim = get.cell.meta.data("cell.xlim")
#          circos.lines(cell.xlim, c(0, 0), lty = 2, col = "#00000040")
#  })
#  pos = get.cell.meta.data("yplot")
#  text(0, mean(pos), "D", adj = c(1.1, 0.5))
#  circos.clear()

## ----genomic_application_rect, echo = FALSE, out.width = "0.8\\textwidth", fig.cap = "Add rectangles under different modes."----
par(mar = c(1, 1, 1, 1))
circos.par("track.height" = 0.1, start.degree = 90,
    canvas.xlim = c(0, 1), canvas.ylim = c(0, 1), gap.degree = 270)
circos.initializeWithIdeogram(chromosome.index = "chr1", plotType = NULL)
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
### track A
bed = generateRandomBed(nr = 100, nc = 4)
circos.genomicTrackPlotRegion(bed, stack = TRUE, 
    panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = f(value[[1]]), border = NA, ...)
})
pos = get.cell.meta.data("yplot")
text(0, mean(pos), "A", adj = c(1.1, 0.5))
### track B
bed1 = generateRandomBed(nr = 100)
bed2 = generateRandomBed(nr = 100)
bed_list = list(bed1, bed2)
circos.genomicTrackPlotRegion(bed_list, stack = TRUE, 
    panel.fun = function(region, value, ...) {
        i = getI(...)
        circos.genomicRect(region, value, ytop = i + 0.4, ybottom = i - 0.4,
            col = f(value[[1]]), ...)
})
pos = get.cell.meta.data("yplot")
text(0, mean(pos), "B", adj = c(1.1, 0.5))
### track C
circos.genomicTrackPlotRegion(bed_list, ylim = c(0, 3), 
    panel.fun = function(region, value, ...) {
        i = getI(...)
        circos.genomicRect(region, value, ytop = i + 0.4, ybottom = i - 0.4, 
            col = f(value[[1]]), ...)
})
pos = get.cell.meta.data("yplot")
text(0, mean(pos), "C", adj = c(1.1, 0.5))
### track D
bed = generateRandomBed(nr = 200)
circos.genomicTrackPlotRegion(bed, 
    panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
            col = ifelse(value[[1]] > 0, "red", "green"), ...)
            
        cell.xlim = get.cell.meta.data("cell.xlim")
        circos.lines(cell.xlim, c(0, 0), lty = 2, col = "#00000040")
})
pos = get.cell.meta.data("yplot")
text(0, mean(pos), "D", adj = c(1.1, 0.5))
circos.clear()

## ----eval=FALSE----------------------------------------------------------
#  circos.genomicTrackPlotRegion(bed, ylim = c(-1, 1),
#      panel.fun = function(region, value, ...) {
#          circos.genomicPoints(region, value, ...)
#  
#          cell.xlim = get.cell.meta.data("cell.xlim")
#          for(h in c(-1, -0.5, 0, 0.5, 1)) {
#              circos.lines(cell.xlim, c(0, 0), lty = 2, col = "grey")
#          }
#          circos.text(x, y, labels)
#          circos.axis("top")
#  })

## ----genomic_links, echo = 2:11, out.width = "\\textwidth", out.height = "0.5\\textwidth", fig.width = 12, fig.height = 6, fig.cap = "Add links from two sets of genomic regions."----
par(mfrow = c(1, 2), mar = c(1, 1, 1, 1))
bed1 = generateRandomBed(nr = 100)
bed1 = bed1[sample(nrow(bed1), 20), ]
bed2 = generateRandomBed(nr = 100)
bed2 = bed2[sample(nrow(bed2), 20), ]
circos.initializeWithIdeogram(plotType = c("axis", "labels"))
circos.genomicLink(bed1, bed2)
circos.clear()

circos.initializeWithIdeogram(plotType = c("axis", "labels"))
circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.5), 
    border = NA)
circos.clear()
par(mfrow = c(1, 1))

## ----genomic_highlight_1, eval = FALSE-----------------------------------
#  circos.par("track.height" = 0.1)
#  circos.initializeWithIdeogram(plotType = c("axis", "labels"))
#  
#  for(i in 1:5) {
#      bed = generateRandomBed(nr = 100)
#      circos.genomicTrackPlotRegion(bed, panel.fun = function(region, value, ...) {
#          circos.genomicPoints(region, value, pch = 16, cex = 0.5, ...)
#      })
#  }

## ----genomic_highlight_2, eval = FALSE-----------------------------------
#  highlight.chromosome(c("chrX", "chrY", "chr1"))
#  highlight.chromosome("chr3", col = "#00FF0040", padding = c(0.05, 0.05, 0.15, 0.05))
#  highlight.chromosome("chr5", col = NA, border = "red", lwd = 2,
#      padding = c(0.05, 0.05, 0.15, 0.05))
#  highlight.chromosome("chr7", col = "#0000FF40", track.index = c(2, 4, 5))
#  highlight.chromosome(c("chr9", "chr10", "chr11"), col = NA, border = "green",
#      lwd = 2, track.index = c(2, 4, 5))
#  highlight.chromosome(paste0("chr", c(1:22, "X", "Y")), col = "#FFFF0040",
#      track.index = 6)
#  circos.clear()

## ----genomic_highlight, echo = FALSE, out.width = "0.8\\textwidth", fig.cap = "Highlight different chromosomes and tracks."----
par(mar = c(1, 1, 1, 1))
circos.par("track.height" = 0.1)
circos.initializeWithIdeogram(plotType = c("axis", "labels"))

for(i in 1:5) {
    bed = generateRandomBed(nr = 100)
    circos.genomicTrackPlotRegion(bed, panel.fun = function(region, value, ...) {
        circos.genomicPoints(region, value, pch = 16, cex = 0.5, ...)
    })
}
highlight.chromosome(c("chrX", "chrY", "chr1"))
highlight.chromosome("chr3", col = "#00FF0040", padding = c(0.05, 0.05, 0.15, 0.05))
highlight.chromosome("chr5", col = NA, border = "red", lwd = 2, 
    padding = c(0.05, 0.05, 0.15, 0.05))
highlight.chromosome("chr7", col = "#0000FF40", track.index = c(2, 4, 5))
highlight.chromosome(c("chr9", "chr10", "chr11"), col = NA, border = "green", 
    lwd = 2, track.index = c(2, 4, 5))
highlight.chromosome(paste0("chr", c(1:22, "X", "Y")), col = "#FFFF0040", 
    track.index = 6)
circos.clear()

## ----eval=FALSE----------------------------------------------------------
#  circos.genomicTrackPlotRegion(data, panel.fun = function(region, value, ...) {
#      circos.genomicPoints(region, value, posTransform = posTransform.default, ...)
#  })

## ----eval=FALSE----------------------------------------------------------
#  circos.genomicPosTransformLines(data, posTransform = posTransform.default)
#  circos.genomicPosTransformLines(data, posTransform = posTransform.default,
#      horizontalLine = "top")
#  circos.genomicPosTransformLines(data, posTransform = posTransform.default,
#      direction = "outside")

## ----genomic_postransform_1, eval = FALSE--------------------------------
#  circos.par(cell.padding = c(0, 0, 0, 0))
#  circos.initializeWithIdeogram()
#  bed = generateRandomBed(nr = 100, nc = 4)
#  
#  # note how 'horizontalLine' works
#  circos.genomicPosTransformLines(bed, posTransform = posTransform.default,
#      horizontalLine = "top", track.height = 0.1)
#  
#  f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
#  circos.genomicTrackPlotRegion(bed, stack = TRUE,
#      panel.fun = function(region, value, ...) {
#          circos.genomicRect(region, value, col = f(value[[1]]),
#              border = f(value[[1]]), posTransform = posTransform.default, ...)
#  }, bg.border = NA)
#  
#  circos.clear()

## ----genomic_postransform_2, eval = FALSE--------------------------------
#  circos.par(cell.padding = c(0, 0, 0, 0))
#  circos.initializeWithIdeogram(plotType = NULL)
#  
#  circos.genomicTrackPlotRegion(bed, stack = TRUE,
#      panel.fun = function(region, value, ...) {
#          circos.genomicRect(region, value, col = f(value[[1]]),
#              border = f(value[[1]]), posTransform = posTransform.default, ...)
#  }, bg.border = NA)
#  
#  circos.genomicPosTransformLines(bed, posTransform = posTransform.default,
#      direction = "outside", horizontalLine = "bottom", track.height = 0.1)
#  
#  cytoband = read.cytoband()$df
#  circos.genomicTrackPlotRegion(cytoband, stack = TRUE,
#      panel.fun = function(region, value, ...) {
#          circos.genomicRect(region, value, col = cytoband.col(value[[2]]), border = NA)
#          xlim = get.cell.meta.data("xlim")
#          ylim = get.cell.meta.data("ylim")
#          circos.rect(xlim[1], ylim[1], xlim[2], ylim[2], border = "black")
#  }, track.height = 0.05, bg.border = NA)
#  
#  circos.clear()

## ----genomic_postransform, echo = FALSE, out.width = "0.6\\textwidth", out.height = "1.2\\textwidth", fig.width = 6, fig.height = 12, fig.cap = "Position transformation with {\\tt posTransform.default}. A) transformation is inside; B) transformation is outside."----
par(mfrow = c(2, 1), mar = c(1, 1, 1, 1))
circos.par(cell.padding = c(0, 0, 0, 0))
circos.initializeWithIdeogram()
bed = generateRandomBed(nr = 100, nc = 4)

# note how 'horizontalLine' works
circos.genomicPosTransformLines(bed, posTransform = posTransform.default, 
    horizontalLine = "top", track.height = 0.1)

f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrackPlotRegion(bed, stack = TRUE, 
    panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = f(value[[1]]), 
            border = f(value[[1]]), posTransform = posTransform.default, ...)
}, bg.border = NA)

circos.clear()
text(-0.9, 0.9, "A", cex = 1.5)
circos.par(cell.padding = c(0, 0, 0, 0))
circos.initializeWithIdeogram(plotType = NULL)

circos.genomicTrackPlotRegion(bed, stack = TRUE, 
    panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = f(value[[1]]), 
            border = f(value[[1]]), posTransform = posTransform.default, ...)
}, bg.border = NA)

circos.genomicPosTransformLines(bed, posTransform = posTransform.default, 
    direction = "outside", horizontalLine = "bottom", track.height = 0.1)

cytoband = read.cytoband()$df
circos.genomicTrackPlotRegion(cytoband, stack = TRUE, 
    panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = cytoband.col(value[[2]]), border = NA)
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        circos.rect(xlim[1], ylim[1], xlim[2], ylim[2], border = "black")
}, track.height = 0.05, bg.border = NA)

circos.clear()
text(-0.9, 0.9, "B", cex = 1.5)
par(mfrow = c(1, 1))

## ----eval=FALSE----------------------------------------------------------
#  bed = generateRandomBed(nr = 400, fun = function(k) rep("text", k))
#  circos.genomicTrackPlotRegion(bed, ylim = c(0, 1),
#      panel.fun = function(region, value, ...) {
#          circos.genomicText(region, value, y = 0, labels.column = 1,
#              facing = "clockwise", adj = c(0, 0.5),
#              posTransform = posTransform.text, cex = 0.8)
#  }, track.height = 0.1, bg.border = NA)

## ----eval=FALSE----------------------------------------------------------
#  i_track = get.cell.meta.data("track.index")  # the nearest track
#  # we put `y`, `labels`, ... into a self-defined function
#  # because these parameters will affect the text position
#  circos.genomicPosTransformLines(bed, direction = "outside",
#      posTransform = function(region, value)
#          posTransform.text(region, y = 0, labels = value[[1]],
#              cex = 0.8, track.index = i_track)
#  )

## ----eval=FALSE----------------------------------------------------------
#  circos.genomicTrackPlotRegion(bed, ylim = c(0, 1), track.height = 0.1, bg.border = NA)
#  i_track = get.cell.meta.data("track.index")  # remember this empty track, we'll come back
#  
#  circos.genomicTrackPlotRegion(bed, ylim = c(0, 1),
#      panel.fun = function(region, value, ...) {
#          circos.genomicText(region, value, y = 1, labels.column = 1,
#              facing = "clockwise", adj = c(1, 0.5),
#              posTransform = posTransform.text, cex = 0.8)
#  }, track.height = 0.1, bg.border = NA)
#  tr_track = get.cell.meta.data("track.index") # position transformation track
#  
#  # because `circos.genomicPosTransformLines` is implemented by
#  # `circos.trackPlotRegion`, it accepts `track.index` argument.
#  circos.genomicPosTransformLines(bed,
#      posTransform = function(region, value)
#          posTransform.text(region, y = 1, labels = value[[1]],
#              cex = 0.8, track.index = tr_track),
#      direction = "inside", track.index = i_track
#  )

## ----eval=FALSE----------------------------------------------------------
#  circos.genomicTrackPlotRegion(bed, ylim = c(0, 1),
#      panel.fun = function(region, value, ...) {
#          circos.genomicText(region, value, y = 0, labels.column = 1,
#              facing = "clockwise", adj = c(0, 0.5), posTransform = posTransform.text,
#              cex = 0.8, padding = 0.2)
#  }, track.height = 0.1, bg.border = NA)
#  
#  i_track = get.cell.meta.data("track.index")  # previous track
#  circos.genomicPosTransformLines(bed,
#      posTransform = function(region, value) posTransform.text(region, y = 0,
#          labels = value[[1]], cex = 0.8, padding = 0.2, track.index = i_track),
#      direction = "outside"
#  )

## ----genomic_text_pos_transformation, echo = FALSE, out.width = "\\textwidth", fig.cap = "Transformation of text positions."----
source("src/genomic-07-posTransformLinesText.R")

## ----genome_more_labels_1, eval = FALSE----------------------------------
#  set.seed(999)
#  bed = generateRandomBed(nr = 800, fun = function(k) rep("text", k))
#  
#  par(mar = c(1, 1, 1, 1))
#  circos.par(start.degree = 75, canvas.xlim = c(0, 1), canvas.ylim = c(0, 1),
#      gap.degree = 300, cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0))
#  circos.initializeWithIdeogram(plotType = NULL, chromosome.index = "chr1")

## ----genome_more_labels_2, eval = FALSE----------------------------------
#  circos.genomicTrackPlotRegion(bed, ylim = c(0, 1),
#      panel.fun = function(region, value, ...) {
#  
#      # original positions
#      circos.genomicPoints(region, data.frame(rep(0, nrow(region))), pch = 16)
#  
#      xlim = get.cell.meta.data("xlim")
#      breaks = seq(xlim[1], xlim[2], length.out = 7)
#      midpoints = (region[[1]] + region[[2]])/2
#      for(i in seq_along(breaks)[-1]) {
#          # index for current interval
#          l = midpoints >= breaks[i - 1] & midpoints < breaks[i]
#          # if there is no data in this interval
#          if(sum(l) == 0) next
#  
#          # sub-regions in this interval
#          sub_region = region[l, , drop = FALSE]
#          sub_value = value[l, , drop = FALSE]
#  
#          # note here i == 2, 4, 6, ... corresponds to odd intervals
#          # text in odd intervals are in lower position
#          if(i %% 2 == 0) {
#              y = 0.2
#          } else {
#              y = 0.8
#          }
#  
#          # get the transformed position and add text with new positions
#          tr_region = posTransform.text(sub_region, y = y, labels = sub_value[[1]],
#              cex = 0.8, adj = c(0, 0.5))
#          circos.genomicText(tr_region, sub_value, labels.column = 1, y = y,
#              adj = c(0, 0.5), facing = "clockwise", niceFacing = TRUE, cex = 0.8)
#  
#          # add position transformation lines for odd intervals
#          if(i %% 2 == 0) {
#              for(i in seq_len(nrow(sub_region))) {
#                  x = c( (sub_region[i, 1] + sub_region[i, 2])/2,
#                         (sub_region[i, 1] + sub_region[i, 2])/2,
#                         (tr_region[i, 1] + tr_region[i, 2])/2,
#                         (tr_region[i, 1] + tr_region[i, 2])/2)
#                  y = c(0, 0.2/3, 0.2/3*2, 0.2)
#                  circos.lines(x, y)
#              }
#          } else { # add position transformation lines for even intervals
#              median_sub_region_midpoint = median(midpoints[l])
#              sub_region_width = max(midpoints[l]) - min(midpoints[l])
#              for(i in seq_len(nrow(sub_region))) {
#                  x = c( (sub_region[i, 1] + sub_region[i, 2])/2,
#                         (sub_region[i, 1] + sub_region[i, 2])/2,
#                         median_sub_region_midpoint +
#                             sub_region_width*(i - nrow(sub_region))/nrow(sub_region) * 0.2,
#                         median_sub_region_midpoint +
#                             sub_region_width*(i - nrow(sub_region))/nrow(sub_region) * 0.2,
#                         (tr_region[i, 1] + tr_region[i, 2])/2,
#                         (tr_region[i, 1] + tr_region[i, 2])/2)
#                  y = c(0, 0.1, 0.2, 0.6, 0.7, 0.8)
#                  circos.lines(x, y)
#              }
#          }
#      }
#  
#  }, track.height = 0.2, bg.border = NA)
#  
#  circos.clear()

## ----genome_more_labels, echo = FALSE, out.width = "\\textwidth", out.height = "0.5\\textwidth", fig.width = 14, fig.height = 7, fig.cap = "Position transformation for a lot of text. A) put text on two layers; B) put text on one layer."----
par(mfrow = c(1, 2), mar = c(1, 1, 1, 1))
set.seed(999)
bed = generateRandomBed(nr = 800, fun = function(k) rep("text", k))

par(mar = c(1, 1, 1, 1))
circos.par(start.degree = 75, canvas.xlim = c(0, 1), canvas.ylim = c(0, 1), 
    gap.degree = 300, cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0))
circos.initializeWithIdeogram(plotType = NULL, chromosome.index = "chr1")
text(0.1, 0.1, "A", cex = 1.5)

circos.genomicTrackPlotRegion(bed, ylim = c(0, 1), 
    panel.fun = function(region, value, ...) {
    
    # original positions
    circos.genomicPoints(region, data.frame(rep(0, nrow(region))), pch = 16)
        
    xlim = get.cell.meta.data("xlim")
    breaks = seq(xlim[1], xlim[2], length.out = 7)
    midpoints = (region[[1]] + region[[2]])/2
    for(i in seq_along(breaks)[-1]) {
        # index for current interval
        l = midpoints >= breaks[i - 1] & midpoints < breaks[i]
        # if there is no data in this interval
        if(sum(l) == 0) next
        
        # sub-regions in this interval
        sub_region = region[l, , drop = FALSE]
        sub_value = value[l, , drop = FALSE]
        
        # note here i == 2, 4, 6, ... corresponds to odd intervals
        # text in odd intervals are in lower position
        if(i %% 2 == 0) {
            y = 0.2
        } else {
            y = 0.8
        }
        
        # get the transformed position and add text with new positions
        tr_region = posTransform.text(sub_region, y = y, labels = sub_value[[1]], 
            cex = 0.8, adj = c(0, 0.5))
        circos.genomicText(tr_region, sub_value, labels.column = 1, y = y, 
            adj = c(0, 0.5), facing = "clockwise", niceFacing = TRUE, cex = 0.8)
        
        # add position transformation lines for odd intervals
        if(i %% 2 == 0) {    
            for(i in seq_len(nrow(sub_region))) {
                x = c( (sub_region[i, 1] + sub_region[i, 2])/2,
                       (sub_region[i, 1] + sub_region[i, 2])/2,
                       (tr_region[i, 1] + tr_region[i, 2])/2,
                       (tr_region[i, 1] + tr_region[i, 2])/2)
                y = c(0, 0.2/3, 0.2/3*2, 0.2)
                circos.lines(x, y)
            }
        } else { # add position transformation lines for even intervals
            median_sub_region_midpoint = median(midpoints[l])
            sub_region_width = max(midpoints[l]) - min(midpoints[l])
            for(i in seq_len(nrow(sub_region))) {
                x = c( (sub_region[i, 1] + sub_region[i, 2])/2,
                       (sub_region[i, 1] + sub_region[i, 2])/2,
                       median_sub_region_midpoint + 
                           sub_region_width*(i - nrow(sub_region))/nrow(sub_region) * 0.2,
                       median_sub_region_midpoint + 
                           sub_region_width*(i - nrow(sub_region))/nrow(sub_region) * 0.2,
                       (tr_region[i, 1] + tr_region[i, 2])/2,
                       (tr_region[i, 1] + tr_region[i, 2])/2)
                y = c(0, 0.1, 0.2, 0.6, 0.7, 0.8)
                circos.lines(x, y)
            }
        }
    }
    
}, track.height = 0.2, bg.border = NA)

circos.clear()
par(mar = c(1, 1, 1, 1))
circos.par(start.degree = 75, canvas.xlim = c(0, 1), canvas.ylim = c(0, 1), gap.degree = 300, cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0))
circos.initializeWithIdeogram(plotType = NULL, chromosome.index = "chr1")
circos.genomicTrackPlotRegion(bed, ylim = c(0, 1), panel.fun = function(region, value, ...) {
    circos.genomicText(region, value, y = 0, labels.column = 1, facing = "clockwise", adj = c(0, 0.5),
        posTransform = posTransform.text, cex = 0.8, niceFacing = F)
}, track.height = 0.1, bg.border = NA)
i_track = get.cell.meta.data("track.index")

circos.genomicPosTransformLines(bed, 
    posTransform = function(region, value) posTransform.text(region, y = 0, labels = value[[1]], cex = 0.8, track.index = i_track),
    direction = "outside"
)

circos.genomicTrackPlotRegion(bed, ylim = c(0, 1), panel.fun = function(region, value, ...) {
    circos.points( (region[[1]] + region[[2]])/2, rep(0.5, nrow(region)), pch = 16)
}, track.height = 0.02, bg.border = NA)

circos.clear()

text(0.1, 0.1, "B", cex = 1.5)

## ----eval = FALSE--------------------------------------------------------
#  circos.genomicDensity(bed)
#  circos.genomicDensity(bed, baseline = 0)
#  circos.genomicDensity(bed, window.size = 1e6)
#  circos.genomicDensity(bedlist, col = c("#FF000080", "#0000FF80"))

## ----eval = FALSE--------------------------------------------------------
#  circos.genoimcRainfall(bed)
#  circos.genoimcRainfall(bedlist, col = c("red", "green"))

## ----genomic_rainfall, out.width = "0.8\\textwidth", fig.cap = "Genomic rainfall plot and densities."----
load(paste0(system.file(package = "circlize"), "/extdata/DMR.RData"))
par(mar = c(1, 1, 1, 1))
circos.initializeWithIdeogram(plotType = c("axis", "labels"))

bed_list = list(DMR_hyper, DMR_hypo)
circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = c("#FF000080", "#0000FF80"))
circos.genomicDensity(DMR_hyper, col = c("#FF000080"), track.height = 0.1)
circos.genomicDensity(DMR_hypo, col = c("#0000FF80"), track.height = 0.1)
circos.clear()

