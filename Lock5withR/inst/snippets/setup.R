require(grDevices)
require(datasets)
require(stats)
require(lattice)
require(grid)
require(mosaic)
require(vcd)
require(fastR)
require(Lock5withR)
trellis.par.set(theme=col.mosaic())
trellis.par.set(fontsize=list(text=9))
options(keep.blank.line=FALSE)
options(width=90)
options(digits=3)
require(knitr)
opts_chunk$set(
tidy=TRUE,
boxedLabel=TRUE,
size='small',
dev='pdf',
fig.width=3, fig.height=2,
out.width=".47\\textwidth",
fig.align='center',
fig.show='hold',
keep.source=TRUE,
comment=NA
)
opts_template$set(fig4 = list(fig.height = 4, fig.width = 6, out.width=".65\\textwidth"))
opts_template$set(fig3 = list(fig.height = 5*.35, fig.width = 8*.35, out.width=".31\\textwidth"))
opts_template$set(fig2 = list(fig.height = 2, fig.width = 3, out.width=".47\\textwidth"))
opts_template$set(fig1 = list(fig.height = 2, fig.width = 8, out.width=".95\\textwidth"))
opts_template$set(figbig = list(fig.height = 9, fig.width = 12, out.width=".95\\textwidth"))
knit_hooks$set(seed = function(before, options, envir) {
    if (before) set.seed(options$seed) 
})

knit_hooks$set(chunk = function (x, options) {
  if ( !is.null(options$boxedLabel) && options$boxedLabel && 
         ! grepl("unnamed-chunk", options$label) &&
		(is.null(options$echo) || options$echo) ) {
		labeling <- paste0( 
			"\\endgraf\\nobreak\\null\\endgraf\\penalty-2\\kern-.5\\baselineskip",
			"\n\n",
			"\\hfill \\makebox[0pt][r]{\\fbox{\\tiny ",
			options$label,
			"}}", 
			"\\endgraf",
			"\\kern-4.5ex\n\n")
	}  else {
		labeling <- ""
	}
    ai = knitr:::output_asis(x, options)
    col = if (!ai)
        paste(knitr:::color_def(options$background), if (!knitr:::is_tikz_dev(options))
            "\\color{fgcolor}", sep = "")
    k1 = paste(col, "\\begin{kframe}\n", sep = "")
    k2 = "\\end{kframe}"
    x = knitr:::.rm.empty.envir(paste(k1, labeling, x, k2, sep = ""))
    size = if (options$size == "normalsize")
        ""
    else sprintf("\\%s", options$size)
    if (!ai)
        x = sprintf("\\begin{knitrout}%s\n%s\n\\end{knitrout}",
            size, x)
    if (options$split) {
        name = knitr:::fig_path(".tex", options)
        if (!file.exists(dirname(name)))
            dir.create(dirname(name))
        cat(x, file = name)
        sprintf("\\input{%s}", name)
}
else x 
}
)

knit_hooks$set(document = function(x) { 
			   sub('\\usepackage[]{color}', '\\usepackage[]{xcolor}', 
			   x, fixed = TRUE) 
}) 

