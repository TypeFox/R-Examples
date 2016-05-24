#' Density and histogram plots from PLOS Article Level Metrics data
#'
#' @import ggplot2
#' @importFrom grid unit grid.text gpar viewport grid.newpage
#' @importFrom reshape2 melt
#' @importFrom plyr round_any
#' @importFrom lubridate year
#' @export
#'
#' @param input A data.frame, usuaally from a call to \code{link{alm}}.
#' @param source Data source (column) to plot. Can be a single element, or a
#'    character vector.
#' @param color Color of the density plot and the title. Can be a hex color or
#'    rgb, etc.
#' @param title Title for the plot, in top matching the color of the density plot.
#' @param description Optional description, subtending the title.
#' @param plot_type Type of plot, one of density (default) or histogram.
#' @author Martin Fenner, mfenner@@plos.org, modified by Scott Chamberlain
#' @references See a tutorial/vignette for alm at
#' \url{http://ropensci.org/tutorials/alm_tutorial.html}
#' @examples \dontrun{
#' library('rplos'); library('plyr')
#' dois <- searchplos(q='*:*', fl="id",
#'    fq=list(
#'      'cross_published_journal_key:PLoSONE',
#'      'doc_type:full',
#'      'publication_date:[2010-01-01T00:00:00Z TO 2010-12-31T23:59:59Z]',
#'      '-article_type:correction'),
#'    limit=50)
#' dois <- dois$data$id[!grepl("annotation", dois$data$id)]
#' alm <- alm_ids(doi=dois, total_details=TRUE)
#' alm <- ldply(alm$data, data.frame)
#' plot_density(alm)
#' plot_density(alm, color="#DCA121")
#' plot_density(alm, title="Scopus citations from 2010")
#' plot_density(alm, title="Scopus citations from 2010", description="Probablity of
#'    X number of citations for a paper")
#' plot_density(alm, description="Probablity of X number of citations for a paper")
#' plot_density(alm, source="crossref_total")
#' plot_density(alm, source="twitter_total")
#' plot_density(alm, source="counter_total")
#' plot_density(alm, source=c("counter_total","facebook_likes"))
#' plot_density(alm, source=c("counter_total","facebook_likes"))
#' plot_density(alm, source=c("counter_total","facebook_likes", "twitter_total"))
#' plot_density(alm, source=c("counter_total","crossref_total","twitter_total"),
#'    color=c("#DBAC6A", "#E09B33", "#A06D34"))
#' plot_density(alm, source=c("counter_total","crossref_total",
#'    "twitter_total"))
#' plot_density(alm, source=c("counter_total","crossref_total"),
#'    title="Counter and Crossref")
#' plot_density(alm, source=c("counter_total","crossref_total",
#'    "twitter_total"), color=c("#83DFB4","#EFA5A5","#CFD470"))
#' }


plot_density <- function(input, source="scopus_total", color="#1447f2", title = "",
                         description = "", plot_type="density") {

  plos_color <- "#1447f2"
  input$date_modified <- as.Date(input$date_modified)
  plot_type <- match.arg(plot_type, choices=c("histogram","density"))
  yvarname <- switch(plot_type,
                      histogram = "No. of articles",
                      density = "Probability")
  plot_type <- switch(plot_type,
                      histogram = "geom_histogram",
                      density = "geom_density")

  # Calculate quantiles and min/max x-axis limits
  minmax <- c(0, round_any(max(input[,source]), 10, f=ceiling))

  # Determine margin at top of plot
  if(nchar(title)==0 && nchar(description)!=0 | nchar(description)==0 && nchar(title)!=0){
    plot_margin <- 2
  } else if(nchar(description)!=0 && nchar(title)!=0){
    plot_margin <- 4
  } else if(nchar(description)==0 && nchar(title)==0){
    plot_margin <- 1
  }

  theme_density <- function(){
    list(theme_bw(base_size=16),
         theme(panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),
               panel.border=element_blank(),
               axis.line = element_line(color = 'black'),
               plot.title = element_text(colour=plos_color, size=24, hjust=0),
               plot.margin = unit(c(plot_margin,1,1,1), "cm")),
         labs(x="", y=""))
  }
  # Plot
  if(length(source)==1){
    plos_xlab <- alm_capwords(gsub("_"," ",source))
    p <- ggplot(input, aes_string(x=source)) +
      theme_density() +
      eval(parse(text=plot_type))(fill=color, colour=color) +
      scale_x_continuous(limits=minmax)
    grid.newpage()
    print(p, vp = viewport(width = 1, height = 1))
    grid.text(yvarname, x=0.03, y=0.5, rot=90, gp=gpar(cex=1, col=plos_color))
    grid.text(plos_xlab, x=0.5, y=0.04, gp=gpar(cex=1, col=plos_color))
    if(nchar(title)==0 & nchar(description)!=0){
      grid.text(paste(strwrap(description,width=90), collapse="\n"), x=0.1, y=0.95, just="left", gp=gpar(cex=1))
    } else if(nchar(description)==0 & nchar(title)!=0){
      grid.text(title, x=0.1, y=0.95, just="left", gp=gpar(cex=2, col=plos_color))
    } else if(nchar(description)!=0 & nchar(title)!=0){
      grid.text(title, x=0.1, y=0.95, just="left", gp=gpar(cex=2, col=plos_color))
      grid.text(paste(strwrap(description,width=90), collapse="\n"), x=0.1, y=0.87, just="left", gp=gpar(cex=1))
    }
  } else
  {
    if(length(color)==1){
      colors <- rep(color,length(source))
    } else
    {
      colors <- color
    }
    df <- input[,source]
    dfm <- melt(df)
    p <- ggplot(dfm, aes(x=value, fill=variable, colour=variable)) +
      theme_density() +
      eval(parse(text=plot_type))() +
      theme(legend.position="none") +
      facet_wrap(~variable, scales="free") +
      scale_colour_manual(values = colors) +
      scale_fill_manual(values = colors)
    grid.newpage()
    print(p, vp = viewport(width = 1, height = 1))
    grid.text(yvarname, x=0.03, y=0.5, rot=90, gp=gpar(cex=1))
#     grid.text("Count", x=0.5, y=0.04, gp=gpar(cex=1))
    if(nchar(title)==0 & nchar(description)!=0){
      grid.text(paste(strwrap(description,width=90), collapse="\n"), x=0.1, y=0.95, just="left", gp=gpar(cex=1))
    } else if(nchar(description)==0 & nchar(title)!=0){
      grid.text(title, x=0.1, y=0.95, just="left", gp=gpar(cex=2))
    } else if(nchar(description)!=0 & nchar(title)!=0){
      grid.text(title, x=0.1, y=0.95, just="left", gp=gpar(cex=2))
      grid.text(paste(strwrap(description,width=90), collapse="\n"), x=0.1, y=0.87, just="left", gp=gpar(cex=1))
    }
  }
}
