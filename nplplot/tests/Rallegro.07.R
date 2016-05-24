library(nplplot)
if(nplplot.old(c("RALLEGRO.07"), col=1, row=1, mode="p", output="Allegro.07.pdf", yline=2.00, ymin=0.00, ymax=3.00, batch=TRUE, yfix=FALSE, titles="Chromosome 7", headerfiles="RALLEGRO.07.hdr"
, lgnd=TRUE)==F) {
   system("/bin/rm Allegro.07.pdf")
}
