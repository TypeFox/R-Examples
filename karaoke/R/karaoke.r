karaoke<-function(infile=NULL, outfile=NULL, sampf=NULL) {
wobj<-readWave(infile)
wl<-mono(wobj, "left")
wr<-mono(wobj, "right")
wobj<-wl-wr
wobj<-stereo(wobj,wobj)
savewav(wobj,f=sampf,filename=outfile)
}
