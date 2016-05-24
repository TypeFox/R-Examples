ANOMgen <- function(mu, n=NULL, gm=NULL, lo, up, names=NULL, alternative="two.sided",
                    xlabel="Group", ylabel="Endpoint", printn=T, p=NULL, bg="white"){
  
  if(is.null(gm)){
    gm <- weighted.mean(mu, n)
  }
  
  if(!(is.null(n)) & !(is.null(gm))){
    check <- abs(weighted.mean(mu, n) - gm)
    if(check > 0.01){
      stop("The 'gm' value you inserted is not the grand mean\n
           computed from your 'mu' and 'n' values!")
    }
    }
  
  grp <- 1:length(mu)
  grpf <- as.factor(grp)
  
  if(!(is.null(names))){
    names <- as.factor(names)
  }else{
    names <- grpf
  }
  
  if(is.null(n)){
    printn <- F
  }
  
  if(!(is.null(p))){
    printp <- T
    p <- paste("p=", round(p, 3), sep="")
    p[p=="p=0"] <- "p<0.001"
  }else{
    printp <- F
  }
  
  dir <- match.arg(alternative, choices=c("two.sided", "greater", "less"))
  
  if(dir=="two.sided"){
    if(all(lo==-Inf) | all(lo==-1) | all(lo==0) | all(up==1) | all(up==Inf)){
      warning("Is your alternative really 'two.sided'? Doesn't seem so!")
    }
  }
  
  ldl <- gm - abs((mu - gm) - up)
  udl <- gm + abs((mu - gm) - lo)
  
  if(printp==T){
    if(printn==T){
      set <- data.frame(mu, gm, lo, up, ldl, udl, grp, grpf, names, n, p)
    }else{
      set <- data.frame(mu, gm, lo, up, ldl, udl, grp, grpf, names, p)
    }
  }else{
    if(printn==T){
      set <- data.frame(mu, gm, lo, up, ldl, udl, grp, grpf, names, n)
    }else{
      set <- data.frame(mu, gm, lo, up, ldl, udl, grp, grpf, names)
    }
  }
  
  bg <- match.arg(bg, choices=c("gray", "grey", "white"))
  
  if(bg=="white"){
    back <- theme_bw()
  }else{
    back <- theme_gray()
  }
  
  if(dir=="two.sided"){
    
    basic <- ggplot(set, aes(x=grpf, y=mu)) +
      scale_x_discrete(labels=names) +
      xlab(xlabel) +
      ylab(ylabel) +
      geom_rect(aes(xmin=grp-0.5, xmax=grp+0.5, ymin=ldl, ymax=udl),
                alpha=0.5, fill="darkgray", linetype=0) +
      geom_segment(aes(x=grp-0.5, xend=grp+0.5, y=udl, yend=udl), size=1, linetype=2) +
      geom_segment(aes(x=grp-0.5, xend=grp+0.5, y=ldl, yend=ldl), size=1, linetype=2) +
      geom_segment(aes(x=0.5, xend=max(grp)+0.5, y=gm, yend=gm), size=1.5) +
      geom_segment(aes(x=grp, xend=grp, y=mu, yend=gm), size=1) +
      geom_point(size=5) +
      annotate("text", label="LDL", x=max(grp)+0.4, y=ldl[max(grp)], size=4, vjust=1.5) +
      annotate("text", label="UDL", x=max(grp)+0.4, y=udl[max(grp)], size=4, vjust=-0.75) +
      ylim(min(min(mu), min(ldl))-(gm-min(min(mu), min(ldl)))/5,
           max(max(mu), max(udl))+(max(max(mu), max(udl))-gm)/5) +
      back +
      theme(axis.text.x=element_text(size=18), axis.text.y=element_text(size=18),
            axis.title.x=element_text(size=25), axis.title.y=element_text(size=25))
    
  }
  
  if(dir=="greater"){
    
    basic <- ggplot(set, aes(x=grpf, y=mu)) +
      scale_x_discrete(labels=names) +
      xlab(xlabel) +
      ylab(ylabel) +
      geom_rect(aes(xmin=grp-0.5, xmax=grp+0.5, ymin=ldl, ymax=udl),
                alpha=0.5, fill="darkgray", linetype=0) +
      geom_segment(aes(x=grp-0.5, xend=grp+0.5, y=udl, yend=udl), size=1, linetype=2) +
      geom_segment(aes(x=0.5, xend=max(grp)+0.5, y=gm, yend=gm), size=1.5) +
      geom_segment(aes(x=grp, xend=grp, y=mu, yend=gm), size=1) +
      geom_point(size=5) +
      annotate("text", label="UDL", x=max(grp)+0.4, y=udl[max(grp)], size=4, vjust=-0.75) +
      ylim((min(mu)-(gm-min(mu))/5)[1],
           (max(max(mu), max(udl))+(max(max(mu), max(udl))-gm)/5)[1]) +
      back +
      theme(axis.text.x=element_text(size=18), axis.text.y=element_text(size=18),
            axis.title.x=element_text(size=25), axis.title.y=element_text(size=25))
    
  }
  
  if(dir=="less"){
    
    basic <- ggplot(set, aes(x=grpf, y=mu)) +
      scale_x_discrete(labels=names) +
      xlab(xlabel) +
      ylab(ylabel) +
      geom_rect(aes(xmin=grp-0.5, xmax=grp+0.5, ymin=ldl, ymax=udl),
                alpha=0.5, fill="darkgray", linetype=0) +
      geom_segment(aes(x=grp-0.5, xend=grp+0.5, y=ldl, yend=ldl), size=1, linetype=2) +
      geom_segment(aes(x=0.5, xend=max(grp)+0.5, y=gm, yend=gm), size=1.5) +
      geom_segment(aes(x=grp, xend=grp, y=mu, yend=gm), size=1) +
      geom_point(size=5) +
      annotate("text", label="LDL", x=max(grp)+0.4, y=ldl[max(grp)], size=4, vjust=1.5) +
      ylim((min(min(mu), min(ldl))-(gm-min(min(mu), min(ldl)))/5)[1],
           (max(mu)+(max(mu)-gm)/5)[1]) +
      back +
      theme(axis.text.x=element_text(size=18), axis.text.y=element_text(size=18),
            axis.title.x=element_text(size=25), axis.title.y=element_text(size=25))
    
  }
  
  if(dir=="greater"){
    npart <- geom_text(aes(y=min(mu), label=paste("n=", n, sep="")), vjust=3, size=5)
  }else{
    npart <- geom_text(aes(y=min(min(mu), min(ldl)), label=paste("n=", n, sep="")), vjust=3, size=5)
  }
  
  if(dir=="less"){
    ppart <- geom_text(aes(y=max(mu), label=p), vjust=-2, size=5)
  }else{
    ppart <- geom_text(aes(y=max(max(mu), max(udl)), label=p), vjust=-2, size=5)
  }
  
  if(printn==T){
    if(printp==T){
      chart <- basic + npart + ppart
    }else{
      chart <- basic + npart
    }
  }else{ 
    if(printp==T){
      chart <- basic + ppart
    }else{
      chart <- basic
    }
  }
  
  print(chart)
  
}