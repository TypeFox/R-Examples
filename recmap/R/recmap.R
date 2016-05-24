
.draw_recmap_us_state_ev <- function(plot=TRUE){
  
  cm <- c("#FF0000", "#FF0505", "#FF0A0A", "#FF1010", "#FF1515", "#FF1A1A", "#FF1F1F",
          "#FF2424", "#FF2A2A", "#FF2F2F", "#FF3434", "#FF3939", "#FF3E3E", "#FF4444",
          "#FF4949", "#FF4E4E", "#FF5353", "#FF5858", "#FF5E5E", "#FF6363", "#FF6868",
          "#FF6D6D", "#FF7272", "#FF7878", "#FF7D7D", "#FF8282", "#FF8787", "#FF8D8D",
          "#FF9292", "#FF9797", "#FF9C9C", "#FFA1A1", "#FFA7A7", "#FFACAC", "#FFB1B1",
          "#FFB6B6", "#FFBBBB", "#FFC1C1", "#FFC6C6", "#FFCBCB", "#FFD0D0", "#FFD5D5",
          "#FFDBDB", "#FFE0E0", "#FFE5E5", "#FFEAEA", "#FFEFEF", "#FFF5F5", "#FFFAFA",
          "#FFFFFF", "#FFFFFF", "#FAFAFF", "#F5F5FF", "#EFEFFF", "#EAEAFF", "#E5E5FF",
          "#E0E0FF", "#DBDBFF", "#D5D5FF", "#D0D0FF", "#CBCBFF", "#C6C6FF", "#C1C1FF",
          "#BBBBFF", "#B6B6FF", "#B1B1FF", "#ACACFF", "#A7A7FF", "#A1A1FF", "#9C9CFF",
          "#9797FF", "#9292FF", "#8D8DFF", "#8787FF", "#8282FF", "#7D7DFF", "#7878FF",
          "#7272FF", "#6D6DFF", "#6868FF", "#6363FF", "#5E5EFF", "#5858FF", "#5353FF",
          "#4E4EFF", "#4949FF", "#4444FF", "#3E3EFF", "#3939FF", "#3434FF", "#2F2FFF",
          "#2A2AFF", "#2424FF", "#1F1FFF", "#1A1AFF", "#1515FF", "#1010FF", "#0A0AFF",
          "#0505FF", "#0000FF")
  
  recmap_us_state_ev.file <- system.file("extdata", 
                                         "recmap_us_state_ev.polygon", 
                                         package = "recmap")
  
  recmap_us_state_ev <- read.table(recmap_us_state_ev.file, sep = '|', 
                                   col.names=c('x', 'y'))
  
  us_state_election_2004.file <- system.file("extdata", 
                                             "us_state_election_2004.csv", 
                                             package = "recmap")
  
  us_state_election_2004 <- read.table(us_state_election_2004.file, 
                                       sep = ',')
  if(plot){
    plot(recmap_us_state_ev$x, recmap_us_state_ev$y, type='n', asp = 1, xlab='', ylab='', axes=FALSE)
    
    polygon(recmap_us_state_ev$x, recmap_us_state_ev$y,
            col=cm[round(length(cm)*(us_state_election_2004$V8/(us_state_election_2004$V8 + us_state_election_2004$V9)))+1])
    
    text(us_state_election_2004$V1,us_state_election_2004$V2-7,
         as.character(us_state_election_2004$V3),cex=round(us_state_election_2004$V5*20)/100,
         lwd=2.5,
         pos=3,
         col="black");
  }
  
  # res <- data.frame()
}


.checker_board <- function(n = 8, ratio = 4){
  xy <- (t(combn(1:n, 2)))
  xy <- rbind(cbind(xy[,1], xy[,2]), cbind(xy[,2], xy[,1]), cbind(1:n, 1:n))
  
  
  z.bool <- (xor(xy[,1] %% 2 == 1 , xy[,2] %% 2 == 0))
  z <- rep(1, length(xy[,1]))
  
  z[which(z.bool)] <- z[which(z.bool)] * ratio
  z[which(!z.bool)] <- z[which(!z.bool)] 
  
  res <- data.frame(x = xy[, 1], 
                    y = xy[,2], 
                    dx=0.5, 
                    dy=0.5, 
                    z=z, 
                    name=paste(letters[1:n][xy[,1]], xy[,2], sep=''))
  

  
  res <- res[with(res, order(x, y)), ]
  row.names(res) <- 1:nrow(res); # paste(letters[1:n][xy[,1]], xy[,2], sep='')
  class(res) = c('data.frame', 'recmapFrame')
  res
}


plot_recmap <- function(S, col='#00000011', col.text = 'grey', ...){
  plot(S$x, S$y, 
       xlim = c(min(S$x - S$dx), max(S$x + S$dx)), 
       ylim = c(min(S$y - S$dy), max(S$y + S$dy)), 
       type = 'n', 
       asp=1,
       xlab = '',
       ylab = '',
       axes = FALSE)
  
  # col.idx <- (length(colormap) -1  ) * (S$z - min(S$z) / (max(S$z) - min(S$z))) + 1
  
  rect(xleft = S$x - S$dx, 
       ybottom = S$y - S$dy,  
       xright = S$x + S$dx, 
       ytop = S$y + S$dy, 
       col = col, 
       border = 'darkgreen', ...)
  
  if (sqrt(length(S$x)) < 10){
    text(S$x, S$y, 
         S$name,
         cex=1.5/sqrt(sqrt(length(S$x))),
         col = col.text)
    
    text(S$x, S$y, 
         S$dfs_num,
         col = col.text, 
         pos=1, 
         cex=1/sqrt(sqrt(length(S$x))))
    }
}


