plot.lba <- function(x, 
                     budget.prop = TRUE, 
                     col.points  = NULL,
                     col.lines   = NULL,# only to K = 2
                     col.budget  = NULL,
                     pch.points  = NULL,
                     pch.budget  = NULL,
                     lty.lines   = NULL,# only to K = 2
                     lty.budget  = NULL,
                     lwd.lines   = NULL,# only to K = 2
                     lwd.budget  = NULL,
                     legend      = TRUE, 
                     with.ml     = c('mix','lat'),
                     type        =c("lba","corr"),  
                     ...)
{

 switch(match.arg(with.ml),
        mix = {
         if(is.integer(grep('fe',class(x))) | is.integer(grep('logit',class(x)))) {
          nrowss <- dim(x[[6]])[1]
          alfas <- x[[6]]
          alfas <- alfas/rowSums(alfas)
          pk <- x[[9]]
         } else {
          nrowss <- dim(x[[4]])[1]
          alfas <- x[[4]]
          alfas <- alfas/rowSums(alfas)
          pk <- x[[6]] 
         }
        },
        lat = {
         if(is.integer(grep('fe',class(x))) | is.integer(grep('logit',class(x)))) { 
          nrowss <- dim(x[[8]])[1]
          alfas <- x[[8]]
          alfas <- alfas/rowSums(alfas)
          pk <- x[[9]]
         } else {
          nrowss <- dim(x[[6]])[1]
          alfas <- x[[6]]
          alfas <- alfas/rowSums(alfas)
          pk <- x[[6]] 
         }
        }
        )

 if(ncol(alfas)==1) stop('Only K between 2 and 3 budgets') 
 switch(match.arg(type),
        #latent budget analysis#
        lba = {

         rlabels <- rownames(alfas) 
         K <- ncol(alfas) 

         if(is.null(col.points)){
          col.points    <- rep(1,nrowss)
         } else { col.points <- col.points } 

         if(is.null(col.lines)){
          col.lines    <- rep(1,nrowss)
         } else { col.lines <- col.lines } 

         if(is.null(col.budget)){
          col.budget    <- 1
         } else { col.budget <- col.budget } 

         if(is.null(pch.points)){
          pch.points    <- as.character(1:nrowss)
         } else { pch.points <- pch.points }

         if(is.null(pch.budget)){
          pch.budget <- 1
         } else {pch.budget <- pch.budget }

         if(is.null(lty.lines)){
          lty.lines <- rep(1,nrowss)
         } else {lty.lines <- lty.lines }

         if(is.null(lty.budget)){
          lty.budget <- 1
         } else {lty.budget <- lty.budget}

         if(is.null(lwd.lines)){
          lwd.lines <- rep(1,nrowss)
         } else {lwd.lines <- lwd.lines }

         if(is.null(lwd.budget)){
          lwd.budget <- 1
         } else {lwd.budget <- lwd.budget }

         if (K == 2){ 
          plot.new()
          plot.window(seq(0, 1), 
                      seq(0, 1))
          axis(1, 
               at = seq(0, 
                        1, 
                        0.1))
          y <- seq(0, 
                   0.8, 
                   length = nrowss)
          for (i in 1:nrowss) {
           segments(alfas[i, 1], 
                    0, 
                    alfas[i, 1], 
                    y[i],
                    col = col.lines[i],
                    lty = lty.lines[i],
                    lwd = lwd.lines[i],
                    ...)
           points(alfas[i, 1], 
                  y[i]+0.02, 
                  pch = pch.points[i],
                  col = col.points[i],
                  ...)
          }
          if(budget.prop){
           segments(pk[1,1], 
                    0, 
                    pk[1,1], 
                    1, 
                    lty = lty.budget, 
                    lwd = lwd.budget, 
                    ...)
           points(pk[1, 1],
                  1,
                  col = col.budget,
                  pch = pch.budget, 
                  ...)
          }
          if(legend){

           legend('topleft',
                  rlabels,
                  pch = pch.points)
          }

         } else {

          if(ncol(alfas)>3)stop('Use corr type to options graphics. See documentation!')

          triax.plot(alfas[,rev(order(colnames(alfas)))], 
                     #show.grid = TRUE,
                     col.symbols = col.points,
                     pch = pch.points,
                     #show.legend = FALSE,
                     ...)
          if(budget.prop){
           triax.points(matrix(pk[,rev(order(colnames(pk)))],
                               ncol = 3),
                        show.legend = FALSE,
                        col.symbols = col.budget,
                        pch = pch.budget)
           segments(x0 = 0, 
                    y0 = 0, 
                    x1 = 1-(pk[,3] + pk[,1] * 0.5), 
                    y1 = pk[,1] * sin(pi/3),
                    lty = lty.budget,
                    lwd = lwd.budget)
           segments(x0 = 1,
                    y0 = 0, 
                    x1 = 1-(pk[,3] + pk[,1] * 0.5), 
                    y1 = pk[,1] * sin(pi/3),
                    lty = lty.budget,
                    lwd = lwd.budget)
           segments(x0 = 1-(pk[,3] + pk[,1] * 0.5), 
                    y0 = pk[,1] * sin(pi/3), 
                    x1 = 0.5, 
                    y1 = 0.865,
                    lty = lty.budget,
                    lwd = lwd.budget)
          }  
         }
        },
        # correspondence analysis #
        corr = { 
         if(ncol(alfas)==2)stop('Use lba type to options graphics. See documentation!')
         plot(ca(alfas),
              ...)
        }
        )   
}
