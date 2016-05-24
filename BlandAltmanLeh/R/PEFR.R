#' PEFR Data from Bland JM and Altman DG 1986
#' 
#' Peak expiratory flow data from 17 members of Bland's family, taken with two 
#' different instruments, each twice. This data is for explanatory use only.
#' Columns 1 and 2 were measured with the "Wright" peak flow meter, columns 3
#' and 4 with the "Mini Wright" peak flow meter. 
#' These are the data behind fig. 1, fig. 2 and fig. 6 of the original paper and
#' these can be easily reconstructed
#' @examples
#' # this is what fig. 1. would have looked like in R:
#' x <- bland.altman.PEFR[["bigger.first"]]
#' y <- bland.altman.PEFR[["smaller.first"]]
#' plot(x,y, xlab="PEFR by large meter",ylab="PEFR by mini meter", 
#'      xlim=c(0,800), ylim=c(0,800))
#' abline(0,1)
#' @export
#' 
bland.altman.PEFR <- data.frame(
    bigger.first = c(494,395,516,434,476, 557,413,442,650, 433,417,656,267, 
                     478,178,423,427),
    bigger.second= c(490,397,512,401,470, 611,415,431,638, 429,420,633,275,
                     492,165,372,421),
    smaller.first= c(512,430,520,428,500, 600,364,380,658, 445,432,626,260,
                     477,259,350,451),
    smaller.second=c(525,415,508,444,500, 625,460,390,642, 432,420,605,227,
                     467,268,370,443))