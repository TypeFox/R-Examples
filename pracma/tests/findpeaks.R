##
##  f i n d p e a k s . R  Test suite
##


findpeaks <- pracma::findpeaks

x <- c(2, 12, 4, 6, 9, 4, 3, 1, 19, 7)
identical(findpeaks(x),
          matrix(c(12,9,19, 2,5,9, 1,3,8, 3,8,10), nrow=3, ncol=4))
identical(findpeaks(x, npeaks = 1, sortstr = TRUE),
          matrix(c(19, 9, 8, 10), nrow = 1))
identical(findpeaks(x, minpeakheight = 15),
          matrix(c(19, 9, 8, 10), nrow = 1))
identical(findpeaks(x, threshold = 10),
          matrix(c(19, 9, 8, 10), nrow = 1))
# Not yet implemented
# identical(findpeaks(x, threshold = 10),
#           c(19, 9, 8, 9))

