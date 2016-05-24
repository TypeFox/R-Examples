water <- rbind(c(14, 15, 8, 4), c(11, 26, 16, 6)) # to combine Tap and Filtered
water
colnames(water) <- c('Aquafina', 'Fiji', 'SamsChoice', 'Tap') # add column names
rownames(water) <- c('Bottled', 'Tap/Filtered') # add row names
water

