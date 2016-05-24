require(ggplot2, quietly=TRUE)
require(scales)
require(bdscale)

nyse <- yahoo()

set.seed(12345)
df <- data.frame(date=nyse, price=cumsum(rnorm(length(nyse))) + 100)
df <- subset(df, as.Date('2014-08-01') <= date & date <= as.Date('2014-10-08')) 

plot <- ggplot(df, aes(x=date, y=price)) + geom_step() + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

plot + ggtitle('calendar dates')

# ggsave(file='man/figures/calendar.PNG', width=6, height=2)

plot + scale_x_bd(business.dates=nyse, labels=date_format("%b '%y")) + 
  ggtitle('business dates, month breaks')

# ggsave(file='man/figures/business.month.PNG', width=6, height=2)

plot + scale_x_bd(business.dates=nyse, max.major.breaks=10, labels=date_format('%b %d')) + 
  ggtitle('business dates, week breaks')

# ggsave(file='man/figures/business.week.PNG', width=6, height=2)

options <- as.Date(c('2014-08-15', '2014-09-19'))

plot + 
  geom_vline(xintercept=as.numeric(options), size=2, alpha=0.25) + 
  ggtitle('calendar dates, option expiry')

# ggsave(file='man/figures/calendar.options.PNG', width=6, height=2)

plot + 
  geom_vline(xintercept=bd2t(options, business.dates=nyse), size=2, alpha=0.25) + 
  scale_x_bd(business.dates=nyse) +
  ggtitle('business dates, option expiry')

# ggsave(file='man/figures/business.options.PNG', width=6, height=2)
