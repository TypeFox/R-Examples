library('cumplyr')

x <- 1:3
y <- 2:4
z <- 3:5

cartesian_product(c('x', 'y', 'z'))
cartesian_product(c('x', 'y'))
cartesian_product(c('x', 'z'))
cartesian_product(c('y', 'z'))
cartesian_product(c('x'))

data(rt.data)

results <- cumddply(rt.data,
                    c('Subject', 'Block'),
                    c('Trial'),
                    function (df) {with(df, mean(RT))})

#all(abs(read.csv('output_data.csv') - results) < 10e-6)
