results <- cumddply(rt.data,
                    c('Subject', 'Block'),
                    c('Trial'),
                    function (df) {with(df, data.frame(RT = RT, MeanRT = mean(RT)))})
# This doesn't work, because RT has many values in each df, whereas MeanRT is a unitary value.
