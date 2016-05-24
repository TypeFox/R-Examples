require("rtfbs")
exampleArchive <- system.file("extdata", "NRSF.zip", package="rtfbs")
pwmFile <- "pwm.meme"
unzip(exampleArchive, pwmFile)
# Read in Position Weight Matrix (PWM) from MEME file from
#  the examples into a Matrix object
pwm <- read.pwm(pwmFile)
# Print PWM as an R matrix
print(pwm)
