# Run accelerometry code from paper in The R Journal

# Load four-column matrix with counts and steps over 7 days
data(tridata)

readline("Press <Enter> to continue")

# Generate basic PA variables using default settings
dailyPA1 <- accel.process.tri(counts.tri = tridata[, 1:3], steps = tridata[, 4])

readline("Press <Enter> to continue")

# Request full set of PA variables, and use triaxial vector magnitude for non-wear
# detection rather than vertical axis, with 90-minute rather than 60-minute window
dailyPA2 <- accel.process.tri(counts.tri = tridata[, 1:3], steps = tridata[, 4],
                              brevity = 3, nonwear.axis = "mag", nonwear.window = 90)

readline("Press <Enter> to continue")

# Check variable names for dailyPA1 and dailyPA2
colnames(dailyPA1)
colnames(dailyPA2)

readline("Press <Enter> to continue")

# Print contents of dailyPA1 and first 15 variables in dailyPA2
dailyPA1
dailyPA2[, 1:15]

readline("Press <Enter> to continue")

# Calculate average for cpm_vert from dailyPA1 and dailyPA2
mean(dailyPA1[, "cpm_vert"])
mean(dailyPA2[, "cpm_vert"])