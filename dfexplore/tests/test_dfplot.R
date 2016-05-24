library(dfexplore)
#test_dfplot
data(example_df)
dfplot(example_df)

# Plot df with matrix
df_with_matrix <- simulate_dataframe(includeMatrix=T)
str(df_with_matrix)
dfplot(df_with_matrix)
