textPrompt <- R.cache:::.textPrompt

ans <- textPrompt("Do you have a minute?")
print(ans)

ans <- textPrompt("Do you have a minute?", options=c("yes", "no"))
print(ans)

## Output to standard error
ans <- textPrompt("Do you have a minute?", type="message")
print(ans)

## Output to standard output
ans <- textPrompt("Do you have a minute?", type="output")
print(ans)
