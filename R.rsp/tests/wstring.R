library("R.rsp")

x <- sqrt(2)
strings <- c(character(0L), "", NA_character_, "Hello world!", "x={{x}}", "x^2={{x^2}}")

for (s in strings) {
  print(s)
  t <- wstring(s)
  print(t)
}

t <- wstring(strings)
print(t)
