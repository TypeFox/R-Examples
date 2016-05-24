print.ICA <-
function(x, ...) {
  with(x, {
    cat(" Call:\n")
    print(call)
    cat("\n Best solution: ", position, "\n",
        "Best value: ", value, "\n",
        "No. of Imperialists: ", nimp, "\n",
        "Timings:\n")
    print(time)
  })
  invisible(x)
}
