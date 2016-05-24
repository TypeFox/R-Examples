sample_seed <- function(n = 1, for_ms = FALSE) {
  max_value <- ifelse(for_ms, 65536, 4294967296) # 16 / 32 bit
  sample.int(max_value, n)
}
