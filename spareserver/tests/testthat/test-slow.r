
context("Slow or unreliable network")

OS <- Sys.info()["sysname"]

tc_init <- function() {
  system("sudo modprobe ifb")
  system("sudo ip link set dev ifb0 up")
  system("sudo tc qdisc add dev eth0 ingress")
  system("sudo tc filter add dev eth0 parent ffff: protocol ip u32 match u32 0 0 flowid 1:1 action mirred egress redirect dev ifb0")
  system("sudo ifconfig ifb0 192.168.1.23 netmask 255.255.255.0 up")
}

tc_reset <- function() {
  system("sudo tc qdisc del dev ifb0 root")
  system("sudo tc qdisc del dev eth0 handle ffff: ingress")
  system("sudo ifconfig ifb0 down")
}

tc_slow <- function(ms_delay = 100) {
  cmd <- sprintf("sudo tc qdisc add dev ifb0 root netem delay %ims", ms_delay)
  system(cmd)
}

tc_unreliable <- function(drop_percent = 25) {
  cmd <- sprintf("sudo tc qdisc add dev ifb0 root netem loss %g%%", drop_percent)
  system(cmd)
}

test_that("slow networks are fine", {

  if (OS != "Linux") skip("Slow network test only works on Linux")

  ## Does not work on Travis :(
  if (identical(Sys.getenv("TRAVIS"), "true")) skip("On Travis")

  assign("spare_services", list(), envir = asNamespace("spareserver"))
  
  tc_init()
  on.exit(tc_reset(), add = TRUE)
  tc_slow(500)

  add_service(
    "test", 
    server("http://google.com", priority = 10, timeout = 5),
    server("http://192.0.2.1/foobar", priority = 5, timeout = 5)
  )

  q <- spare_q("test", "", httr::GET)
  expect_equal(httr::status_code(q), 200)
  
})

test_that("unreliable networks are fine, too", {

  if (OS != "Linux") skip("Unreliable network test only works on Linux")

  ## Does not work on Travis :(
  if (identical(Sys.getenv("TRAVIS"), "true")) skip("On Travis")

  assign("spare_services", list(), envir = asNamespace("spareserver"))
  
  tc_init()
  on.exit(tc_reset(), add = TRUE)
  tc_unreliable(10)

  add_service(
    "test", 
    server("http://google.com", priority = 10),
    server("http://192.0.2.1/foobar", priority = 5)
  )

  q <- spare_q("test", "", httr::GET)
  expect_equal(httr::status_code(q), 200)

})
