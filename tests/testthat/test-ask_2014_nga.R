test_that("ask_2014_nga", {
  ask_2014_nga(M = 5, T = 1000, Rrup = 90, Rjb = 85, Rx = 85, Ry0 = 80, dip = 80,
               lambda = 75, after_shock = 0, HW = 1, W = 10, region = 1, Vs30 = 350, Vs30_code = 0)
})
