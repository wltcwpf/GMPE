test_that("test bssa_2014_nga vs NGA-West2 spreadsheet", {
  expect_equal(bssa_2014_nga(M = 5, T = 1000, Rjb = 85, Fault_Type = 3, region = 1,
                             z1 = 999, Vs30 = 350)$med[c(1, 2, 3, 4, 8, 18, 24, 29, 35, 40, 43, 47,
                                                         53, 59, 65, 70, 75, 80, 86, 92, 97, 102,
                                                         107)],
               c(0.269344818, 0.006198101, 0.006281779, 0.006117412, 0.006158664, 0.006901768,
                 0.009132362, 0.011669154, 0.013643745, 0.013382022, 0.012548628, 0.011509203,
                 0.009160364, 0.007262497, 0.003952132, 0.002409222, 0.001069659, 0.000579513,
                 0.000248159, 0.000143274, 9.44772E-05, 4.33956E-05, 2.25331E-05))

  expect_equal(round(bssa_2014_nga(M = 5, T = 1000, Rjb = 85, Fault_Type = 3, region = 1,
                             z1 = 999, Vs30 = 350)$sigm[c(1, 2, 3, 4, 8, 18, 24, 29, 35, 40, 43, 47,
                                                         53, 59, 65, 70, 75, 80, 86, 92, 97, 102,
                                                         107)], digits = 4),
               c(0.7051, 0.7022, 0.7055, 0.7106, 0.7382, 0.7928, 0.7969, 0.7740, 0.7336, 0.7051,
                 0.6945, 0.6852, 0.6792, 0.6848, 0.7018, 0.7109, 0.7125, 0.7159, 0.7255, 0.7285,
                 0.7201, 0.6934, 0.6648))

})
