test_that("test zea_2016_gmpe", {

  zea_2017_SC_UM(mag = 7.5, T = 1000, fDepth = 0.0, dist = 40.0, DistV = 0.0,
                 eqkType = 2, FaultMech = 2, SiteClass = 4)

  zea_2017_SC_UM(mag = 7.5, T = c(0.025, 0.5), fDepth = 0.0, dist = 40.0, DistV = 0.0,
                 eqkType = 2, FaultMech = 2, SiteClass = 4)

  zea_2016_Sub_Interf(mag = 7.5, T = 1000, fDepth = 0.0, dist = 40.0, DistV = 0.0, SiteClass = 4)

  zea_2016_Sub_Interf(mag = 7.5, T = c(0.025, 0.5), fDepth = 0.0, dist = 40.0, DistV = 0.0, SiteClass = 4)

  zea_2016_Sub_Slab(mag = 7.5, T = 1000, fDepth = 0.0, dist = 40.0, DistV = 0.0, SiteClass = 4)

  zea_2016_Sub_Slab(mag = 7.5, T = c(0.025, 0.5), fDepth = 0.0, dist = 40.0, DistV = 0.0, SiteClass = 4)

})
