#' The GMPE of BCHydro 2018 for subduction earthquakes (incoporated regional updates for Japan)
#'
#' This function calculates the ground motion median values and standard deviations
#' @param M Moment magnitude, a numeric value
#' @param T Period (sec) (range should be between 0 s and 10 s). PSA at 0.01 s will be return for
#' any T <= 0.01 s.
#' Use 1000 (by default) for output the array of Sa with original BCHyrdo per-defined periods
#' @param Rrup Closest distance (km) to the ruptured plane. Use for Interface.
#' @param Rhypo Hypocentral distance(km). Use for Intraslab.
#' @param F_faba Flag for BackArc sites: 0 for Forearc or unknown sites; 1 for Backarc sites
#' @param Vs30 The time averaged shear wave velocity on the top 30 m.
#' @param F_event Flag for interface event or intraslab events:
#' 0 for Interface - use rupture distance; 1 for Intraslab - use hypocentral distance.
#' @param Zh The depth of hypocenter (km)
#' @return A list of five elements is returned: Sa - median spectral acceleration prediction (in g);
#' sigma - totla logarithmic standard deviation (log);
#' phi - within event residuals logarithmic standard deviation (log);
#' tau - between event residuals logarithmic standard deviation (log);
#' period - the corresponding oscillator periods
#' @examples bchydro_2018(M = 7.8, T = 1000, Rrup = 85, Rhypo = 100, F_faba = 0, Vs30 = 300, F_event = 0,
#' Zh = 25)
#'
#' bchydro_2018(M = 7.8, T = c(0.010, 0.020, 0.030, 0.050, 0.075, 0.100, 0.150, 0.200, 0.250, 0.300,
#' 0.400, 0.500, 0.600, 0.750, 1.000, 1.500, 2.000, 2.500, 3.000, 4.000, 5.000), Rrup = 85,
#' Rhypo = 100, F_faba = 1, Vs30 = 300, F_event = 1, Zh = 65)
#' @references Abrahamson, N., Gregor, N., and Addo, K. (2016). BC Hydro Ground Motion Prediction Equations
# for Subduction Earthquakes. Earthquake Spectra. 32(1), 23-44.
#' @export
#' @importFrom stats approx
bchydro_2018 <- function(M, T = 1000, Rrup, Rhypo, F_faba, Vs30, F_event, Zh){

  if(F_event == 0){
    R = Rrup  # interface
  }else{
    R = Rhypo  # intraslab
  }

  period <- c(0.010, 0.020, 0.030, 0.050, 0.075, 0.100, 0.150, 0.200, 0.250, 0.300,
              0.400, 0.500, 0.600, 0.750, 1.000, 1.500, 2.000, 2.500, 3.000, 4.000,
              5.000, 6.000, 7.500, 10.000)
  Vlin <- c(865.1, 865.1, 907.8, 1053.5, 1085.7, 1032.5, 877.6, 748.2, 654.3, 587.1,
            503.0, 456.6, 430.3, 410.5, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0,
            400.0, 400.0, 400.0)
  b <- c(-1.186, -1.219, -1.273, -1.346, -1.471, -1.624, -1.931, -2.188, -2.381, -2.518,
         -2.657, -2.669, -2.599, -2.401, -1.955, -1.025, -0.299, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
  n <- 1.18
  c <- 1.88
  C4 = 10
  CAS_inter = 1.056
  CAS_slab = 0.829
  theta1 = c(3.78779, 3.89949, 4.19326, 4.60231, 5.09669, 5.17563, 4.90983, 4.5761, 4.32152, 4.11569,
             3.71588, 3.33302, 3.02101, 2.54275, 1.9678, 1.09764, 0.43351, -0.05315, -0.45424, -0.91168,
             -1.34684, -1.70854, -2.09138, -2.43297)
  theta2 = c(-1.017, -1.044, -1.112, -1.11, -1.11, -1.11, -1.084, -1.027, -0.983, -0.947, -0.89, -0.845,
             -0.809, -0.76, -0.698, -0.612, -0.55, -0.501, -0.46, -0.455, -0.45, -0.45, -0.45, -0.45)
  theta3 = rep(0.1, (length(period)))
  theta4 = c(0.59, 0.59, 0.59, 0.59, 0.59, 0.59, 0.59, 0.62, 0.64, 0.66, 0.68, 0.68, 0.68, 0.68, 0.68,
             0.68, 0.68, 0.68, 0.68, 0.68, 0.73, 0.78, 0.84, 0.93)
  theta5 = rep(0, length(period))
  theta6 = c(-0.0082, -0.00816, -0.00819, -0.00902, -0.00938, -0.0092, -0.00865, -0.00838, -0.0081,
             -0.00791, -0.00721, -0.00663, -0.00631, -0.0059, -0.00562, -0.00545, -0.00501, -0.00512,
             -0.00528, -0.00542, -0.00532, -0.00483, -0.00468, -0.0052)
  theta7 = rep(0, length(period))
  theta8 = rep(0, length(period))
  theta9 = rep(0.4, length(period))
  theta10 = rep(1.73, length(period))
  theta11 = c(0.017, 0.017, 0.017, 0.018, 0.018, 0.018, 0.0175, 0.017, 0.016, 0.0152, 0.014, 0.013,
              0.0122, 0.0113, 0.01, 0.0082, 0.007, 0.006, 0.0052, 0.004, 0.003, 0.0022, 0.0013, 0)
  theta12 = c(0.95809, 1.00907, 1.11174, 1.33512, 1.60945, 1.69169, 1.76427, 1.88192, 2.01389, 2.09567,
              2.21036, 2.19004, 2.08569, 1.79636, 1.28803, 0.31333, -0.45672, -0.76367, -0.74076, -0.67785,
              -0.63657, -0.61146, -0.57038, -0.49203)
  theta13 = c(-0.0135, -0.0135, -0.0135, -0.0138, -0.0142, -0.0145, -0.0153, -0.0162, -0.0172, -0.0183,
              -0.0206, -0.0231, -0.0256, -0.0296, -0.0363, -0.0493, -0.061, -0.0711, -0.0798, -0.0935,
              -0.098, -0.098, -0.098, -0.098)
  theta14 = c(-0.223, -0.196, -0.128, -0.13, -0.13, -0.13, -0.156, -0.172, -0.184, -0.194, -0.21, -0.223,
              -0.233, -0.245, -0.261, -0.285, -0.301, -0.313, -0.323, -0.282, -0.25, -0.25, -0.25, -0.25)
  theta15 = rep(0, length(period))
  theta16 = rep(0, length(period))
  C1_inter = c(8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.15, 8.1, 8.05, 8, 7.95,
               7.9, 7.85, 7.8, 7.8, 7.8, 7.8)
  C1_slab = rep(7.2, length(period))
  C1 = F_event*C1_slab+(1-F_event)*C1_inter
  phi = 0.6
  tau = 0.43
  sigma_ss = 0.6
  base = rep(0, length(period))
  f_mag = rep(0, length(period))
  f_dep = rep(0, length(period))
  f_faba = rep(0, length(period))
  f_site = rep(0, length(period))
  if(Vs30 > 1000){
    Vs_star = 1000
  }else{
    Vs_star = Vs30
  }
  for(i in seq(1, length(period))){
    base[i] = theta1[i] + theta4[i]*(C1_slab[i] - C1_inter[i])*F_event +
      (theta2[i] + theta14[i]*F_event + theta3[i]*(M - 7.8))*log(R + C4*exp((M - 6)*theta9[i])) +
      theta6[i]*R + theta10[i]*F_event
    if(M <= C1[i]){
      f_mag[i] = theta4[i]*(M - (C1[i])) + theta13[i]*(10 - M)^2
    }else{
      f_mag[i] = theta5[i]*(M - (C1[i])) + theta13[i]*(10 - M)^2
    }
    if(Zh <= 100){
      f_dep[i] = theta11[i]*(Zh - 60)*F_event
    }else{
      f_dep[i] = theta11[i]*(100-60)*F_event
    }
    if(F_event){
      f_faba[i] = (theta7[i] + theta8[i]*log(max(R, 85)/40))*F_faba
    }else{
      f_faba[i] = (theta15[i] + theta16[i]*log(max(R, 100)/40))*F_faba
    }
    if(i == 1){
      if(F_event){
        PGA1000 = exp(base[1] + f_mag[1] + f_dep[1] + f_faba[1] +
                        theta12[1]*log(1000/Vlin[1]) + b[1]*n*log(1000/Vlin[1])+ CAS_slab)
      }else{
        PGA1000 = exp(base[1] + f_mag[1] + f_dep[1] + f_faba[1] + theta12[1]*log(1000/Vlin[1]) +
                        b[1]*n*log(1000/Vlin[1])+ CAS_inter)
      }
    }
    if(Vs30 < Vlin[i]){
      f_site[i] = theta12[i]*log(Vs_star/Vlin[i]) - b[i]*log(PGA1000 + c) +
        b[i]*log(PGA1000 + c*(Vs_star/Vlin[i])^n)
    }else{
      f_site[i] = theta12[i]*log(Vs_star/Vlin[i]) + b[i]*n*log(Vs_star/Vlin[i])
    }
  }
  sa = exp(base + f_mag + f_dep + f_faba + f_site)

  # output results at pre-defined periods
  if (length(T) == 1 & T[1] == 1000) {

    res <- list()

    res$Sa <- sa

    res$sigma <- rep(sqrt(phi^2 + tau^2), length(period))

    res$phi <- rep(phi, length(period))

    res$tau <- rep(tau, length(period))

    res$period <- period

    return(res)

  } else {

    # handle with periods other than pre-defined periods
    I = (T <= 0.01) # <= 0.01 s

    sa_int = approx(log(c(1e-10, period[-1])), sa, log(T), rule = 2)$y

    sa_int[I] = sa[1]   # Set small periods to PGA

    res <- list()

    res$Sa <- sa_int

    res$sigma <- rep(sqrt(phi^2 + tau^2), length(T))

    res$phi <- rep(phi, length(T))

    res$tau <- rep(tau, length(T))

    res$period <- T

    return(res)
  }
}


