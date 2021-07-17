#' The GMPE of CB 2014
#'
#' This function calculates the ground motion median values and standard deviations
#' @param M Moment magnitude, a numeric value
#' @param T Period (sec); Use Period = -1 for PGV computation.
#' Use 1000 (by default) for output the array of Sa with original NGA West2 periods
#' @param Rrup Closest distance (km) to the ruptured plane
#' @param Rjb Joyner-Boore distance (km); closest distance (km) to surface
#' projection of rupture plane
#' @param Rx Site coordinate (km) measured perpendicular to the fault strike
#' from the fault line with down-dip direction to be positive
#' @param W Down-dip rupture width (km). 999 if unknown
#' @param Ztor Depth(km) to the top of ruptured plane. 999 if unknown
#' @param Zbot Depth(km) to the bottom of the seismogenic crust. It is needed only when W is unknown.
#' @param dip Fault dip angle (in degree)
#' @param lambda Rake angle (in degree)
#' @param Fhw Falg for hanging wall sites. 1 for sites on the hanging wall side of the fault, 0 otherwise.
#' @param Vs30 Shear wave velocity averaged over top 30 m (in m/s)
#' @param Z25 Basin depth (km): depth from the groundsurface to the 2.5km/s shear-wave horizon.
#' 999 if unknown.
#' @param Zhyp Hypocentral depth of the earthquake measured from sea level. 999 in unknown
#' @param region Region indicator: 0 for global (incl. Taiwan); 1 for California; 2 for Japan;
#' 3 for China or Turkey; 4 for Italy
#' @return A list of five elements is returned: med - median spectral acceleration prediction (in g);
#' sigma - logarithmic standard deviation of spectral acceleration prediction; phi - logarithmic
#' standard deviation of within event residuals; tau - logarithmic standard deviation of between
#' event residuals; period - the corresponding oscillator periods
#' @examples cb_2014_nga(M = 5, T = 1000, Rrup = 90, Rjb = 85, Rx = 85, W = 10, Zbot = 15, dip = 80,
#' lambda = 75, Fhw = 1, Vs30 = 350, region = 1, Zhyp = 8, Z25 = 999)
#' @references Campbell, K. W., and Bozorgnia, Y. (2014). NGA-West2 Ground Motion Model
#' for the Average Horizontal Components of PGA, PGV, and 5% Damped Linear
#' Acceleration Response Spectra. Earthquake Spectra, 30(3), 1087-1115.
#' @export
#' @importFrom stats approx
cb_2014_nga <- function(M, T = 1000, Rrup, Rjb, Rx, W = 999, Ztor = 999, Zbot, dip, lambda, Fhw, Vs30, Z25 = 999, Zhyp = 999,
                        region) {

  period = c(0.010, 0.020, 0.030, 0.050, 0.075, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.0, 1.5,
             2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 0, -1)

  # Set initial A1100 value to 999. A1100: median estimated value of PGA on rock with Vs30 = 1100m/s
  A1100 = 999

  # Style of faulting
  Frv = (lambda > 30 & lambda < 150)

  Fnm = (lambda > -150 & lambda < -30)

  # if Ztor is unknown
  if (Ztor == 999)
    ifelse (Frv == 1, Ztor <- max(2.704 - 1.226 * max(M - 5.849, 0), 0)^2,
            Ztor <- max(2.673 - 1.136 * max(M - 4.970, 0), 0)^2)

  # if W is unknown
  if (W == 999) {

    ifelse(Frv == 1, Ztori <- max(2.704 - 1.226 * max(M - 5.849, 0), 0)^2,
           Ztori <- max(2.673 - 1.136 * max(M - 4.970, 0), 0)^2)

    W = min(sqrt(10^((M - 4.07) / 0.98)), (Zbot - Ztori) / sin(pi / 180 * dip))

    Zhyp = 9

  }

  # if Zhyp is unknown
  if (Zhyp == 999 & W != 999) {

    ifelse (M < 6.75, fdZM <- -4.317 + 0.984 * M, fdZM <- 2.325)

    ifelse (dip <= 40, fdZD <- 0.0445 * (dip - 40), fdZD <- 0)

    ifelse (Frv == 1, Ztori <- max(2.704 - 1.226 * max(M - 5.849, 0), 0)^2,
            Ztori <- max(2.673 - 1.136 * max(M - 4.970, 0), 0)^2)

    # The depth to the bottom of the rupture plane
    Zbor = Ztori + W * sin(pi / 180 * dip)

    d_Z = exp(min(fdZM + fdZD, log(0.9 * (Zbor - Ztori))))

    Zhyp = d_Z + Ztori

  }

  # Compute Sa and sigma with pre-defined period
  if (length(T) == 1 & T[1] == 1000) {

    Sa = rep(0, length(period))

    Sigma = rep(0, length(period))

    Phi = rep(0, length(period))

    Tau = rep(0, length(period))

    for (ip in 1:length(period)) {

      res <- cb_2014_subroutine(M, ip, Rrup, Rjb, Rx, W, Ztor, Zbot, dip, Fhw, Vs30, Z25,
                                Zhyp, lambda, Frv, Fnm, region, A1100)

      Sa[ip] <- res[1]

      Sigma[ip] <- res[2]

      Phi[ip] <- res[3]

      Tau[ip] <- res[4]

      res_PGA <- cb_2014_subroutine(M, 22, Rrup, Rjb, Rx, W, Ztor, Zbot, dip, Fhw, Vs30, Z25,
                                    Zhyp, lambda, Frv, Fnm, region, A1100)

      if (Sa[ip] < res_PGA[1] & period[ip] < 0.25 & period[ip] > 0)
        Sa[ip] <- res_PGA[1]

    }

    T <- period

  } else {

    Sa = rep(0, length(T))

    Sigma = rep(0, length(T))

    Phi = rep(0, length(T))

    Tau = rep(0, length(T))

    for (i in 1:length(T)) {

      Ti = T[i]

      if (min(abs(period - Ti)) > 0.0001) {

        # user defined period needs interpolation
        T_low = max(period(period < Ti))

        T_high = min(period(period > Ti))

        ip_low  = which(period == T_low)

        ip_high = which(period == T_high)

        res_low <- cb_2014_subroutine(M, ip_low, Rrup, Rjb, Rx, W, Ztor, Zbot, dip, Fhw, Vs30,
                                      Z25, Zhyp, lambda, Frv, Fnm, region, A1100)

        res_high <- cb_2014_subroutine(M, ip_high, Rrup, Rjb, Rx, W, Ztor, Zbot, dip, Fhw, Vs30,
                                       Z25, Zhyp, lambda, Frv, Fnm, region, A1100)

        res_PGA <- cb_2014_subroutine(M, 22, Rrup, Rjb, Rx, W, Ztor, Zbot, dip, Fhw, Vs30, Z25, Zhyp, lambda,
                                      Frv, Fnm, region, A1100)

        Sa[i] = exp(approx(x = c(log(T_low), log(T_high)), y = c(log(res_low[1]), log(res_high[1])),
                           xout = log(Ti), rule = 2)$y)

        Sigma[i] = exp(approx(x = c(log(T_low), log(T_high)), y = c(log(res_low[2]), log(res_high[2])),
                              xout = log(Ti), rule = 2)$y)

        Phi[i] = exp(approx(x = c(log(T_low), log(T_high)), y = c(log(res_low[3]), log(res_high[3])),
                            xout = log(Ti), rule = 2)$y)

        Tau[i] = exp(approx(x = c(log(T_low), log(T_high)), y = c(log(res_low[4]), log(res_high[4])),
                            xout = log(Ti), rule = 2)$y)

        if (Sa[i] < res_PGA[1] & Ti < 0.25 & Ti > 0)
          Sa[i] = res_PGA[1]

      } else {

        ip_T = which(abs(period- Ti) < 0.0001)

        res <- cb_2014_subroutine(M, ip_T, Rrup, Rjb, Rx, W, Ztor,Zbot,dip, Fhw, Vs30, Z25,
                                  Zhyp, lambda, Frv, Fnm, region, A1100)

        res_PGA <- cb_2014_subroutine(M, 22, Rrup, Rjb, Rx, W, Ztor, Zbot, dip, Fhw, Vs30, Z25,
                                      Zhyp, lambda, Frv, Fnm, region, A1100)

        Sa[i] <- res[1]

        Sigma[i] <- res[2]

        Phi[i] <- res[3]

        Tau[i] <- res[4]

        if (Sa[i] < res_PGA[1] & Ti < 0.25 & Ti > 0)
          Sa[i] <- res_PGA[1]

      }
    }
  }

  res <- list()

  res$med <- Sa

  res$sigma <- Sigma

  res$phi <- Phi

  res$tau <- Tau

  res$period <- T

  return(res)
}

#' The subroutine of GMPE of CB 2014
#'
#' This is a subroutine of CB 2014
#' @param M Moment magnitude, a numeric value
#' @param ip The index of working period in the per-defined periods: (0.010, 0.020, 0.030, 0.050, 0.075,
#' 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 0, -1)
#' @param Rrup Closest distance (km) to the ruptured plane
#' @param Rjb Joyner-Boore distance (km); closest distance (km) to surface
#' projection of rupture plane
#' @param Rx Site coordinate (km) measured perpendicular to the fault strike
#' from the fault line with down-dip direction to be positive
#' @param W Down-dip rupture width (km). 999 if unknown
#' @param Ztor Depth(km) to the top of ruptured plane. 999 if unknown
#' @param Zbot Depth(km) to the bottom of the seismogenic crust. It is needed only when W is unknown.
#' @param dip Fault dip angle (in degree)
#' @param Fhw Falg for hanging wall sites. 1 for sites on the hanging wall side of the fault, 0 otherwise.
#' @param Vs30 Shear wave velocity averaged over top 30 m (in m/s)
#' @param Z25in Basin depth (km): depth from the groundsurface to the 2.5km/s shear-wave horizon.
#' 999 if unknown.
#' @param Zhyp Hypocentral depth of the earthquake measured from sea level. 999 if unknown.
#' @param lambda Rake angle (in degree)
#' @param Frv Flag of reverse fault. 1 for reverse fault, 0 otherwise.
#' @param Fnm Flag of normal fault. 1 for normal fault, 0 otherwise.
#' @param region Region indicator: 0 for global (incl. Taiwan); 1 for California; 2 for Japan;
#' 3 for China or Turkey; 4 for Italy
#' @param A1100 Median estimated value of PGA on rock with Vs30 = 1100m/s. 999 if unknown
#' @return An array of four values: (1) median spectral acceleration prediction (in g); (2)
#' logarithmic standard deviation of spectral acceleration prediction; (3) logarithmic
#' standard deviation of within event residuals; and (4) logarithmic standard deviation of between
#' event residuals at the period of ip.
#' @references Campbell, K. W., and Bozorgnia, Y. (2014). NGA-West2 Ground Motion Model
#' for the Average Horizontal Components of PGA, PGV, and 5% Damped Linear
#' Acceleration Response Spectra. Earthquake Spectra, 30(3), 1087-1115.
#' @export
#' @importFrom stats approx
cb_2014_subroutine <- function(M, ip, Rrup, Rjb, Rx, W, Ztor, Zbot, dip, Fhw, Vs30, Z25in,
                        Zhyp, lambda, Frv, Fnm, region, A1100 = 999) {

  # coefficients
  c0 = c(-4.365,	-4.348,	-4.024,	-3.479,	-3.293,	-3.666,	-4.866,	-5.411,	-5.962,	-6.403,	-7.566,	-8.379,
         -9.841,	-11.011,	-12.469,	-12.969,	-13.306,	-14.02,	-14.558,	-15.509,	-15.975,	-4.416,
         -2.895)

  c1 = c(0.977,	0.976,	0.931,	0.887,	0.902,	0.993,	1.267,	1.366,	1.458,	1.528,	1.739,	1.872,
         2.021,	2.180,	2.270,	2.271,	2.150,	2.132,	2.116,	2.223,	2.132,	0.984,	1.510)

  c2 = c(0.533,	0.549,	0.628,	0.674,	0.726,	0.698,	0.510,	0.447,	0.274,	0.193,	-0.020,	-0.121,
         -0.042,	-0.069,	0.047,	0.149,	0.368,	0.726,	1.027,	0.169,	0.367,	0.537,	0.270)

  c3 = c(-1.485,	-1.488,	-1.494,	-1.388,	-1.469,	-1.572,	-1.669,	-1.750,	-1.711,	-1.770,	-1.594,	-1.577,
         -1.757,	-1.707,	-1.621,	-1.512,	-1.315,	-1.506,	-1.721,	-0.756,	-0.800,	-1.499,	-1.299)

  c4 = c(-0.499,	-0.501,	-0.517,	-0.615,	-0.596,	-0.536,	-0.490,	-0.451,	-0.404,	-0.321,	-0.426,	-0.440,
         -0.443,	-0.527,	-0.630,	-0.768,	-0.890,	-0.885,	-0.878,	-1.077,	-1.282,	-0.496,	-0.453)

  c5 = c(-2.773,	-2.772,	-2.782,	-2.791,	-2.745,	-2.633,	-2.458,	-2.421,	-2.392,	-2.376,	-2.303,	-2.296,
         -2.232,	-2.158,	-2.063,	-2.104,	-2.051,	-1.986,	-2.021,	-2.179,	-2.244,	-2.773,	-2.466)

  # c6	= c(0.248,	0.247,	0.246,	0.240,	0.227,	0.210,	0.183,	0.182,	0.189,	0.195,	0.185,	0.186,
  #        0.186,	0.169,	0.158,	0.158,	0.148,	0.135,	0.140,	0.178,	0.194,	0.248,	0.204) # these coefficients are in the PEER report

  c6 = c(0.248,	0.247,	0.246,	0.240,	0.227,	0.210,	0.183,	0.182,	0.189,	0.195,	0.185,	0.186,
         0.186,	0.169,	0.158,	0.158,	0.148,	0.135,	0.135,	0.165,	0.180,	0.248,	0.204) # these coefficients are in the Earthquake Spectra paper

  c7 = c(6.753,	6.502,	6.291,	6.317,	6.861,	7.294,	8.031,	8.385,	7.534,	6.990,	7.012,	6.902,
         5.522,	5.650,	5.795,	6.632,	6.759,	7.978,	8.538,	8.468,	6.564,	6.768,	5.837)

  c8 = c(0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0)

  c9 = c(-0.214,	-0.208,	-0.213,	-0.244,	-0.266,	-0.229,	-0.211,	-0.163,	-0.150,	-0.131,	-0.159,	-0.153,
         -0.090,	-0.105,	-0.058,	-0.028,	0,	0,	0,	0,	0,	-0.212,	-0.168)

  c10	= c(0.720,	0.730,	0.759,	0.826,	0.815,	0.831,	0.749,	0.764,	0.716,	0.737,	0.738,	0.718,
          0.795,	0.556,	0.480,	0.401,	0.206,	0.105,	0,	0,	0,	0.720,	0.305)

  c11	= c(1.094,	1.149,	1.290,	1.449,	1.535,	1.615,	1.877,	2.069,	2.205,	2.306,	2.398,	2.355,
          1.995,	1.447,	0.330,	-0.514,	-0.848,	-0.793,	-0.748,	-0.664,	-0.576,	1.090,	1.713)

  c12	= c(2.191,	2.189,	2.164,	2.138,	2.446,	2.969,	3.544,	3.707,	3.343,	3.334,	3.544,	3.016,
          2.616,	2.470,	2.108,	1.327,	0.601,	0.568,	0.356,	0.075,	-0.027,	2.186,	2.602)

  c13 = c(1.416,	1.453,	1.476,	1.549,	1.772,	1.916,	2.161,	2.465,	2.766,	3.011,	3.203,	3.333,
          3.054,	2.562,	1.453,	0.657,	0.367,	0.306,	0.268,	0.374,	0.297,	1.420,	2.457)

  c14	= c(-0.0070,	-0.0167,	-0.0422,	-0.0663,	-0.0794,	-0.0294,	0.0642,	0.0968,	0.1441,	0.1597,
          0.1410,	0.1474,	0.1764,	0.2593,	0.2881,	0.3112,	0.3478,	0.3747,	0.3382,	0.3754,	0.3506,	-0.0064,
          0.1060)

  c15 = c(-0.207,	-0.199,	-0.202,	-0.339,	-0.404,	-0.416,	-0.407,	-0.311,	-0.172,	-0.084,	0.085,	0.233,
          0.411,	0.479,	0.566,	0.562,	0.534,	0.522,	0.477,	0.321,	0.174,	-0.202,	0.332)

  c16	= c(0.390,	0.387,	0.378,	0.295,	0.322,	0.384,	0.417,	0.404,	0.466,	0.528,	0.540,	0.638,
          0.776,	0.771,	0.748,	0.763,	0.686,	0.691,	0.670,	0.757,	0.621,	0.393,	0.585)

  c17	= c(0.0981,	0.1009,	0.1095,	0.1226,	0.1165,	0.0998,	0.0760,	0.0571,	0.0437,	0.0323,	0.0209,	0.0092,
          -0.0082,	-0.0131,	-0.0187,	-0.0258,	-0.0311,	-0.0413,	-0.0281,	-0.0205,	0.0009,	0.0977,
          0.0517)

  c18	= c(0.0334,	0.0327,	0.0331,	0.0270,	0.0288,	0.0325,	0.0388,	0.0437,	0.0463,	0.0508,	0.0432,	0.0405,
          0.0420,	0.0426,	0.0380,	0.0252,	0.0236,	0.0102,	0.0034,	0.0050,	0.0099,	0.0333,	0.0327)

  c19	= c(0.00755,	0.00759,	0.00790,	0.00803,	0.00811,	0.00744,	0.00716,	0.00688,	0.00556,	0.00458,
          0.00401,	0.00388,	0.00420,	0.00409,	0.00424,	0.00448,	0.00345,	0.00603,	0.00805,	0.00280,
          0.00458,	0.00757,	0.00613)

  c20	= c(-0.0055,	-0.0055,	-0.0057,	-0.0063,	-0.0070,	-0.0073,	-0.0069,	-0.0060,	-0.0055,	-0.0049,
          -0.0037,	-0.0027,	-0.0016,	-0.0006,	0,	0,	0,	0,	0,	0,	0,	-0.0055,	-0.0017)

  Dc20 = c(0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0)

  Dc20_JI = c(-0.0035,	-0.0035,	-0.0034,	-0.0037,	-0.0037,	-0.0034,	-0.0030,	-0.0031,	-0.0033,
              -0.0035,	-0.0034,	-0.0034,	-0.0032,	-0.0030,	-0.0019,	-0.0005,	0,	0,	0,	0,	0,
              -0.0035,	-0.0006)

  Dc20_CH = c(0.0036,	0.0036,	0.0037,	0.0040,	0.0039,	0.0042,	0.0042,	0.0041,	0.0036,	0.0031,	0.0028,
              0.0025,	0.0016,	0.0006,	0,	0,	0,	0,	0,	0,	0,	0.0036,	0.0017)

  a2 = c(0.168,	0.166,	0.167,	0.173,	0.198,	0.174,	0.198,	0.204,	0.185,	0.164,	0.160,	0.184,
         0.216,	0.596,	0.596,	0.596,	0.596,	0.596,	0.596,	0.596,	0.596,	0.167,	0.596)

  h1 = c(0.242,	0.244,	0.246,	0.251,	0.260,	0.259,	0.254,	0.237,	0.206,	0.210,	0.226,	0.217,
         0.154,	0.117,	0.117,	0.117,	0.117,	0.117,	0.117,	0.117,	0.117,	0.241,	0.117)

  h2 = c(1.471,	1.467,	1.467,	1.449,	1.435,	1.449,	1.461,	1.484,	1.581,	1.586,	1.544,	1.554,
         1.626,	1.616,	1.616,	1.616,	1.616,	1.616,	1.616,	1.616,	1.616,	1.474,	1.616)

  h3 = c(-0.714,	-0.711,	-0.713,	-0.701,	-0.695,	-0.708,	-0.715,	-0.721,	-0.787,	-0.795,	-0.770,	-0.770,
         -0.780,	-0.733,	-0.733,	-0.733,	-0.733,	-0.733,	-0.733,	-0.733,	-0.733,	-0.715,	-0.733)

  h4 = c(1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,
         1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000)

  h5 = c(-0.336,	-0.339,	-0.338,	-0.338,	-0.347,	-0.391,	-0.449,	-0.393,	-0.339,	-0.447,	-0.525,	-0.407,
         -0.371,	-0.128,	-0.128,	-0.128,	-0.128,	-0.128,	-0.128,	-0.128,	-0.128,	-0.337,	-0.128)

  h6 = c(-0.270,	-0.263,	-0.259,	-0.263,	-0.219,	-0.201,	-0.099,	-0.198,	-0.210,	-0.121,	-0.086,	-0.281,
         -0.285,	-0.756,	-0.756,	-0.756,	-0.756,	-0.756,	-0.756,	-0.756,	-0.756,	-0.270,	-0.756)

  k1 = c(865, 865, 908, 1054, 1086, 1032, 878, 748, 654, 587, 503, 457, 410, 400, 400, 400, 400, 400, 400,
         400, 400, 865, 400)

  k2 = c(-1.186,	-1.219,	-1.273,	-1.346,	-1.471,	-1.624,	-1.931,	-2.188,	-2.381,	-2.518,	-2.657,	-2.669,
         -2.401,	-1.955,	-1.025,	-0.299,	0.000,	0.000,	0.000,	0.000,	0.000,	-1.186,	-1.955)

  k3 = c(1.839,	1.840,	1.841,	1.843,	1.845,	1.847,	1.852,	1.856,	1.861,	1.865,	1.874,	1.883,
         1.906,	1.929,	1.974,	2.019,	2.110,	2.200,	2.291,	2.517,	2.744,	1.839,	1.929)

  c	= 1.88

  n	= 1.18

  f1 = c(0.734,	0.738,	0.747,	0.777,	0.782,	0.769,	0.769,	0.761,	0.744,	0.727,	0.690,	0.663,
         0.606,	0.579,	0.541,	0.529,	0.527,	0.521,	0.502,	0.457,	0.441,	0.734,	0.655)

  f2 = c(0.492,	0.496,	0.503,	0.520,	0.535,	0.543,	0.543,	0.552,	0.545,	0.568,	0.593,	0.611,
         0.633,	0.628,	0.603,	0.588,	0.578,	0.559,	0.551,	0.546,	0.543,	0.492,	0.494)

  t1 = c(0.404,	0.417,	0.446,	0.508,	0.504,	0.445,	0.382,	0.339,	0.340,	0.340,	0.356,	0.379,
         0.430,	0.470,	0.497,	0.499,	0.500,	0.543,	0.534,	0.523,	0.466,	0.409,	0.317)

  t2 = c(0.325,	0.326,	0.344,	0.377,	0.418,	0.426,	0.387,	0.338,	0.316,	0.300,	0.264,	0.263,
         0.326,	0.353,	0.399,	0.400,	0.417,	0.393,	0.421,	0.438,	0.438,	0.322,	0.297)

  flnAF	= c(0.300,	0.300,	0.300,	0.300,	0.300,	0.300,	0.300,	0.300,	0.300,	0.300,	0.300,	0.300,
            0.300,	0.300,	0.300,	0.300,	0.300,	0.300,	0.300,	0.300,	0.300,	0.300,	0.300)

  rlnPGA_lnY = c(1.000, 0.998, 0.986, 0.938, 0.887, 0.870, 0.876,0.870, 0.850, 0.819, 0.743, 0.684, 0.562,
                 0.467, 0.364, 0.298, 0.234, 0.202, 0.184, 0.176, 0.154, 1.000, 0.684)

  # Adjustment factor based on region
  if (region == 2) {

    Dc20 = Dc20_JI

  } else if (region == 4) {

    Dc20 = Dc20_JI

  } else if (region == 3) {

    Dc20 = Dc20_CH

  }

  # if region is in Japan
  Sj = (region == 2)

  # if Z2.5 is unknown
  if (Z25in == 999) {

    if (region != 2) {

      Z25 = exp(7.089 - 1.144 * log(Vs30))

      Z25A = exp(7.089 - 1.144 * log(1100))

    } else if (region == 2) {

      Z25 = exp(5.359 - 1.102 * log(Vs30))

      Z25A = exp(5.359 - 1.102 * log(1100))

    }

  } else {

    # Assign Z2.5 from user input into Z25 and calc Z25A for Vs30=1100m/s
    if (region != 2) {

      Z25 = Z25in

      Z25A = exp(7.089 - 1.144 * log(1100))

    } else if (region == 2) {

      Z25 = Z25in

      Z25A = exp(5.359 - 1.102 * log(1100))

    }
  }

  # Magnitude dependence
  if (M <= 4.5) {

    fmag = c0[ip] + c1[ip] * M

  } else if (M <= 5.5) {

    fmag = c0[ip] + c1[ip] * M + c2[ip] * (M - 4.5)

  } else if (M <= 6.5) {

    fmag = c0[ip] + c1[ip] * M + c2[ip] * (M - 4.5) + c3[ip] * (M - 5.5)

  } else {

    fmag = c0[ip] + c1[ip] * M + c2[ip] * (M - 4.5) + c3[ip] * (M - 5.5) + c4[ip] * (M - 6.5)

  }

  # Geometric attenuation term
  fdis = (c5[ip] + c6[ip] * M) * log(sqrt(Rrup^2 + c7[ip]^2))

  # Style of faulting
  if (M <= 4.5) {

    F_fltm = 0

  } else if (M <= 5.5) {

    F_fltm = M - 4.5

  } else {

    F_fltm = 1

  }

  fflt = ((c8[ip] * Frv) + (c9[ip] * Fnm)) * F_fltm

  # Hanging-wall effects

  R1 = W * cos(pi / 180 * dip) # W - downdip width

  R2 = 62 * M - 350

  f1_Rx = h1[ip] + h2[ip] * (Rx / R1) + h3[ip] * (Rx / R1)^2

  f2_Rx = h4[ip] + h5[ip] * ((Rx - R1) / (R2 - R1)) + h6[ip] * ((Rx - R1) / (R2 - R1))^2

  if (Fhw == 0) {

    f_hngRx = 0

  } else if (Rx < R1 & Fhw == 1) {

    f_hngRx = f1_Rx

  } else if (Rx >= R1 & Fhw == 1) {

    f_hngRx = max(f2_Rx, 0)

  }

  ifelse (Rrup == 0, f_hngRup <- 1, f_hngRup <- (Rrup - Rjb) / Rrup)

  if (M <= 5.5) {

    f_hngM = 0

  } else if (M <= 6.5) {

    f_hngM = (M - 5.5) * (1 + a2[ip] * (M - 6.5))

  } else {

    f_hngM = 1 + a2[ip] * (M - 6.5)

  }

  ifelse (Ztor <= 16.66, f_hngZ <- 1 - 0.06 * Ztor, f_hngZ <- 0)

  f_hngdelta = (90 - dip) / 45

  fhng = c10[ip] * f_hngRx * f_hngRup * f_hngM * f_hngZ * f_hngdelta

  # Site conditions
  if (Vs30 <= k1[ip]) {

    if (A1100 == 999)
      A1100 = cb_2014_nga(M, 0, Rrup, Rjb, Rx, W, Ztor, Zbot, dip, lambda, Fhw, 1100, Z25A, Zhyp, region)$med

    f_siteG = c11[ip] * log(Vs30 / k1[ip]) + k2[ip] * (log(A1100 + c * (Vs30 / k1[ip])^n) - log(A1100 + c))

  } else if (Vs30 > k1[ip]) {

    f_siteG = (c11[ip] + k2[ip] * n) * log(Vs30 / k1[ip])

  }

  ifelse (Vs30 <= 200, f_siteJ <- (c12[ip] + k2[ip] * n) * (log(Vs30 / k1[ip]) - log(200 / k1[ip])) * Sj,
          f_siteJ <- (c13[ip] + k2[ip] * n) * log(Vs30 / k1[ip]) * Sj)

  fsite = f_siteG + f_siteJ

  # Basin Response Term - Sediment effects
  if (Z25 <= 1) {

    fsed = (c14[ip] + c15[ip] * Sj) * (Z25 - 1)

  } else if (Z25 <= 3) {

    fsed = 0

  } else if (Z25 > 3) {

    fsed = c16[ip] * k3[ip] * exp(-0.75) * (1 - exp(-0.25 * (Z25 - 3)))

  }

  # Hypocenteral Depth term
  if (Zhyp <= 7) {

    f_hypH = 0

  } else if (Zhyp <= 20) {

    f_hypH = Zhyp - 7

  } else {

    f_hypH = 13

  }

  if (M <= 5.5) {

    f_hypM = c17[ip]

  } else if (M <= 6.5) {

    f_hypM = c17[ip] + (c18[ip] - c17[ip]) * (M - 5.5)

  } else {

    f_hypM = c18[ip]

  }

  fhyp = f_hypH * f_hypM

  # Fault Dip term
  if (M <= 4.5) {

    f_dip = c19[ip] * dip

  } else if (M <= 5.5) {

    f_dip = c19[ip] * (5.5 - M) * dip

  } else {

    f_dip = 0

  }

  # Anelastic Attenuation Term
  ifelse (Rrup > 80, f_atn <- (c20[ip] + Dc20[ip]) * (Rrup - 80),
          f_atn <- 0)

  # Median value
  Sa = exp(fmag + fdis + fflt + fhng + fsite + fsed + fhyp + f_dip + f_atn)

  # Standard deviation computations
  if (M <= 4.5) {

    tau_lny = t1[ip]

    tau_lnPGA = t1[22] # PGA

    phi_lny = f1[ip]

    phi_lnPGA = f1[22]

  } else if (M < 5.5) {

    tau_lny = t2[ip] + (t1[ip] - t2[ip]) * (5.5 - M)

    tau_lnPGA = t2[22] + (t1[22] - t2[22]) * (5.5 - M) #  ip = PGA

    phi_lny = f2[ip] + (f1[ip] - f2[ip]) * (5.5 - M)

    phi_lnPGA = f2[22] + (f1[22] - f2[22]) * (5.5 - M)

  } else {

    tau_lny = t2[ip]

    tau_lnPGA = t2[22] # PGA

    phi_lny = f2[ip]

    phi_lnPGA = f2[22]

  }

  tau_lnyB = tau_lny

  tau_lnPGAB = tau_lnPGA

  phi_lnyB = sqrt(phi_lny^2 - flnAF[ip]^2)

  phi_lnPGAB = sqrt(phi_lnPGA^2 - flnAF[ip]^2)

  ifelse (Vs30 < k1[ip], alpha <- k2[ip] * A1100 * ((A1100 + c * (Vs30 / k1[ip])^n)^(-1) - (A1100 + c)^(-1)),
          alpha <- 0)

  tau = sqrt(tau_lnyB^2 + alpha^2 * tau_lnPGAB^2 + 2 * alpha * rlnPGA_lnY[ip] * tau_lnyB * tau_lnPGAB)

  phi = sqrt(phi_lnyB^2 + flnAF[ip]^2 + alpha^2 * phi_lnPGAB^2 +
               2 * alpha * rlnPGA_lnY[ip] * phi_lnyB * phi_lnPGAB)

  Sigma = sqrt(tau^2 + phi^2)

  return(c(Sa, Sigma, phi, tau))
}


















