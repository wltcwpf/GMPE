#' The GMPE of NGA Central and Eastern North America (CENA) Model
#'
#' This function calculates Vs30-specific ground motion median values and standard deviations
#'
#' @param M Moment magnitude, a numeric value. The applicable range is between 4 and 8.2
#' @param T Period (sec); Use Period = -1 for PGV computation.
#' Use 1000 (by default) for output the array of Sa with pre-defined CENA periods: (-1,
#' 0, 0.010, 0.020, 0.030, 0.040, 0.050, 0.075, 0.080, 0.100, 0.110, 0.112202, 0.114815,
#' 0.11749, 0.12, 0.125893, 0.13, 0.134896, 0.14, 0.15, 0.2, 0.25, 0.3, 0.4,
#' 0.5, 0.75, 0.8, 1.000, 1.500, 2.000, 3.000, 4.000, 5.000, 7.500, 10.000)
#' @param Rrup The closest distance (km) to the fault rupture. The applicable range is less than 1500 km.
#' @param Vs30 Shear wave velocity averaged over top 30 m (in m/s). The applicable range is between 200 and 3000 m/s.
#' @return A list of seven elements is returned: period - the input period of interest;
#' med - median spectral acceleration prediction (in g);
#' sigma - logarithmic standard deviation of spectral acceleration prediction;
#' med_ref - median spectral acceleration prediction on reference site condition;
#' Fs - the site amplification (in arithmetic); lnflin - the linear site amplification (in log);
#' lnfnl - the nonlinear site amplification (in log).
#' @examples nga_cena_gmm(M = 6, T = 1000, Rrup = 30, Vs30 = 350)
#' @references
#' Report PEER 2018/08 - Central and Eastern North America
#' Ground-Motion Characterization - NGA-East Final Report
#' Christine Goulet, Yousef Bozorgnia, Norman Abrahamson,
#' Nicolas Kuehn, Linda Al Atik, Robert Youngs, and Robert Graves
#' https://apps.peer.berkeley.edu/publications/peer_reports/reports_2018/NGA-East_1MainReport_FINAL.pdf
#'
#' Hashash YMA, Ilhan O, Harmon JA, et al. Nonlinear site amplification model for ergodic seismic hazard analysis in Central and Eastern North America. Earthquake Spectra. 2020;36(1):69-86. doi:10.1177/8755293019878193
#'
#' Stewart JP, Parker GA, Atkinson GM, Boore DM, Hashash YMA, Silva WJ. Ergodic site amplification model for central and eastern North America. Earthquake Spectra. 2020;36(1):42-68. doi:10.1177/8755293019878185
#' @export
#' @importFrom stats approx
#' @importFrom pracma interp2
nga_cena_gmm <- function(M, T = 1000, Rrup, Vs30){

  if (M < 4 | M > 8.2) {
    stop('Your input M is outside of the applicable range (between 4 and 8.2). Please adjust!')
  }

  if (Rrup > 1500 | Rrup < 1) {
    stop('Your input Rrup is outside of the applicable range (between 1 km and 1500 km). Please adjust!')
  }

  if (Vs30 < 200 | Vs30 > 3000) {
    stop('Your input Vs30 is outside of the applicable range (between 200 m/s and 3000 m/s). Please adjust!')
  }

  if (length(T) == 1 & T[1] == 1000) {
    T <- c(-1, 0, 0.010, 0.020, 0.030, 0.040, 0.050, 0.075, 0.080, 0.100, 0.110, 0.112202, 0.114815,
           0.11749, 0.12, 0.125893, 0.13, 0.134896, 0.14, 0.15, 0.2, 0.25, 0.3, 0.4,
           0.5, 0.75, 0.8, 1.000, 1.500, 2.000, 3.000, 4.000, 5.000, 7.500, 10.000)
  }

  med <- rep(NA, length(T))
  sigma <- rep(NA, length(T))
  med_ref <- rep(NA, length(T))
  Fs <- rep(NA, length(T))
  lnflin <- rep(NA, length(T))
  lnfnl <- rep(NA, length(T))

  PGAr <- nga_cena_gmm_ref(M = M, thisT = 0, Rrup = Rrup, SigmaType = 'Ergodic')$med

  for (i in 1:length(T)) {
    gmm_ref <- nga_cena_gmm_ref(M = M, thisT = T[i], Rrup = Rrup, SigmaType = 'Ergodic')
    med_ref[i] <- gmm_ref$med
    sigma[i] <- gmm_ref$sigma
    gmm_fs <- nga_cena_gmmsiteamp(thisT = T[i], Vs30 = Vs30, PGAr = PGAr,
                                  flin_coeffs = sea_2020_cena_flin_coeffs,
                                  fnl_coeffs = hea_2020_cena_fnl_coeffs)
    Fs[i] <- gmm_fs$Fs
    lnflin[i] <- gmm_fs$lnflin
    lnfnl[i] <- gmm_fs$lnfnl
    med[i] <- gmm_ref$med * gmm_fs$Fs
  }

  res <- list()
  res$period <- T
  res$med <- med
  res$sigma <- sigma
  res$med_ref <- med_ref
  res$Fs <- Fs
  res$lnflin <- lnflin
  res$lnfnl <- lnfnl

  return(res)
}




#' The subroutine of CENA ground motion prediction on reference condition
#'
#' This is a subroutine of CENA GMM for reference condition
#' @param M Moment magnitude, a numeric value
#' @param thisT One single input period. -1 is for PGV, 0 is for PGA, and 0.01-10 for PSA at these oscillator period range.
#' @param Rrup The closest distance (km) to the fault rupture
#' @param SigmaType Type of Sigma: Ergodic or SingleStation. By default, it is Ergodic.
#' @return a list of two elements, median GMM and standard deviation.
#' @references
#' Report PEER 2018/08 - Central and Eastern North America
#' Ground-Motion Characterization - NGA-East Final Report
#' Christine Goulet, Yousef Bozorgnia, Norman Abrahamson,
#' Nicolas Kuehn, Linda Al Atik, Robert Youngs, and Robert Graves
#' https://apps.peer.berkeley.edu/publications/peer_reports/reports_2018/NGA-East_1MainReport_FINAL.pdf
#' @export
nga_cena_gmm_ref <- function(M, thisT, Rrup, SigmaType = 'Ergodic'){

  # get Sigma
  Flist = c(100,50,40,33.333,25,20,13.333,10,6.667,5,4,3.333,2.5,2,1.333,1,0.667,0.5,0.333,0.25,0.2,0.133,0.1,0,-1)

  Tlist <- ifelse(Flist > 0, 1/Flist, Flist)

  if (SigmaType == 'Ergodic') {

    if(length(which(abs(Tlist - thisT) < 0.0001)) == 0) { # The user defined period requires interpolation
      idx <- which(Tlist < thisT)
      T_low = max(Tlist[idx])
      idx <- which(Tlist > thisT)
      T_high = min(Tlist[idx])
      ip_low  = which(Tlist == T_low)
      ip_high = which(Tlist == T_high)

      sigma_low <- getSigmaCompErg(ip_low, M)
      sigma_high <- getSigmaCompErg(ip_high, M)

      Sigma = approx(c(log(T_low), log(T_high)), c(sigma_low, sigma_high), log(thisT), rule = 2)$y
    } else {
      ip = which(abs((Tlist - thisT)) < 0.0001)
      Sigma <- getSigmaCompErg(ip, M)
    }

  } else if (SigmaType == 'SingleStation') {

    if(length(which(abs(Tlist - thisT) < 0.0001)) == 0) { # The user defined period requires interpolation
      idx <- which(Tlist < thisT)
      T_low = max(Tlist[idx])
      idx <- which(Tlist > thisT)
      T_high = min(Tlist[idx])
      ip_low  = which(Tlist == T_low)
      ip_high = which(Tlist == T_high)

      sigma_low <- getSigmaCompSS(ip_low, M)
      sigma_high <- getSigmaCompSS(ip_high, M)

      Sigma = approx(c(log(T_low), log(T_high)), c(sigma_low, sigma_high), log(thisT), rule = 2)$y
    } else {
      ip = which(abs((Tlist - thisT)) < 0.0001)
      Sigma <- getSigmaCompSS(ip, M)
    }

  }

  # Median GMM
  Flist = unique(MedianLnYrefDataList$MedianLnYrefFrequencyValueColumn)

  Tlist <- ifelse(Flist > 0, 1/Flist, Flist)

  if(length(which(abs(Tlist - thisT) < 0.0001)) == 0) { # The user defined period requires interpolation
    idx <- which(Tlist < thisT)
    T_low = max(Tlist[idx])
    idx <- which(Tlist > thisT)
    T_high = min(Tlist[idx])
    ip_low  = which(Tlist == T_low)
    ip_high = which(Tlist == T_high)

    LnY_low <- getLnY(ip_low, M, Rrup)
    LnY_high <- getLnY(ip_high, M, Rrup)

    LnY = approx(c(log(T_low), log(T_high)), c(LnY_low, LnY_high), log(thisT), rule = 2)$y
  } else {
    ip = which(abs((Tlist - thisT)) < 0.0001)
    LnY <- getLnY(ip, M, Rrup)
  }

  res <- list()
  res$med <- exp(LnY)
  res$sigma <- Sigma
  return(res)
}

getSigmaCompErg <- function(ip, M){

  WeightSigmaCompErgC <- 0.63
  WeightSigmaCompErgH <- 0.185
  WeightSigmaCompErgL <- 0.185
  WeightSigmaCompErgSum <- 1

  # Central Model
  SigCompErg1C <- c(0.7618,0.761,0.7646,0.7676,0.781,0.792,0.7881,0.7726,0.7597,0.7556,0.7523,0.7488,0.7412,
                    0.7335,0.7151,0.6992,0.6742,0.6562,0.6337,0.6212,0.6139,0.6053,0.6016,0.7618,0.7045)
  SigCompErg2C <- c(0.7471,0.7463,0.75,0.753,0.7667,0.7779,0.7739,0.7581,0.745,0.7408,0.7374,0.7338,0.7261,
                    0.7182,0.6995,0.6832,0.6576,0.6392,0.616,0.6032,0.5956,0.5867,0.5829,0.7471,0.6995)
  SigCompErg3C <- c(0.6819,0.6814,0.6856,0.6891,0.7044,0.7169,0.7133,0.6969,0.685,0.6844,0.6844,0.6842,0.6825,
                    0.6786,0.6661,0.6528,0.6299,0.6127,0.5885,0.5749,0.567,0.5576,0.5535,0.6819,0.6569)
  SigCompErg4C <- c(0.5809,0.5808,0.5861,0.5905,0.609,0.6239,0.6214,0.6039,0.5941,0.6001,0.6065,0.6129,0.6231,
                    0.6274,0.6283,0.6215,0.605,0.5907,0.5651,0.5509,0.5426,0.5327,0.5285,0.5809,0.5825)
  S1 = SigCompErg1C[ip]
  S2 = SigCompErg2C[ip]
  S3 = SigCompErg3C[ip]
  S4 = SigCompErg4C[ip]
  SigmaCompErgC = getSigmaCompErgi(M, S1, S2, S3, S4)

  # High Model
  SigCompErg1H <- c(0.8747,0.8737,0.8767,0.8791,0.8909,0.9011,0.899,0.8858,0.8725,0.8668,0.8621,0.8575,0.8485,
                 0.8402,0.8228,0.8096,0.7916,0.7807,0.7687,0.7627,0.7591,0.7545,0.7519,0.8747,0.8134)
  SigCompErg2H <- c(0.8635,0.8625,0.8655,0.8679,0.8799,0.8902,0.888,0.8747,0.8613,0.8555,0.8508,0.8461,0.837,
                    0.8286,0.811,0.7976,0.7792,0.7681,0.7558,0.7496,0.746,0.7412,0.7386,0.8635,0.809)
  SigCompErg3H <- c(0.7963,0.7958,0.7993,0.8022,0.8158,0.8275,0.8267,0.8137,0.8023,0.8004,0.7994,0.7983,0.7956,
                    0.7918,0.7814,0.7717,0.7568,0.7471,0.7343,0.7276,0.7238,0.7188,0.716,0.7963,0.7661)
  SigCompErg4H <- c(0.7076,0.7078,0.7121,0.7157,0.7315,0.7452,0.7457,0.7328,0.7233,0.7246,0.727,0.7295,0.734,
                    0.736,0.7361,0.7314,0.7203,0.7139,0.7025,0.6967,0.6932,0.6883,0.6855,0.7076,0.7163)
  S1 = SigCompErg1H[ip]
  S2 = SigCompErg2H[ip]
  S3 = SigCompErg3H[ip]
  S4 = SigCompErg4H[ip]
  SigmaCompErgH = getSigmaCompErgi(M, S1, S2, S3, S4)

  # Low Model
  SigCompErg1L <- c(0.6539,0.6533,0.6575,0.661,0.6758,0.6874,0.682,0.6645,0.6521,0.6496,0.6476,0.6452,0.6392,
                 0.6323,0.6145,0.5977,0.5703,0.5497,0.5229,0.5076,0.4983,0.4871,0.4826,0.6539,0.6034)
  SigCompErg2L <- c(0.6361,0.6355,0.6398,0.6434,0.6586,0.6705,0.665,0.6471,0.6343,0.6317,0.6296,0.6271,0.621,
                    0.6138,0.5954,0.578,0.5496,0.5282,0.5004,0.4844,0.4747,0.463,0.4583,0.6361,0.598)
  SigCompErg3L <- c(0.5736,0.5731,0.5779,0.5819,0.5987,0.6118,0.606,0.5867,0.5743,0.575,0.5761,0.5767,0.576,
                    0.5723,0.5586,0.543,0.5156,0.4945,0.4648,0.448,0.4379,0.4256,0.4207,0.5736,0.5587)
  SigCompErg4L <- c(0.4807,0.4805,0.4864,0.4913,0.5113,0.5264,0.52,0.4982,0.4876,0.4958,0.504,0.5119,0.5244,
                    0.5293,0.5299,0.5227,0.508,0.4926,0.4624,0.4446,0.4337,0.4205,0.4153,0.4807,0.4936)
  S1 = SigCompErg1L[ip]
  S2 = SigCompErg2L[ip]
  S3 = SigCompErg3L[ip]
  S4 = SigCompErg4L[ip]

  SigmaCompErgL = getSigmaCompErgi(M, S1, S2, S3, S4)
  SigmaCompErg = SigmaCompErgC * WeightSigmaCompErgC / WeightSigmaCompErgSum +
    SigmaCompErgH * WeightSigmaCompErgH / WeightSigmaCompErgSum +
    SigmaCompErgL * WeightSigmaCompErgL / WeightSigmaCompErgSum

  return(SigmaCompErg)
}

getSigmaCompErgi <- function(M, S1, S2, S3, S4){

  if(M <= 4.5 ){
    SigmaCompErgi = S1
  } else if(M <= 5.0 ){
    SigmaCompErgi = S1 + (S2 - S1) * (M - 4.5) / 0.5
  } else if(M <= 5.5 ){
    SigmaCompErgi = S2 + (S3 - S2) * (M - 5.0) / 0.5
  } else if(M <= 6.5 ){
    SigmaCompErgi = S3 + (S4 - S3) * (M - 5.5) / 1.0
  } else {
    SigmaCompErgi = S4
  }

  return(SigmaCompErgi)
}

getSigmaCompSS <- function(ip, M){

  WeightSigmaCompSSC <- 0.63
  WeightSigmaCompSSH <- 0.185
  WeightSigmaCompSSL <- 0.185
  WeightSigmaCompSSSum <- 1

  # Central Model
  SigCompSS1C <- c(0.7006,0.6998,0.6994,0.699,0.6981,0.6974,0.6955,0.6935,0.6898,0.6862,
                   0.6827,0.6794,0.673,0.6671,0.6541,0.6432,0.6262,0.6143,0.5998,0.5925,
                   0.5889,0.5858,0.5853,0.7006,0.6303)
  SigCompSS2C <- c(0.6846,0.6838,0.6833,0.6829,0.682,0.6813,0.6793,0.6773,0.6735,0.6698,
                   0.6662,0.6628,0.6563,0.6502,0.6369,0.6256,0.6082,0.5959,0.5811,0.5736,
                   0.5698,0.5666,0.5661,0.6846,0.6247)
  SigCompSS3C <- c(0.6127,0.6122,0.6119,0.6116,0.611,0.6105,0.6092,0.6079,0.6063,0.6066,
                   0.6068,0.6072,0.6075,0.606,0.5999,0.5921,0.5781,0.5674,0.5517,0.5438,
                   0.5397,0.5364,0.5357,0.6127,0.5765)
  SigCompSS4C <- c(0.4973,0.4974,0.4974,0.4974,0.4975,0.4975,0.4977,0.4978,0.5007,0.5092,
                   0.5171,0.5254,0.5399,0.5482,0.5577,0.5576,0.551,0.5436,0.5267,0.5181,
                   0.5138,0.5102,0.5095,0.4973,0.4896)
  S1 <- SigCompSS1C[ip]
  S2 <- SigCompSS2C[ip]
  S3 <- SigCompSS3C[ip]
  S4 <- SigCompSS4C[ip]
  SigmaCompSSC = getSigmaCompSSi(M, S1, S2, S3, S4)

  # High Model
  SigCompSS1H <- c(0.8193,0.8183,0.8178,0.8173,0.8162,0.8153,0.8129,0.8106,0.8061,0.8019,
                   0.798,0.7943,0.7876,0.7816,0.7697,0.7611,0.7502,0.7444,0.7391,0.7372,
                   0.7364,0.7358,0.7357,0.8193,0.748)
  SigCompSS2H <- c(0.8075,0.8064,0.8059,0.8054,0.8043,0.8034,0.801,0.7986,0.7941,0.7899,
                   0.7859,0.7822,0.7753,0.7693,0.7572,0.7484,0.7373,0.7313,0.7258,0.7237,
                   0.7228,0.7222,0.7221,0.8075,0.7431)
  SigCompSS3H <- c(0.7348,0.7343,0.734,0.7338,0.7333,0.7329,0.7317,0.7306,0.7292,0.7295,
                   0.7297,0.73,0.7305,0.7295,0.7256,0.7209,0.7137,0.7093,0.7033,0.7009,
                   0.6999,0.699,0.6989,0.7348,0.696)
  SigCompSS4H <- c(0.6377,0.6379,0.6381,0.6382,0.6384,0.6385,0.639,0.6395,0.6413,0.6453,
                   0.6492,0.6537,0.6623,0.6681,0.6762,0.6772,0.6743,0.6738,0.6698,0.6687,
                   0.6682,0.6679,0.6679,0.6377,0.6419)
  S1 <- SigCompSS1H[ip]
  S2 <- SigCompSS2H[ip]
  S3 <- SigCompSS3H[ip]
  S4 <- SigCompSS4H[ip]
  SigmaCompSSH <- getSigmaCompSSi(M, S1, S2, S3, S4)

  # Low Model
  SigCompSS1L <- c(0.5879,0.5873,0.587,0.5867,0.586,0.5855,0.584,0.5825,0.5795,0.5765,
                   0.5736,0.5707,0.565,0.5595,0.5467,0.5354,0.5171,0.5037,0.487,0.4784,
                   0.474,0.4702,0.4695,0.5879,0.522)
  SigCompSS2L <- c(0.5682,0.5676,0.5672,0.5669,0.5662,0.5657,0.5641,0.5626,0.5595,0.5564,
                   0.5533,0.5503,0.5444,0.5387,0.5253,0.5135,0.4944,0.4803,0.4628,0.4538,
                   0.4491,0.4451,0.4444,0.5682,0.5159)
  SigCompSS3L <- c(0.4981,0.4976,0.4973,0.497,0.4964,0.496,0.4946,0.4933,0.4917,0.4921,
                   0.4924,0.4928,0.4931,0.4913,0.4837,0.4741,0.4566,0.4433,0.4243,0.4148,
                   0.4101,0.4061,0.4053,0.4981,0.4703)
  SigCompSS4L <- c(0.3901,0.3899,0.3898,0.3897,0.3895,0.3893,0.3887,0.3882,0.3904,0.3994,
                   0.4079,0.4168,0.4326,0.4411,0.4504,0.4508,0.4482,0.4416,0.4226,0.4124,
                   0.4071,0.4025,0.4017,0.3901,0.3928)
  S1 <- SigCompSS1L[ip]
  S2 <- SigCompSS2L[ip]
  S3 <- SigCompSS3L[ip]
  S4 <- SigCompSS4L[ip]

  SigmaCompSSL <- getSigmaCompSSi(M, S1, S2, S3, S4)
  SigmaCompSS <- SigmaCompSSC * WeightSigmaCompSSC / WeightSigmaCompSSSum +
    SigmaCompSSH * WeightSigmaCompSSH / WeightSigmaCompSSSum +
    SigmaCompSSL * WeightSigmaCompSSL / WeightSigmaCompSSSum

  return(SigmaCompSS)
}

getSigmaCompSSi <- function(M, S1, S2, S3, S4){

  if(M <= 4.5 ){
    SigmaCompSSi = S1
  } else if(M <= 5.0 ){
    SigmaCompSSi = S1 + (S2 - S1) * (M - 4.5) / 0.5
  } else if(M <= 5.5 ){
    SigmaCompSSi = S2 + (S3 - S2) * (M - 5.0) / 0.5
  } else if(M <= 6.5 ){
    SigmaCompSSi = S3 + (S4 - S3) * (M - 5.5) / 1.0
  } else {
    SigmaCompSSi = S4
  }

  return(SigmaCompSSi)
}

getLnY <-function(ip, M, Rrup){

  idx_row <- c(1:11) + (ip - 1)*11

  x_Rup <- MedianLnYrefDataList$MedianLnYrefRrupRow

  y_M <- MedianLnYrefDataList$MedianLnYrefMomentColumn[idx_row]

  # linear interpolation for M and Rrup
  LnY = interp2(x = log(x_Rup[-1]), y = y_M, Z = MedianLnYrefDataList$MedianLnYrefDataTable[idx_row,-1],
                xp = log(Rrup), yp = M)

  return(LnY)
}


#' The subroutine of CENA site amplification function
#'
#' This is a subroutine of CENA GMM for site amplification
#' @param thisT One single input period. -1 is for PGV, 0 is for PGA, and 0.01-10 for PSA at these oscillator period range.
#' @param Vs30 Shear wave velocity averaged over top 30 m (in m/s)
#' @param PGAr The PGA on reference rock condition (Vs30 = 3000 m/s)
#' @param flin_coeffs The coefficient for linear site amplification from Stewart et al 2020 CENA Flin coefficients.
#' You can use the internal saved data object, sea_2020_cena_flin_coeffs
#' @param fnl_coeffs The coefficient for nonlinear site amplification from Hashash et al 2020 CENA Fnl coefficients.
#' You can use the internal saved data object, hea_2020_cena_fnl_coeffs
#' @references
#' Hashash YMA, Ilhan O, Harmon JA, et al. Nonlinear site amplification model for ergodic seismic hazard analysis in Central and Eastern North America. Earthquake Spectra. 2020;36(1):69-86. doi:10.1177/8755293019878193
#'
#' Stewart JP, Parker GA, Atkinson GM, Boore DM, Hashash YMA, Silva WJ. Ergodic site amplification model for central and eastern North America. Earthquake Spectra. 2020;36(1):42-68. doi:10.1177/8755293019878185
#' @export
nga_cena_gmmsiteamp <- function(thisT, Vs30, PGAr, flin_coeffs = sea_2020_cena_flin_coeffs,
                                fnl_coeffs = hea_2020_cena_fnl_coeffs){

  # linear site amplifcation
  # ----------------
  Tlist <- flin_coeffs$Period..s.
  w1 <- 0.767
  w2 <- 0.1
  Vw1 <- 600
  Vw2 <- 400

  if(length(which(abs(Tlist - thisT) < 0.0001)) == 0) { # The user defined period requires interpolation
    idx <- which(Tlist < thisT)
    T_low = max(Tlist[idx])
    idx <- which(Tlist > thisT)
    T_high = min(Tlist[idx])
    ip_low  = which(Tlist == T_low)
    ip_high = which(Tlist == T_high)

    c = approx(c(log(flin_coeffs$Period..s.[ip_low]), log(flin_coeffs$Period..s.[ip_high])),
               flin_coeffs$c[c(ip_low, ip_high)], log(thisT), rule = 2)$y
    Vref = approx(c(log(flin_coeffs$Period..s.[ip_low]), log(flin_coeffs$Period..s.[ip_high])),
                  flin_coeffs$Vref[c(ip_low, ip_high)], log(thisT), rule = 2)$y
    V1 = approx(c(log(flin_coeffs$Period..s.[ip_low]), log(flin_coeffs$Period..s.[ip_high])),
                flin_coeffs$V1[c(ip_low, ip_high)], log(thisT), rule = 2)$y
    V2 = approx(c(log(flin_coeffs$Period..s.[ip_low]), log(flin_coeffs$Period..s.[ip_high])),
                flin_coeffs$V2[c(ip_low, ip_high)], log(thisT), rule = 2)$y
    Vf = approx(c(log(flin_coeffs$Period..s.[ip_low]), log(flin_coeffs$Period..s.[ip_high])),
                flin_coeffs$Vf[c(ip_low, ip_high)], log(thisT), rule = 2)$y
    Vl = approx(c(log(flin_coeffs$Period..s.[ip_low]), log(flin_coeffs$Period..s.[ip_high])),
                flin_coeffs$Vl[c(ip_low, ip_high)], log(thisT), rule = 2)$y
    Vu = approx(c(log(flin_coeffs$Period..s.[ip_low]), log(flin_coeffs$Period..s.[ip_high])),
                flin_coeffs$Vu[c(ip_low, ip_high)], log(thisT), rule = 2)$y
    F760_imp = approx(c(log(flin_coeffs$Period..s.[ip_low]), log(flin_coeffs$Period..s.[ip_high])),
                      flin_coeffs$F760_imp[c(ip_low, ip_high)], log(thisT), rule = 2)$y
    F760_grad = approx(c(log(flin_coeffs$Period..s.[ip_low]), log(flin_coeffs$Period..s.[ip_high])),
                       flin_coeffs$F760_grad[c(ip_low, ip_high)], log(thisT), rule = 2)$y
  } else {
    ip = which(abs((Tlist - thisT)) < 0.0001)
    c = flin_coeffs$c[ip]
    Vref = flin_coeffs$Vref[ip]
    V1 = flin_coeffs$V1[ip]
    V2 = flin_coeffs$V2[ip]
    Vf = flin_coeffs$Vf[ip]
    Vl = flin_coeffs$Vl[ip]
    Vu = flin_coeffs$Vu[ip]
    F760_imp = flin_coeffs$F760_imp[ip]
    F760_grad = flin_coeffs$F760_grad[ip]
  }

  # Compute Impedance Weight
  if (Vs30 >= Vw1) {
    w_imp = w1
  } else if (Vs30 < Vw2) {
    w_imp = w2
  } else {
    w_imp = (w1 - w2)*log(Vs30/Vw2)/log(Vw1/Vw2) + w2
  }

  # Compute Gradient Weight
  w_grad = 1 - w_imp

  # F760 model
  F760 = w_imp*F760_imp + w_grad*F760_grad

  # VS30-Scaling Model
  if ((Vs30 >= Vl) & (Vs30 <= V1)) {
    Fv = c*log(V1/Vref)
  } else if ((Vs30 > V1) & (Vs30 <= V2)) {
    Fv = c*log(Vs30/Vref)
  } else if ((Vs30 > V2) & (Vs30 <= Vu)) {
    Fv = c*log(V2/Vref)
  } else if ((Vs30 > Vu) & (Vs30 <= 3000)) {
    Fv = c*log(V2/Vref) - (c*log(V2/Vref) + F760)*(log(Vs30/Vu)/log(3000/Vu))
  } else {
    stop('Error: Vs30 outside of usable range, between 200 m/s and 3000 m/s')
  }

  # Combine VS30-Scaling and F760 models:
  lnFlin = F760 + Fv


  # nonlinear site amplifcation
  # ----------------
  Vref <- 760

  Tlist = fnl_coeffs$period

  if(length(which(abs(Tlist - thisT) < 0.0001)) == 0) { # The user defined period requires interpolation
    idx <- which(Tlist < thisT)
    T_low = max(Tlist[idx])
    idx <- which(Tlist > thisT)
    T_high = min(Tlist[idx])
    ip_low  = which(Tlist == T_low)
    ip_high = which(Tlist == T_high)

    f3 = approx(c(log(fnl_coeffs$period[ip_low]), log(fnl_coeffs$period[ip_high])),
                fnl_coeffs$f3[c(ip_low, ip_high)], log(thisT), rule = 2)$y
    f4_original = approx(c(log(fnl_coeffs$period[ip_low]), log(fnl_coeffs$period[ip_high])),
                         fnl_coeffs$f4[c(ip_low, ip_high)], log(thisT), rule = 2)$y
    f4_modified = approx(c(log(fnl_coeffs$period[ip_low]), log(fnl_coeffs$period[ip_high])),
                         fnl_coeffs$f4_modified[c(ip_low, ip_high)], log(thisT), rule = 2)$y
    f5 = approx(c(log(fnl_coeffs$period[ip_low]), log(fnl_coeffs$period[ip_high])),
                fnl_coeffs$f5[c(ip_low, ip_high)], log(thisT), rule = 2)$y
    Vc = approx(c(log(fnl_coeffs$period[ip_low]), log(fnl_coeffs$period[ip_high])),
                fnl_coeffs$Vc[c(ip_low, ip_high)], log(thisT), rule = 2)$y
  } else {
    ip = which(abs((Tlist - thisT)) < 0.0001)
    f3 = fnl_coeffs$f3[ip]
    f4_original = fnl_coeffs$f4[ip]
    f4_modified = fnl_coeffs$f4_modified[ip]
    f5 = fnl_coeffs$f5[ip]
    Vc = fnl_coeffs$Vc[ip]
  }

  # Compute "Equally Weighted Nonlinear Amplification Model" f4 coefficient:
  f4 = 0.5*(f4_original + f4_modified)

  # Compute f2:
  f2 = f4*(exp(f5*(min(Vs30,Vref) - 360)) - exp(f5*(Vref - 360)))

  # Nonlinear Effects Model:
  if (Vs30 < Vc) {
    lnFnl = f2*log((PGAr+f3)/f3)
  } else {
    lnFnl = 0
  }

  # return results
  # ----------------
  res <- list()
  res$lnflin <- lnFlin
  res$lnfnl <- lnFnl
  res$Fs <- exp(lnFlin + lnFnl)
  return(res)
}
