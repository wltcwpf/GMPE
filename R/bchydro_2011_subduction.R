#' The GMPE of BCHydro 2011 for subduction earthquakes
#'
#' This function calculates the ground motion median values and standard deviations
#' @param M Moment magnitude, a numeric value
#' @param T Period (sec) (range should be between 0 s and 10 s).
#' Use 1000 (by default) for output the array of Sa with original BCHyrdo per-defined periods
#' @param Vs30 The time averaged shear wave velocity on the top 30 m.
#' @param Rrup Closest distance (km) to the ruptured plane
#' @param Rhypo Hypocentral distance(km)
#' @param Ftype Flag for interface event or intraslab events:
#' 0 for Interface - use rupture distance; 1 for Intraslab - use hypocentral distance.
#' @param faba Flag for BackArc sites: 0 for Non-Backarc sites; 1 for Backarc sites
#' @param depth The depth of hypocenter (km)
#' @param flag_deltaC1 Flag for Recommended deltaC1 Values:
#' 0 for Recommended deltaC1 Values NOT period dependent
#' (-0.5, 0.0, 0.5) as on the 06/2010 model;
#' 1 for Recommended deltaC1 Values Period dependent as on the 08/2011 model.
#' Note: here we implement only the CENTRAL deltaC1 values!
#' @return A list of four elements is returned: Sa - median spectral acceleration prediction (in g);
#' RockSa - median spectral acceleration prediction on Rock reference site (in g);
#' sigma - totla standard deviation (log); period - the corresponding oscillator periods
#' @examples bchydro_2011(M = 6, T = 1000, Vs30 = 450, Rrup = 85, Rhypo = 100,
#' Ftype = 0, faba = 0, depth = 25, flag_deltaC1 = 1)
#'
#' bchydro_2011(M = 7.5, T = c(0.01, 2.5), Vs30 = 450, Rrup = 85, Rhypo = 100,
#' Ftype = 1, faba = 1, depth = 65, flag_deltaC1 = 1)
#' @references BC Hydro subduction Ground motion model.
#' A report to BC Hydro Engineering No. E658 - Vol. 3 May 2009 by Abrahamson, et al.
#' @export
#' @importFrom stats approx
bchydro_2011 <- function(M, T = 1000, Vs30, Rrup, Rhypo, Ftype, faba, depth, flag_deltaC1 = 1){

  Vs30_rock = 1000
  Td = 0

  res = BCHydro_2011_subroutine(M, Vs30_rock, Td, Rrup, Rhypo, Ftype, faba, depth, flag_deltaC1, 0)
  PGA_rock <- res$Sa

  res <- BCHydro_2011_subroutine(M, Vs30, T, Rrup, Rhypo, Ftype, faba, depth, flag_deltaC1, PGA_rock)
  Sa <- res$Sa
  Sa_noamp <- res$Sa_noamp
  tsigma <- res$tsigma
  period1 <- res$period1

  bc <- list()
  bc$Sa <- Sa
  bc$RockSa <- Sa_noamp
  bc$sigma <- tsigma
  bc$period <- period1

  return(bc)

}


BCHydro_2011_subroutine <- function(M, Vs30, T, Rrup, Rhypo, Ftype, faba, depth, flag_deltaC1, PGA_rock){

  # The subroutine of BCHydro 2011 for subduction earthquakes

  # for the given period T, get the index for the constants
  period = c(0.0,0.02,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.75,
             1.00,1.5,2.00,2.5,3.00,4.0,5.0,6.0,7.5,10)

  Td = 0  # PGA case
  V_rock = get_BCHydro_constants(1,flag_deltaC1)  # return a list
  PGA_rock = exp(calc_val(M, Rrup, Rhypo, 0, 1000, Ftype, faba, depth, V_rock))

  nT = length(T)
  iflg = 0

  if(nT == 1){
    if(T == 1000){
      iflg = 1
      nperi = length(period)
      Sa = rep(0, nperi)
      Sa_noamp = rep(0, nperi)
      LSa_noamp <- rep(0, nperi)
      tsigma = rep(0, nperi)
      sigma = rep(0, nperi)
      tau = rep(0, nperi)
      period1 = seq(1, nperi)
      for(index in 1:nperi){
        # get constants for the given index value
        V = get_BCHydro_constants(index,flag_deltaC1)  # return a list
        if(period[index] == 0){  # PGA case
          Sa_noamp[index] = exp(calc_val(M, Rrup, Rhypo, PGA_rock, Vs30, Ftype, faba, depth, V))
          Sa[index] = Sa_noamp[index];
        }else{  # other spectral periods, implies amplification!
          Sa_noamp[index] = exp(calc_val2(M, Rrup, Rhypo, PGA_rock, Vs30, Ftype, faba, depth, V))
          specT = period[index]   # cambio!!!
          LSa_noamp[index] = log(Sa_noamp[index])
          AmpFac = BCHHR2Vs760(specT,LSa_noamp[index])
          # Sa[index]=Sa_noamp[index] + AmpFac + 6.89
          Sa[index] = exp(LSa_noamp[index]) * exp(AmpFac)
        }
        res <- BCHydro_sigma(V)
        tsigma[index] <- res$tsigma
        sigma[index] <- res$sigma
        tau[index] <- res$tau
      }
      period1 <- period
    }
  }

  if(iflg == 0){
    Sa <- rep(0, nT)
    tsigma = rep(0, nT)
    sigma <- rep(0, nT)
    tau = rep(0, nT)
    Sa_noamp <- rep(0, nT)
    LSa_noamp <- rep(0, nT)
    period1 = T
    for(it in 1:nT){
      Teach = T[it]
      if(Teach > period[length(period)]){
        Teach = period[length(period)]
        period1 = Teach
      }
      # interpolate between periods if neccesary
      if(length(which(period == Teach)) == 0){
        T_low <- max(period[which(period < Teach)])
        T_hi <- min(period[which(period > Teach)])
        V_low =  get_BCHydro_constants(max(which(period < Teach)), flag_deltaC1)  # return a list
        V_hi =  get_BCHydro_constants(min(which(period > Teach)), flag_deltaC1)   # return a list

        sa_low = exp(calc_val2(M, Rrup, Rhypo, PGA_rock, Vs30, Ftype, faba, depth, V_low))
        sa_hi = exp(calc_val2(M, Rrup, Rhypo, PGA_rock, Vs30, Ftype, faba, depth, V_hi))

        x = c(T_low, T_hi)
        Y_sa = c(sa_low, sa_hi)

        Sa_noamp[it] = approx(x = x, y = Y_sa, xout = Teach, rule = 2)$y
        LSa_noamp[it] = log(Sa_noamp[it])
        AmpFac = BCHHR2Vs760(T[it], LSa_noamp[it])
        # Sa[index] = Sa_noamp[index] + AmpFac + 6.89
        Sa[it] = exp(LSa_noamp[it]) * exp(AmpFac)

        res <- BCHydro_sigma(V_low)   ## sigma values do not change with period!!!!
        tsigma[it] <- res$tsigma
        sigma[it] <- res$sigma
        tau[it] <- res$tau

      }else{
        index = which(abs((period - Teach)) < 0.0001) # Identify the period
        # get constants for the given index value
        V = get_BCHydro_constants(index, flag_deltaC1)  # return a list
        if(period[index] == 0){  # PGA case
          Sa_noamp[it] = exp(calc_val(M, Rrup, Rhypo, PGA_rock, Vs30, Ftype, faba, depth, V))
          Sa[it] = Sa_noamp[it]
        }else{
          Sa_noamp[it] = exp(calc_val2(M, Rrup, Rhypo, PGA_rock, Vs30, Ftype, faba, depth, V))
          LSa_noamp[it]=log(Sa_noamp[it])
          AmpFac = BCHHR2Vs760(T[it], LSa_noamp[it])
          # Sa[index] = Sa_noamp[index] + AmpFac + 6.89
          Sa[it] = exp(LSa_noamp[it]) * exp(AmpFac)
        }
        res <- BCHydro_sigma(V)
        tsigma[it] <- res$tsigma
        sigma[it] <- res$sigma
        tau[it] <- res$tau
      }
    }
    rock <- PGA_rock
  }
  res <- list()
  res$Sa <- Sa
  res$Sa_noamp <- Sa_noamp
  res$tsigma <- tsigma
  res$period1 <- period1
  return(res)
}


D_term <- function(M, Rrup, Rhypo, Ftype, V){
  # The distance subroutine function of BCHydro 2011 for subduction earthquakes

  # Notice: input argument V is a list
  if(Ftype == 0){
    R = Rrup + V$c4*exp((M-6.0)*V$a9)
    DTerm = V$a1 + V$a4*V$DeltaC1 + (V$a2 + V$a14*Ftype + V$a3*(M - 7.8))*log(R) + V$a6*Rrup + V$a10*Ftype
  }else if(Ftype == 1){
    R = Rhypo + V$c4*exp((M-6.0)*V$a9)
    DTerm = V$a1 + V$a4*V$DeltaC1 + (V$a2 + V$a14*Ftype + V$a3*(M - 7.8))*log(R) + V$a6*Rhypo + V$a10*Ftype
  }
  return(DTerm)
}


M_term <- function(M, V){

  # The magnitude subroutine function of BCHydro 2011 for subduction earthquakes

  testmag = (7.8 + V$DeltaC1)
  if(M < testmag){
    MTerm = V$a4*(M-testmag) + V$a13*(10.0-M)^2
  }else{
    MTerm = V$a5*(M-testmag) + V$a13*(10.0-M)^2
  }
  return(MTerm)
}

Depth_term <- function(depth, Ftype, V){

  # The depth subroutine function of BCHydro 2011 for subduction earthquakes

  # Notice: input argument V is a list
  DepthTerm = V$a11*(depth - 60.0)*Ftype
  return(DepthTerm)
}


Arc_term <- function(Rrup, Rhypo, faba, Ftype, V){

  # The arc effect subroutine function of BCHydro 2011 for subduction earthquakes

  # Notice: input argument V is a list
  if(Ftype == 1){
    ArcTerm = (V$a7 + V$a8*log((max(Rhypo, 85.0)/40.0)))*faba
  }else if(Ftype== 0){
    ArcTerm = (V$a15 + V$a16*log((max(Rrup, 100.0)/40.0)))*faba
  }
  return(ArcTerm)
}

Site_term <- function(Vs30,PGA_rock,V){

  # The site effect subroutine function of BCHydro 2011 for subduction earthquakes

  # Notice: input argument V is a list
  if(Vs30 > 1000.0){
    VsStar = 1000.0
  }else{
    VsStar = Vs30
  }
  if(Vs30 > V$vLin){
    SiteTerm = V$a12*log(VsStar/V$vLin) + V$b_soil*V$n*log(VsStar/V$vLin)
  }else{
    SiteTerm = V$a12*log(VsStar/V$vLin) - V$b_soil*log(PGA_rock + V$c) +
      V$b_soil*log(PGA_rock + V$c*(VsStar/V$vLin)^V$n)
  }
  return(SiteTerm)
}


calc_val <- function(M, Rrup, Rhypo, PGA_rock, Vs30, Ftype, faba, depth, V){

  # The subroutine function of BCHydro 2011 for PGA for subduction earthquakes

  # Notice: input argument V is a list
  # calculate predicted value
  # Compute the R term and base model based on either Rupture Distance (Interface events)
  # of Hypocentral distance (Intraslab events).
  if(Ftype == 0){
    R = Rrup + V$c4*exp((M-6.0)*V$a9)
  }else if(Ftype == 1){
    R = Rhypo + V$c4*exp((M-6.0)*V$a9)
  }
  X = D_term(M, Rrup, Rhypo, Ftype, V) +  M_term(M, V) + Depth_term(depth, Ftype, V) +
    Arc_term(Rrup, Rhypo, faba, Ftype, V) + Site_term(Vs30, PGA_rock, V)
  return(X)
}

calc_val2 <- function(M, Rrup, Rhypo, PGA_rock, Vs30, Ftype, faba, depth, V){

  # The subroutine function of BCHydro 2011 for PSA for subduction earthquakes

  # Notice: input argument V is a list
  # calculate predicted value
  # Compute the R term and base model based on
  # either Rupture Distance (Interface events)
  # of Hypocentral distance (Intraslab events).
  if(Ftype == 0){
    R = Rrup + V$c4*exp((M-6.0)*V$a9)
  }else if(Ftype == 1){
    R = Rhypo + V$c4*exp((M-6.0)*V$a9)
  }
  X = D_term(M, Rrup, Rhypo, Ftype, V) +  M_term(M, V) + Depth_term(depth, Ftype, V) +
    Arc_term(Rrup, Rhypo, faba, Ftype, V) + Site_term(Vs30, PGA_rock, V)
  return(X)
}

get_BCHydro_constants <- function(index, flag_deltaC1){

  # The constant terms subroutine function of BCHydro 2011 for subduction earthquakes

  ##### return a list
  # get relevant constants
  # arrays with values by index
  period  = c(0.0, 0.02, 0.05, 0.075, 0.1, 0.15, 0.2,    0.25,     0.3,     0.4,     0.5,     0.6,    0.75,
              1.00,     1.5,    2.00,     2.5,    3.00,     4.0,     5.0,     6.0,     7.5,     10)
  vLin = c(865.1,   865.1,  1053.5,  1085.7,  1032.5,   877.6,   748.2,   654.3,   587.1,   503.0,   456.6,
           430.3,   410.5,   400.0,   400.0,   400.0,   400.0,   400.0,   400.0,   400.0,   400.0,   400.0, 400.0)
  b_soil = c(-1.186,  -1.186,  -1.346,  -1.471,  -1.624,  -1.931,  -2.188,  -2.381,  -2.518,  -2.657,  -2.669,
             -2.599,  -2.401,  -1.955,  -1.025,  -0.299,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 0.0)
  a1 = c(4.2203,  4.2203,  4.5371,  5.0733,  5.2892,  5.4563,  5.2684,  5.0594,  4.7945,  4.4644,  4.0181,  3.6055,
         3.2174,  2.7981,  2.0123,  1.4128,  0.9976,  0.6443,  0.0657, -0.4624, -0.9809, -1.6017, -2.2937)
  a2 = c(-1.35,     -1.35,    -1.4,   -1.45,   -1.45,   -1.45,    -1.4,   -1.35,   -1.28,   -1.18,   -1.08,
         -0.99,   -0.91,   -0.85,   -0.77,   -0.71,   -0.67,   -0.64,   -0.58,   -0.54,    -0.5,   -0.46,   -0.4)
  a6 = c(-0.0012, -0.0012, -0.0012, -0.0012, -0.0012, -0.0014, -0.0018, -0.0023, -0.0027, -0.0035, -0.0044, -0.0050,
         -0.0058, -0.0062, -0.0064, -0.0064, -0.0064, -0.0064, -0.0064, -0.0064, -0.0064, -0.0064, -0.0064)
  a7 = c(1.0988,  1.0988,  1.2536,  1.4175,  1.3997,  1.3582,  1.1648,  0.994,   0.8821,  0.7046,  0.5799,  0.5021,
         0.3687,  0.1746,  -0.082, -0.2821, -0.4108, -0.4466, -0.4344, -0.4368, -0.4586, -0.4433, -0.4828)
  a8 = c(-1.42,   -1.42,   -1.65,   -1.80,   -1.80,   -1.69,   -1.49,  -1.30,    -1.18,   -0.98,   -0.82,   -0.70,
         -0.54,   -0.34,   -0.05,    0.12,    0.25,    0.30,    0.30,    0.30,    0.30,    0.30,    0.30)
  a10 = c(3.12,    3.12,    3.37,    3.37,    3.33,    3.25,    3.03,    2.8,     2.59,     2.2,    1.92,     1.7,
          1.42,     1.1,     0.7,     0.7,     0.7,     0.7,     0.7,     0.7,     0.7,     0.7,     0.7)
  a11 = c(0.0130,  0.0130,  0.0130,  0.0130,  0.0130,  0.0130,  0.0129, 0.0129,   0.0128,  0.0127,  0.0125,  0.0124,
          0.0120,  0.0114,  0.0100,  0.0085,  0.0069,  0.0054,  0.0027,  0.0005, -0.0013, -0.0033, -0.0060)
  a12 = c(0.980,   0.980,   1.288,   1.483,   1.613,   1.882,   2.076,  2.248,    2.348,   2.427,   2.399,   2.273,
          1.993,   1.470,   0.408,  -0.401,  -0.723,  -0.673,  -0.627,  -0.596,  -0.566,  -0.528,  -0.504)
  a13 = c(-0.0135, -0.0135, -0.0138, -0.0142, -0.0145, -0.0153, -0.0162,-0.0172,  -0.0183, -0.0206, -0.0231,
          -0.0256, -0.0296, -0.0363, -0.0493, -0.0610, -0.0711, -0.0798, -0.0935, -0.0980, -0.0980, -0.0980, -0.0980)
  a14 = c(-0.40,   -0.40,   -0.40,   -0.40,   -0.40,   -0.40,   -0.35,  -0.31,    -0.28,   -0.23,   -0.19,   -0.16,
          -0.12,   -0.07,    0.00,    0.00,    0.00,    0.00,    0.00,    0.00,    0.00,    0.00,    0.00)
  a15 = c(0.9969,  0.9969,  1.1030,  1.2732,  1.3042,  1.2600,  1.2230, 1.1600,   1.0500,  0.8000,  0.6620,
          0.5800,  0.4800,  0.3300,  0.3100,  0.3000,  0.3000,  0.3000,  0.3000,  0.3000,  0.3000,  0.3000,  0.3000)
  a16 = c(-1.00,   -1.00,   -1.18,   -1.36,   -1.36,   -1.30,   -1.25,  -1.17,    -1.06,   -0.78,   -0.62,
          -0.50,   -0.34,   -0.14,    0.00,    0.00,    0.00,    0.00,    0.00,    0.00,    0.00,    0.00,    0.00)
  # intra-event variability
  sigs = c(0.603,   0.603,   0.603,   0.603,   0.603,   0.603,   0.603,  0.603,    0.603,   0.603,   0.603,   0.603,
           0.603, 0.603,   0.603,   0.603,   0.603,   0.603,   0.603,   0.603,   0.603,   0.603,   0.603)
  # inter-event variability
  sigt = c(0.482,   0.482,   0.482,   0.482,   0.482,   0.482,   0.482,  0.482,    0.482,   0.482,   0.482, 0.482,
           0.482, 0.482,   0.482,   0.482,   0.482,   0.482,   0.482,   0.482,   0.482,   0.482,   0.482)
  deltaC1_d = c(0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,    0.2,      0.2,  0.1437,     0.1,
                0.0737,  0.0415,     0.0, -0.0585,    -0.1, -0.1550,    -0.2,    -0.2,    -0.2,    -0.2,    -0.2,
                -0.2)
  deltaC1_i = c(0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,    0.0,      0.0,     0.0,
                0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
                0.0,     0.0,     0.0)
  # Constant parameters
  n = 1.18
  c = 1.88
  a3 = 0.1
  a4 = 0.9
  a5 = 0.0
  a9 = 0.4
  c4 = 10.0
  c1 = 7.8

  periodc = period[index]
  vLin = vLin[index]
  b_soil = b_soil[index]
  a1 = a1[index]
  a2 = a2[index]
  a3 = a3
  a4 = a4
  a5 = a5
  a6 = a6[index]
  a7 = a7[index]
  a8 = a8[index]
  a9 = a9
  a10 = a10[index]
  a11 = a11[index]
  a12 = a12[index]
  a13 = a13[index]
  a14 = a14[index]
  a15 = a15[index]
  a16 = a16[index]

  n = n
  c = c
  c1 = c1
  c4 = c4
  sigs = sigs[index]
  sigt = sigt[index]

  if(flag_deltaC1 == 0){  # DeltaC1 independent from period
    DeltaC1 = deltaC1_i[index]
  }else if(flag_deltaC1 == 1){  # DeltaC1 dependent from period (Tohoku and Maule example)
    DeltaC1 = deltaC1_d[index]
  }
  constants <- list()
  constants$period <- periodc
  constants$vLin <- vLin
  constants$b_soil <- b_soil
  constants$a1 <- a1
  constants$a2 <- a2
  constants$a3 <- a3
  constants$a4 <- a4
  constants$a5 <- a5
  constants$a6 <- a6
  constants$a7 <- a7
  constants$a8 <- a8
  constants$a9 <- a9
  constants$a10 <- a10
  constants$a11 <- a11
  constants$a12 <- a12
  constants$a13 <- a13
  constants$a14 <- a14
  constants$a15 <- a15
  constants$a16 <- a16
  constants$n <- n
  constants$c <- c
  constants$c1 <- c1
  constants$c4 <- c4
  constants$sigs <- sigs
  constants$sigt <- sigt
  constants$DeltaC1 <- DeltaC1
  return(constants)  # return a list
}

get_BCHydro_constants_760 <- function(index_760){

  # The constant terms subroutine function of BCHydro 2011 for Vs760 for subduction earthquakes

  # return a list
  # get relevant constants for computing the amplification factors
  # arrays with values by index
  # Paramenters used in the conversion from Hard Rock to Vs760

  period_2Vs760 = c(0.00,    0.01,   0.02,   0.03,    0.05,   0.075,     0.1,    0.15,     0.2,     0.3,   0.5,
                    1.00,   2.00,    3.3,    5.0,  10.00)
  b1_2Vs760 = c(0.044,   0.044, -0.219, -0.434,  -0.424,  -0.267,  -0.102,   0.237,    0.29,   0.333, 0.379,
                0.481,  0.418,   0.37,  0.307,  0.241)
  b2_2Vs760 = c(-0.16,   -0.16, -0.354, -0.273,  -0.126,  -0.031,   0.015,  -0.006,     0.0,  -0.006, 0.019,
                0.018, -0.002,    0.0,    0.0,    0.0)
  b3_2Vs760 = c(0.0072, 0.0072,  0.031, 0.0179, -0.0042, -0.0238, -0.0294, -0.0333, -0.0384, -0.0289,  0.00,
                0.0174,  0.011, 0.0054, 0.0028, 0.0042)
  c1_2Vs760 = c(-3.0,     -3.0,   -2.5,   -2.5,    -2.5,    -2.0,    -2.0,    -1.5,    -1.0,    -0.5,  -2.0,
                -2.0,   -2.5,   -3.0,   -4.0,   -5.0)
  c2_2Vs760 = c(1.0,      1.0,    1.8,    2.0,     1.9,     1.7,     1.6,     1.4,     1.3,     1.2,   0.9,
                0.4,   -0.2,   -0.8,   -1.4,   -2.9)

  period_2Vs760 = period_2Vs760[index_760]
  b1_2Vs760 = b1_2Vs760[index_760]
  b2_2Vs760 = b2_2Vs760[index_760]
  b3_2Vs760 = b3_2Vs760[index_760]
  c1_2Vs760 = c1_2Vs760[index_760]
  c2_2Vs760 = c2_2Vs760[index_760]

  constants760 <- list()
  constants760$period_2Vs760 <- period_2Vs760
  constants760$b1_2Vs760 <- b1_2Vs760
  constants760$b2_2Vs760 <- b2_2Vs760
  constants760$b3_2Vs760 <- b3_2Vs760
  constants760$c1_2Vs760 <- c1_2Vs760
  constants760$c2_2Vs760 <- c2_2Vs760

  return(constants760)
}

BCHHR2Vs760 <- function(T_760,LSa_noamp){

  # The amplifcation subroutine function of BCHydro 2011 for subduction earthquakes

  # Now compute the amp factor based on the given hard rock spectral acceleration.
  period_2Vs760 = c(0.00,0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.3,0.5,1.00,2.00,3.3,5.0,10.00)
  # Find the requested spectral period and corresponding coefficients
  nT760 = length(T_760) # the size of T_760 shoud always be one!
  # specT starts from 0.02!!!
  nperi760 = length(period_2Vs760)
  period1_760 = T_760
  # interpolate the coefficients between periods if neccesary
  if(length(which(period_2Vs760 == period1_760)) == 0){
    index_low =  max(which(period_2Vs760 < period1_760))
    index_hi =  min(which(period_2Vs760 > period1_760))

    V760_low = get_BCHydro_constants_760(index_low)   # return a list
    V760_hi = get_BCHydro_constants_760(index_hi)    # return a list
    X_period760 = c(period_2Vs760[index_low], period_2Vs760[index_hi])
    B1_constants = c(V760_low$b1_2Vs760, V760_hi$b1_2Vs760)
    B2_constants = c(V760_low$b2_2Vs760, V760_hi$b2_2Vs760)
    B3_constants = c(V760_low$b3_2Vs760, V760_hi$b3_2Vs760)
    C1_constants = c(V760_low$c1_2Vs760, V760_hi$c1_2Vs760)
    C2_constants = c(V760_low$c2_2Vs760, V760_hi$c2_2Vs760)

    b1_2Vs760 = approx(x = X_period760, y = B1_constants, xout = period1_760, rule = 2)$y   ## V760 is a list
    b2_2Vs760 = approx(x = X_period760, y = B2_constants, xout = period1_760, rule = 2)$y
    b3_2Vs760 = approx(x = X_period760, y = B3_constants, xout = period1_760, rule = 2)$y
    c1_2Vs760 = approx(x = X_period760, y = C1_constants, xout = period1_760, rule = 2)$y
    c2_2Vs760 = approx(x = X_period760, y = C2_constants, xout = period1_760, rule = 2)$y

    V760 <- list()  ## make V760 as a list
    V760$b1_2Vs760 <- b1_2Vs760
    V760$b2_2Vs760 <- b2_2Vs760
    V760$b3_2Vs760 <- b3_2Vs760
    V760$c1_2Vs760 <- c1_2Vs760
    V760$c2_2Vs760 <- c2_2Vs760
  }else{
    index_760 = which(abs((period_2Vs760 - period1_760)) < 0.0001)  # Identify the period
    # get constants for the given index value
    V760 = get_BCHydro_constants_760(index_760)  # return a list
  }
  # Check for LsaRock value between bounding values C1 and C2
  # Convert Hard Rock SA from Gals to G

  LsaRock = LSa_noamp  # - 6.89 %from fortran
  # LsaRock = exp(Sa_noamp)
  if(LsaRock < V760$c1_2Vs760){
    AmpFac = V760$b1_2Vs760
  }else if(LsaRock > V760$c2_2Vs760){
    AmpFac = V760$b1_2Vs760 + V760$b2_2Vs760*(V760$c2_2Vs760 - V760$c1_2Vs760) +
      V760$b3_2Vs760*(V760$c2_2Vs760 - V760$c1_2Vs760)*(V760$c2_2Vs760 - V760$c1_2Vs760)
  }else{
    AmpFac = V760$b1_2Vs760 + V760$b2_2Vs760*(LsaRock - V760$c1_2Vs760) +
      V760$b3_2Vs760*(LsaRock - V760$c1_2Vs760)*(LsaRock - V760$c1_2Vs760)
  }
  return(AmpFac)
}

BCHydro_sigma <- function(V){

  # The sigma subroutine function of BCHydro 2011 for subduction earthquakes

  # return a lit
  # calculate the sigma
  # use the published coefficients for the geometric mean
  sigma = V$sigs
  tau = V$sigt
  tsigma = sqrt(sigma^2 + tau^2)

  res <- list()
  res$tsigma <- tsigma
  res$sigma <- sigma
  res$tau <- tau
  return(res)
}




