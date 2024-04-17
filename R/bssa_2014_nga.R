#' The GMPE of BSSA 2014 and updated versions with corrected CA subregional anelastic attenuation and basin effects
#'
#' This function calculates the ground motion median values and standard deviations
#' @param M Moment magnitude, a numeric value
#' @param T Period (sec); Use Period = -1 for PGV computation.
#' Use 1000 (by default) for output the array of Sa with original NGA West2 periods
#' @param Rjb Joyner-Boore distance (km); closest distance (km) to surface
#' projection of rupture plane
#' @param Fault_Type The indicator of fault type: 0 for unspecified fault;
#' 1 for strike-slip fault; 2 for normal fault; 3 for reverse fault
#' @param region Region indicator: 0 for global (incl. Taiwan); 1 for California;
#' 2 for Japan; 3 for China or Turkey; 4 for Italy
#' @param z1 Basin depth (km): depth from the groundsurface to the 1km/s shear-wave horizon.
#' 999 if unknown.
#' @param Vs30 Shear wave velocity averaged over top 30 m (in m/s)
#' @param return.type The indicator specifies which type of data return:
#' 1 for the med (in g)/sigma/phi/tau/period; 2 for F_E/F_P/F_P_GS(geomteric spreading)/F_P_AA(anelastic attenuation)/F_S;
#' 3 for r/Ln_Flin/Ln_Fnlin/f2/F_dz1/PGAr
#' @param CA_subreg_path The flag indicating if a CA subregional anelastic attenuation model (Buckreis et al., 2023) is applied.
#' Buckreis, Stewart, Brandenberg, and Wang (2023) "Subregional Anelastic Attenuation Model for California".
#' @param site_lon The longitude of site (in WGS84 degree). Required if CA_subreg_path is TRUE.
#' @param site_lat The latitude of site (in WGS84 degree). Required if CA_subreg_path is TRUE.
#' @param clst_lon Longitude (in WGS84 degrees) of closest point on the surface projection of the rupture surface to the site.
#' Required if CA_subreg_path is TRUE.
#' @param clst_lat Latitude (in WGS84 degrees) of closest point on the surface projection of the rupture surface to the site.
#' Required if CA_subreg_path is TRUE.
#' @return See return.type for the returned results
#' @examples bssa_2014_nga(M = 5, T = 1000, Rjb = 85, Fault_Type = 1, region = 1, z1 = 999, Vs30 = 350)
#' bssa_2014_nga(M = 5, T = 1000, Rjb = 85, Fault_Type = 1, region = 1, z1 = 999, Vs30 = 350, CA_subreg_path = TRUE, site_lon = -118.8713, site_lat = 35.80407, clst_lon = -116.9274, clst_lat = 33.6204)
#' @references Boore, D. M., Stewart, J. P., Seyhan, E., and Atkinson, G. M. (2014).
#' NGA-West2 Equations for Predicting PGA, PGV, and 5% Damped PSA for
#' Shallow Crustal Earthquakes. Earthquake Spectra, 30(3), 1057-1085.
#' @export
#' @importFrom stats approx
bssa_2014_nga <- function(M, T = 1000, Rjb, Fault_Type, region, z1 = 999, Vs30, return.type = 1,
                          CA_subreg_path = FALSE, site_lon = NULL, site_lat = NULL,
                          clst_lon = NULL, clst_lat = NULL){

  coeffs <- bssa_2014_coeffs

  period = c(-1, 0, 0.010, 0.020, 0.022, 0.025, 0.029, 0.030, 0.032, 0.035, 0.036, 0.040, 0.042,
             0.044, 0.045, 0.046, 0.048, 0.050, 0.055, 0.060, 0.065, 0.067, 0.070, 0.075,
             0.080, 0.085, 0.090, 0.095, 0.100, 0.110, 0.120, 0.130, 0.133, 0.140, 0.150,
             0.160, 0.170, 0.180, 0.190, 0.200, 0.220, 0.240, 0.250, 0.260, 0.280, 0.290,
             0.300, 0.320, 0.340, 0.350, 0.360, 0.380, 0.400, 0.420, 0.440, 0.450, 0.460,
             0.480, 0.500, 0.550, 0.600, 0.650, 0.667, 0.700, 0.750, 0.800, 0.850, 0.900,
             0.950, 1.000, 1.100, 1.200, 1.300, 1.400, 1.500, 1.600, 1.700, 1.800, 1.900,
             2.000, 2.200, 2.400, 2.500, 2.600, 2.800, 3.000, 3.200, 3.400, 3.500, 3.600,
             3.800, 4.000, 4.200, 4.400, 4.600, 4.800, 5.000, 5.500, 6.000, 6.500, 7.000,
             7.500, 8.000, 8.500, 9.000, 9.500, 10.000)

  U = (Fault_Type == 0)
  SS = (Fault_Type == 1)
  NS = (Fault_Type == 2)
  RS = (Fault_Type == 3)

  # calculate segment distances for the updated anelastic attenuation model
  if (CA_subreg_path) {
    if (is.null(site_lon) | is.null(site_lat) | is.null(clst_lat) | is.null(clst_lon)) {
      print('Invaid input lon/lat. Use default BSSA14 path!')
      res_dist <- NA
    } else {
      res_dist <- CA_subreg_path_calc(site_lon = site_lon, site_lat = site_lat, clst_lon = clst_lon,
                                      clst_lat = clst_lat)
      if (!((abs(res_dist$site_source_dist - sum(res_dist$sub_dist)) < 500) & !is.na(res_dist$source_region_idx))) {
        print('Significant distance segment outside of CA subregions or the source is outside of CA. Use default BSSA14 path!')
        res_dist <- NA
      }
    }
  }

  if(length(T) == 1 & T[1] == 1000){ # Compute median and sigma with pre-defined period
    med <- rep(0,length(period))
    sig <- rep(0,length(period))
    phi <- rep(0,length(period))
    tau <- rep(0,length(period))
    F_E <- rep(0,length(period))
    F_P <- rep(0,length(period))
    F_P_GS <- rep(0,length(period))
    F_P_AA <- rep(0,length(period))
    F_S <- rep(0,length(period))
    for(ip in 1:length(period)){

      # calculate updated anelastic attenuation
      r <- sqrt(Rjb^2 + coeffs$h..km.[ip]^2)
      if (CA_subreg_path) {
        if (length(res_dist) > 1) {
          w_r <- res_dist$sub_dist / sum(res_dist$sub_dist)
          tmp_dc3_coef <- CA_subreg_coeffs[ip, 3:12]
          tmp_deltac3 <- sum(w_r * tmp_dc3_coef)
          tmp_deltac0 <- CA_subreg_coeffs[ip, res_dist$source_region_idx + 12]
          F_AA <- (coeffs$c3[ip] + tmp_deltac3) * (r - 1) + tmp_deltac0
        } else {
          F_AA <- NA
        }
      } else {
        F_AA <- NA
      }

      res1 <- bssa_2014_subroutine(M, ip, Rjb, U, SS, NS, RS, region, z1, Vs30, coeffs = coeffs, CA_F_AA = F_AA)
      med[ip] <- res1$med
      sig[ip] <- res1$sigma
      phi[ip] <- res1$phi
      tau[ip] <- res1$tau
      res2 <- bssa_2014_subroutine(M, ip, Rjb, U, SS, NS, RS, region, z1, Vs30, return.type = 2, coeffs = coeffs, CA_F_AA = F_AA)
      F_E[ip] <- res2$F_E
      F_P[ip] <- res2$F_P
      F_P_GS[ip] <- res2$F_P_GS
      F_P_AA[ip] <- res2$F_P_AA
      F_S[ip] <- res2$F_S
    }
    T <- period
  }else{  # Compute median and sigma with user-defined period
    med  <- rep(0, length(T))
    sig <- rep(0, length(T))
    phi <- rep(0, length(T))
    tau <- rep(0, length(T))
    F_E <- rep(0, length(T))
    F_P <- rep(0, length(T))
    F_P_GS <- rep(0, length(T))
    F_P_AA <- rep(0, length(T))
    F_S <- rep(0, length(T))
    r <- rep(0, length(T))
    ln_Flin <- rep(0, length(T))
    ln_Fnlin <- rep(0, length(T))
    f2 <- rep(0, length(T))
    F_dz1 <- rep(0, length(T))
    PGAr <- rep(0, length(T))
    period1 <- T
    for(i in 1:length(T)){
      Ti = T[i]
      if(length(which(abs(period-Ti) < 0.0001)) == 0){ # The user defined period requires interpolation
        idx <- which(period < Ti)
        T_low = max(period[idx])
        idx <- which(period > Ti)
        T_high = min(period[idx])
        ip_low  = which(period == T_low)
        ip_high = which(period==T_high)

        # calculate updated anelastic attenuation
        r <- sqrt(Rjb^2 + coeffs$h..km.[ip_low]^2)
        if (CA_subreg_path) {
          if (length(res_dist) > 1) {
            w_r <- res_dist$sub_dist / sum(res_dist$sub_dist)
            tmp_dc3_coef <- CA_subreg_coeffs[ip_low, 3:12]
            tmp_deltac3 <- sum(w_r * tmp_dc3_coef)
            tmp_deltac0 <- CA_subreg_coeffs[ip_low, res_dist$source_region_idx + 12]
            F_AA <- (coeffs$c3[ip_low] + tmp_deltac3) * (r - 1) + tmp_deltac0
          } else {
            F_AA <- NA
          }
        } else {
          F_AA <- NA
        }

        res_low <- bssa_2014_subroutine(M, ip_low, Rjb, U, SS, NS, RS, region, z1, Vs30, coeffs = coeffs, CA_F_AA = F_AA)
        Sa_low <- res_low$med
        sigma_low <- res_low$sigma
        phi_low <- res_low$phi
        tau_low <- res_low$tau
        res_low2 <- bssa_2014_subroutine(M, ip_low, Rjb, U, SS, NS, RS, region, z1, Vs30, return.type = 2, coeffs = coeffs, CA_F_AA = F_AA)
        FE_low <- res_low2$F_E
        FP_low <- res_low2$F_P
        FP_GS_low <- res_low2$F_P_GS
        FP_AA_low <- res_low2$F_P_AA
        FS_low <- res_low2$F_S
        res_low3 <- bssa_2014_subroutine(M, ip_low, Rjb, U, SS, NS, RS, region, z1, Vs30, return.type = 3, coeffs = coeffs, CA_F_AA = F_AA)
        r_low <- res_low3$r
        ln_Flin_low <- res_low3$ln_Flin
        ln_Fnlin_low <- res_low3$ln_Fnlin
        f2_low <- res_low3$f2
        F_dz1_low <- res_low3$F_dz1
        PGAr_low <- res_low3$PGAr


        # calculate updated anelastic attenuation
        r <- sqrt(Rjb^2 + coeffs$h..km.[ip_high]^2)
        if (CA_subreg_path) {
          if (length(res_dist) > 1) {
            w_r <- res_dist$sub_dist / sum(res_dist$sub_dist)
            tmp_dc3_coef <- CA_subreg_coeffs[ip_high, 3:12]
            tmp_deltac3 <- sum(w_r * tmp_dc3_coef)
            tmp_deltac0 <- CA_subreg_coeffs[ip_high, res_dist$source_region_idx + 12]
            F_AA <- (coeffs$c3[ip_high] + tmp_deltac3) * (r - 1) + tmp_deltac0
          } else {
            F_AA <- NA
          }
        } else {
          F_AA <- NA
        }

        res_high <- bssa_2014_subroutine(M, ip_high, Rjb, U, SS, NS, RS, region, z1, Vs30, coeffs = coeffs, CA_F_AA = F_AA)
        Sa_high <- res_high$med
        sigma_high <- res_high$sigma
        phi_high <- res_high$phi
        tau_high <- res_high$tau
        res_high2 <- bssa_2014_subroutine(M, ip_high, Rjb, U, SS, NS, RS, region, z1, Vs30, return.type = 2, coeffs = coeffs, CA_F_AA = F_AA)
        FE_high <- res_high2$F_E
        FP_high <- res_high2$F_P
        FP_GS_high <- res_high2$F_P_GS
        FP_AA_high <- res_high2$F_P_AA
        FS_high <- res_high2$F_S
        res_high3 <- bssa_2014_subroutine(M, ip_high, Rjb, U, SS, NS, RS, region, z1, Vs30, return.type = 3, coeffs, CA_F_AA = F_AA)
        r_high <- res_high3$r
        ln_Flin_high <- res_high3$ln_Flin
        ln_Fnlin_high <- res_high3$ln_Fnlin
        f2_high <- res_high3$f2
        F_dz1_high <- res_high3$F_dz1
        PGAr_high <- res_high3$PGAr

        x = c(T_low,T_high)
        Y_sa = c(Sa_low,Sa_high)
        Y_sigma = c(sigma_low,sigma_high)
        Y_phi <- c(phi_low, phi_high)
        Y_tau <- c(tau_low, tau_high)
        FE_temp <- c(FE_low, FE_high)
        FP_temp <- c(FP_low, FP_high)
        FP_GS_temp <- c(FP_GS_low, FP_GS_high)
        FP_AA_temp <- c(FP_AA_low, FP_AA_high)
        FS_temp <- c(FS_low, FS_high)
        r_temp <- c(r_low, r_high)
        ln_Flin_temp <- c(ln_Flin_low, ln_Flin_high)
        ln_Fnlin_temp <- c(ln_Fnlin_low, ln_Fnlin_high)
        f2_temp <- c(f2_low, f2_high)
        F_dz1_temp <- c(F_dz1_low, F_dz1_high)
        PGAr_temp <- c(log(PGAr_low), log(PGAr_high))

        temp_med <- approx(x, Y_sa, Ti, rule = 2)$y
        med[i] = temp_med    # linear interpolation
        temp_sig <- approx(x, Y_sigma, Ti, rule = 2)$y
        sig[i] <- temp_sig
        temp_phi <- approx(x, Y_phi, Ti, rule = 2)$y
        phi[i] <- temp_phi
        temp_tau <- approx(x, Y_tau, Ti, rule = 2)$y
        tau[i] <- temp_tau

        FE_med <- approx(x, FE_temp, Ti, rule = 2)$y
        FP_med <- approx(x, FP_temp, Ti, rule = 2)$y
        FP_GS_med <- approx(x, FP_GS_temp, Ti, rule = 2)$y
        FP_AA_med <- approx(x, FP_AA_temp, Ti, rule = 2)$y
        FS_med <- approx(x, FS_temp, Ti, rule = 2)$y
        F_E[i] <- FE_med
        F_P[i] <- FP_med
        F_P_GS[i] <- FP_GS_med
        F_P_AA[i] <- FP_AA_med
        F_S[i] <- FS_med

        r_med <- approx(x, r_temp, Ti, rule = 2)$y
        ln_Flin_med <- approx(x, ln_Flin_temp, Ti, rule = 2)$y
        ln_Fnlin_med <- approx(x, ln_Fnlin_temp, Ti, rule = 2)$y
        f2_med <- approx(x, f2_temp, Ti, rule = 2)$y
        F_dz1_med <- approx(x, F_dz1_temp, Ti, rule = 2)$y
        PGAr_med <- approx(x, PGAr_temp, Ti, rule = 2)$y
        r[i] <- r_med
        ln_Flin[i] <- ln_Flin_med
        ln_Fnlin[i] <- ln_Fnlin_med
        f2[i] <- f2_med
        F_dz1[i] <- F_dz1_med
        PGAr[i] <- exp(PGAr_med)

      }else{
        ip_T = which(abs((period - Ti)) < 0.0001)

        # calculate updated anelastic attenuation
        r <- sqrt(Rjb^2 + coeffs$h..km.[ip_T]^2)
        if (CA_subreg_path) {
          if (length(res_dist) > 1) {
            w_r <- res_dist$sub_dist / sum(res_dist$sub_dist)
            tmp_dc3_coef <- CA_subreg_coeffs[ip_T, 3:12]
            tmp_deltac3 <- sum(w_r * tmp_dc3_coef)
            tmp_deltac0 <- CA_subreg_coeffs[ip_T, res_dist$source_region_idx + 12]
            F_AA <- (coeffs$c3[ip_T] + tmp_deltac3) * (r - 1) + tmp_deltac0
          } else {
            F_AA <- NA
          }
        } else {
          F_AA <- NA
        }

        res1 <- bssa_2014_subroutine(M, ip_T, Rjb, U, SS, NS, RS, region, z1, Vs30, coeffs = coeffs, CA_F_AA = F_AA)
        med[i] <- res1$med
        sig[i] <- res1$sigma
        phi[i] <- res1$phi
        tau[i] <- res1$tau
        res2 <- bssa_2014_subroutine(M, ip_T, Rjb, U, SS, NS, RS, region, z1, Vs30, return.type = 2, coeffs = coeffs, CA_F_AA = F_AA)
        F_E[i] <- res2$F_E
        F_P[i] <- res2$F_P
        F_P_GS[i] <- res2$F_P_GS
        F_P_AA[i] <- res2$F_P_AA
        F_S[i] <- res2$F_S
        res3 <- bssa_2014_subroutine(M, ip_T, Rjb, U, SS, NS, RS, region, z1, Vs30, return.type = 3, coeffs = coeffs, CA_F_AA = F_AA)
        r[i] <- res3$r
        ln_Flin[i] <- res3$ln_Flin
        ln_Fnlin[i] <- res3$ln_Fnlin
        f2[i] <- res3$f2
        F_dz1[i] <- res3$F_dz1
        PGAr[i] <- res3$PGAr
      }
    }
  }
  if(return.type == 1){
    res <- list()
    res$med <- med
    res$sigma <- sig
    res$phi <- phi
    res$tau <- tau
    res$period <- T
    return(res)
  }else if(return.type == 2){
    effs <- list()
    effs$F_E <- F_E
    effs$F_P <- F_P
    effs$F_P_GS <- F_P_GS
    effs$F_P_AA <- F_P_AA
    effs$F_S <- F_S
    effs$period <- T
    return(effs)
  }else if(return.type == 3){
    res <- list()
    res$r <- r
    res$ln_Flin <- ln_Flin
    res$ln_Fnlin <- ln_Fnlin
    res$f2 <- f2
    res$F_dz1 <- F_dz1
    res$PGAr <- PGAr
    res$period <- T
    return(res)
  }
}

#' The subroutine of GMPE of BSSA 2014
#'
#' This is a subroutine of BSSA 2014
#' @param M Moment magnitude, a numeric value
#' @param ip The index of working period in the per-defined periods: (-1,
#' 0, 0.010, 0.020, 0.022, 0.025, 0.029, 0.030, 0.032, 0.035, 0.036, 0.040, 0.042,
#' 0.044, 0.045, 0.046, 0.048, 0.050, 0.055, 0.060, 0.065, 0.067, 0.070, 0.075,
#' 0.080, 0.085, 0.090, 0.095, 0.100, 0.110, 0.120, 0.130, 0.133, 0.140, 0.150,
#' 0.160, 0.170, 0.180, 0.190, 0.200, 0.220, 0.240, 0.250, 0.260, 0.280, 0.290,
#' 0.300, 0.320, 0.340, 0.350, 0.360, 0.380, 0.400, 0.420, 0.440, 0.450, 0.460,
#' 0.480, 0.500, 0.550, 0.600, 0.650, 0.667, 0.700, 0.750, 0.800, 0.850, 0.900,
#' 0.950, 1.000, 1.100, 1.200, 1.300, 1.400, 1.500, 1.600, 1.700, 1.800, 1.900,
#' 2.000, 2.200, 2.400, 2.500, 2.600, 2.800, 3.000, 3.200, 3.400, 3.500, 3.600,
#' 3.800, 4.000, 4.200, 4.400, 4.600, 4.800, 5.000, 5.500, 6.000, 6.500, 7.000,
#' 7.500, 8.000, 8.500, 9.000, 9.500, 10.000)
#' @param Rjb Joyner-Boore distance (km); closest distance (km) to surface
#' projection of rupture plane
#' @param U The indicator of fault type: 1 for unspecified fault, 0 for otherwise
#' @param SS The indicator of fault type: 1 for strike-slip fault, 0 for otherwise
#' @param NS The indicator of fault type: 1 for normal fault, 0 for otherwise
#' @param RS The indicator of fault type: 1 for reverse fault, 0 for otherwise
#' @param region Region indicator: 0 for global (incl. Taiwan); 1 for California;
#' 2 for Japan; 3 for China or Turkey; 4 for Italy
#' @param z1 Basin depth (km): depth from the groundsurface to the 1km/s shear-wave horizon.
#' 999 if unknown.
#' @param Vs30 Shear wave velocity averaged over top 30 m (in m/s)
#' @param return.type The indicator specifies which type of data return:
#' 1 for the med (in g)/sigma/phi/tau; 2 for F_E/F_P/F_S;
#' 3 for r/Ln_Flin/Ln_Fnlin/f2/F_dz1/PGAr
#' @param coeffs The coefficient table of BSSA 2014. You can use the internal saved data object, bssa_2014_coeffs
#' @param CA_F_AA The updated anelastic attenuation value by the Buckreis et al 2023 CA subregional anelastic attenuation model.
#' If unknown, input NA and the default BSSA14 is used.
#' @return See return.type for the returned results
#' @references Boore, D. M., Stewart, J. P., Seyhan, E., and Atkinson, G. M. (2014).
#' NGA-West2 Equations for Predicting PGA, PGV, and 5% Damped PSA for
#' Shallow Crustal Earthquakes. Earthquake Spectra, 30(3), 1057-1085.
#' @export
bssa_2014_subroutine <- function(M, ip, Rjb, U, SS, NS, RS, region, z1, Vs30, return.type = 1,
                                 coeffs = bssa_2014_coeffs, CA_F_AA = NA){

  mref <- coeffs$Mref[1]
  rref <- coeffs$Rref..km.[1]
  v_ref <- coeffs$Vref...m.s.[1]
  f1 <- coeffs$f1[1]
  f3 <- coeffs$f3[1]
  v1 <- coeffs$V1..m.s.[1]
  v2 <- coeffs$V2..m.s.[1]
  period <- coeffs$Period.sec.
  mh <- coeffs$Mh
  e0 <- coeffs$e0
  e1 <- coeffs$e1
  e2 <- coeffs$e2
  e3 <- coeffs$e3
  e4 <- coeffs$e4
  e5 <- coeffs$e5
  e6 <- coeffs$e6   #### coefficient different at period 10s
  c1 <- coeffs$c1
  c2 <- coeffs$c2
  c3 <- coeffs$c3
  h <- coeffs$h..km.

  deltac3_gloCATW <- coeffs$Dc3..globalCATWNZ.
  deltac3_CHTU <- coeffs$Dc3..ChinaTurkey.
  deltac3_ITJA <- coeffs$Dc3..ItalyJapan.

  coff_c <- coeffs$c
  vc <- coeffs$Vc..m.s.
  f4 <- coeffs$f4
  f5 <- coeffs$f5
  f6 <- coeffs$f6..1.km.
  f7 <- coeffs$f7
  tau1 <- coeffs$t1
  tau2 <- coeffs$t2
  phi1 <- coeffs$phi1
  phi2 <- coeffs$phi2
  dphiR <- coeffs$DfR
  dphiV <- coeffs$DfV
  R1 <- coeffs$R1..km.
  R2 <- coeffs$R2..km.

  ## The source(event function):
  if(M <= mh[ip]){
    F_E <- e0[ip] * U + e1[ip] * SS + e2[ip] * NS + e3[ip] * RS + e4[ip] * (M - mh[ip]) + e5[ip] * (M - mh[ip])^2
  }else{
    F_E = e0[ip] * U + e1[ip] * SS + e2[ip] * NS + e3[ip] * RS + e6[ip] * (M - mh[ip])
  }

  ## The path function:
  if(region == 0 || region == 1){
    deltac3 <- deltac3_gloCATW
  }else if(region == 3){
    deltac3 = deltac3_CHTU
  }else if(region == 2 || region == 4){
    deltac3 = deltac3_ITJA
  }

  r <- sqrt(Rjb^2+h[ip]^2)

  if (!is.na(CA_F_AA)) {
    F_P <- (c1[ip] + c2[ip] * (M - mref)) * log (r / rref) + CA_F_AA
    F_P_GS <- (c1[ip] + c2[ip] * (M - mref)) * log (r / rref)
    F_P_AA <- CA_F_AA
  } else {
    F_P <- (c1[ip] + c2[ip] * (M - mref)) * log (r / rref) + (c3[ip] + deltac3[ip]) * (r - rref)
    F_P_GS <- (c1[ip] + c2[ip] * (M - mref)) * log (r / rref)
    F_P_AA <- (c3[ip] + deltac3[ip]) * (r - rref)
  }

  ## FIND PGAr
  if(Vs30!=v_ref || ip!=2){
    pgar <- bssa_2014_subroutine(M = M, ip = 2, Rjb = Rjb, U = U, SS = SS, NS = NS, RS = RS, region = region,
                                 z1 = z1, Vs30 = v_ref, return.type = 1,
                                 coeffs = bssa_2014_coeffs, CA_F_AA = CA_F_AA)
    PGA_r <- pgar$med
    sigma_r <- pgar$sigma
    phi_r <- pgar$phi
    tau_r <- pgar$tau

    ## The site function:
    # Linear component
    if(Vs30<= vc[ip]){
      ln_Flin = coff_c[ip] * log(Vs30 / v_ref)
    }else{
      ln_Flin = coff_c[ip] * log(vc[ip]/ v_ref)
    }

    # Nonlinear component
    f2 <- f4[ip]*(exp(f5[ip]*(min(Vs30,760)-360))-exp(f5[ip]*(760-360)))

    ln_Fnlin <- f1+f2*log((PGA_r+f3)/f3)

    # Effect of basin depth
    if(z1 != 999){
      if(region == 1){  # if in California
        mu_z1 <- exp(-7.15/4*log((Vs30^4+570.94^4)/(1360^4+570.94^4)))/1000
      }else if(region == 2){  # if in Japan
        mu_z1 <- exp(-5.23/2*log((Vs30^2+412.39^2)/(1360^2+412.39^2)))/1000
      }else{
        mu_z1 <- exp(-7.15/4*log((Vs30^4+570.94^4)/(1360^4+570.94^4)))/1000
      }
      dz1 <- z1-mu_z1
    }else{
      dz1 <- 0
    }

    if(z1 != 999){
      if(period[ip] < 0.65){
        F_dz1 <- 0
      }else if(period[ip] >= 0.65 && abs(dz1) <= f7[ip]/f6[ip]){
        F_dz1 <- f6[ip]*dz1
      }else{
        F_dz1 <- f7[ip]
      }
    }else{
      F_dz1 <- 0
    }

    F_S <- ln_Flin+ln_Fnlin+F_dz1

    ln_Y=F_E + F_P + F_S
    med=exp(ln_Y)

  }else{

    F_S <- 0
    ln_y <- F_E + F_P + F_S
    med <- exp(ln_y)
    PGA_r <- med

    # Linear component
    ln_Flin <- 0
    # Nonlinear component
    f2 <- f4[ip]*(exp(f5[ip]*(min(Vs30,760)-360))-exp(f5[ip]*(760-360)))
    ln_Fnlin <- f1+f2*log((PGA_r+f3)/f3)

    # Effect of basin depth
    if(z1 != 999){
      if(region == 1){  # if in California
        mu_z1 <- exp(-7.15/4*log((Vs30^4+570.94^4)/(1360^4+570.94^4)))/1000
      }else if(region == 2){  # if in Japan
        mu_z1 <- exp(-5.23/2*log((Vs30^2+412.39^2)/(1360^2+412.39^2)))/1000
      }else{
        mu_z1 <- exp(-7.15/4*log((Vs30^4+570.94^4)/(1360^4+570.94^4)))/1000
      }
      dz1 <- z1-mu_z1
    }else{
      dz1 <- 0
    }

    if(z1 != 999){
      if(period[ip] < 0.65){
        F_dz1 <- 0
      }else if(period[ip] >= 0.65 && abs(dz1) <= f7[ip]/f6[ip]){
        F_dz1 <- f6[ip]*dz1
      }else{
        F_dz1 <- f7[ip]
      }
    }else{
      F_dz1 <- 0
    }

  }
  ## Aleatory - uncertainty function
  if(M <= 4.5){
    tau <- tau1[ip]
    phi_M <- phi1[ip]
  }else if(4.5 < M && M < 5.5){
    tau <- tau1[ip]+(tau2[ip]-tau1[ip])*(M-4.5)
    phi_M <- phi1[ip]+(phi2[ip]-phi1[ip])*(M-4.5)
  }else if(M >= 5.5){
    tau <- tau2[ip]
    phi_M <- phi2[ip]
  }

  if(Rjb <= R1[ip]){
    phi_MR <- phi_M
  }else if(R1[ip] < Rjb && Rjb <= R2[ip]){
    phi_MR = phi_M + dphiR[ip]*(log(Rjb/R1[ip])/log(R2[ip]/R1[ip]))
  }else if(Rjb > R2[ip]){
    phi_MR = phi_M + dphiR[ip]
  }

  if(Vs30 >= v2){
    phi_MRV = phi_MR
  }else if(v1 <= Vs30 && Vs30 <= v2){
    phi_MRV = phi_MR - dphiV[ip]*(log(v2/Vs30)/log(v2/v1))
  }else if(Vs30 <= v1){
    phi_MRV = phi_MR - dphiV[ip]
  }

  sig = sqrt(phi_MRV^2 + tau^2)


  ### results return
  if(return.type == 1){
    pgar <- list()
    pgar$med <- med
    pgar$sigma <- sig
    pgar$phi <- phi_MRV
    pgar$tau <- tau
    return(pgar)
  }else if(return.type == 2){
    effs <- list()
    effs$F_E <- F_E
    effs$F_P <- F_P
    effs$F_P_GS <- F_P_GS
    effs$F_P_AA <- F_P_AA
    effs$F_S <- F_S
    return(effs)
  }else if(return.type == 3){
    res <- list()
    res$r <- r
    res$ln_Flin <- ln_Flin
    res$ln_Fnlin <- ln_Fnlin
    res$f2 <- f2
    res$F_dz1 <- F_dz1
    res$PGAr <- PGA_r
    return(res)
  }
}


#' The subroutine of CA subregional path correction
#'
#' This is a subroutine of CA subregional path correction
#' @param site_lon The longitude of site (in WGS84 degree).
#' @param site_lat The latitude of site (in WGS84 degree).
#' @param clst_lon Longitude (in WGS84 degrees) of closest point on the surface projection of the rupture surface to the site.
#' @param clst_lat Latitude (in WGS84 degrees) of closest point on the surface projection of the rupture surface to the site.
#' @return A list of three elements: 1) the site-to-source distance;
#' 2) the distances of the path within each CA subregion; and 3) the CA subregion index where the source is.
#' @examples CA_subreg_path_calc(site_lon = -118.8713, site_lat = 35.80407, clst_lon = -116.9274, clst_lat = 33.6204)
#' @importFrom sf st_linestring st_sfc st_length st_point st_within st_intersection
#' @export
CA_subreg_path_calc <- function(site_lon, site_lat, clst_lon, clst_lat) {

  site_source_path <- st_sfc(st_linestring(x = matrix(c(site_lon, clst_lon, site_lat, clst_lat), nrow = 2)),
                             crs = 4326)

  source_point <- st_sfc(st_point(x = c(clst_lon, clst_lat)), crs = 4326)

  source_region_idx <- st_within(source_point, CA_subreg)

  site_source_dist <- as.numeric(st_length(site_source_path))

  sub_dist <- rep(0, length(CA_subreg$geometry))

  for (i in 1:length(CA_subreg$geometry)) {
    tmp <- st_intersection(site_source_path, CA_subreg$geometry[i])
    if (length(tmp) > 0) {
      sub_dist[i] <- as.numeric(st_length(tmp))
    }
  }

  res <- list()

  res$site_source_dist <- site_source_dist

  res$sub_dist <- sub_dist

  res$source_region_idx <- ifelse(length(source_region_idx) == 0, NA, source_region_idx[[1]])

  return(res)
}



