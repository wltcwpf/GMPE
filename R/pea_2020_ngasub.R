#' The GMPE of Parker at al 2020
#'
#' This function calculates the ground motion median values and standard deviations
#' @param M Moment magnitude, a numeric value
#' @param T Period (sec); Use Period = -1 for PGV computation.
#' Use 1000 (by default) for output the array of Sa with original Pea et al pre-defined periods
#' @param Rrup Closest distance (km) to the ruptured plane
#' @param VS30 Shear wave velocity averaged over top 30 m (in m/s)
#' @param event.type The indicator of subduction type: 0 for interface; 1 for slab
#' @param region A string that corresponds to options in the DatabaseRegion column of the flatfile,
#' plus global. Must be a string in any of: "Global", "Alaska", "Cascadia", "CAM", "Japan", "SA" or "Taiwan".
#' If no matches, default will be global model.
#' @param saturation.region A string that corresponds to regions defined by C. Ji and R. Archuleta (2018):
#' "Global", "Aleutian","Alaska","Cascadia","CAM_S", "CAM_N", "Japan_Pac","Japan_Phi",
#' "SA_N","SA_S", "Taiwan_W","Taiwan_E". If unknown, use "Global".
#' @param Z2.5 The basin depth to the layer where the shear wave is 2.5 km/s.
#' Only used if DatabaseRegion == "Japan" or "Cascadia".
#' Can also specify "default" to get no basin term (e.g. delta_Z2.5 = 0)
#' @param basin Flag of baisn type. It is only used if DatabaseRegion == "Cascadia".
#' Value can be 0, 1, or 2,
#' where 0 == having an estimate of Z2.5 outside mapped basin, 1 == Seattle basin,
#' and 2 == other mapped basin (Tacoma, Everett, Georgia, etc.)
#' @param hypocentral.depth The depth of hypocenter (km).
#' To use Ztor value to estimate hypocentral depth, see Ch. 4.3.3/Eqs. 4.16 & 4.17
#' of Parker et al. (2020) PEER report
#' @param Ts A single period ranges from 0.01-10s for epistemic uncertainty.
#' Ts = NA will take the same periods as T.
#' @param Reg A string to specify region for epistemic uncertainty. Options include
#' "Global", "Alaska", "Aleutian", "Cascadia", "CAM_N", "CAM_S", "Japan_Pac", "Japan_Phi",
#' "SA_N", "SA_S", "Taiwan".
#' @param coefficients1 The model coefficients for interface. Use internal saved data, no need to specify.
#' @param coefficients2 The model coefficients for slab. Use internal saved data, no need to specify.
#' @param coefficients3 The model coefficients for aleatory. Use internal saved data, no need to specify.
#' @param coefficients4_IF The model coefficients for epistemic for interface. Use internal saved data, no need to specify.
#' @param coefficients4_Slab The model coefficients for epistemic for slab Use internal saved data, no need to specify.
#' @return A dataframe of 15 variables. 1. period, the same as input; 2. VS30, the same as input;
#' 3. Rrup, the same as input; 4. Tau, even-to-event Std; 5. Phi_total, total within-event Std;
#' 6. Phi_SS, single station Std; 7. PhiS2S, site-to-site Std;
#' 8. Phi_partitioned, remaining Std - Phi_total excludes Phi_SS and PhiS2S;
#' 9. Sigma_total, total Std; 10. Sigma_partitioned, Sigma_total excludes Tau and Phi_total;
#' 11. Mu, median prediction (g - PSA/PGA, cm/s - PGV); 12. Flin, linear amplification;
#' 13. Fnl, nonlinear amplification; 14. Rock, median prediction on reference rock site (g - PSA/PGA,
#' cm/s - PGV); 15. Epi, epistemic uncertainty
#' @examples pea_2020_ngasub(M = 7.8, T = 1000, Rrup = 50, VS30 = 350, event.type = 1,
#' hypocentral.depth = 45)
#'
#' pea_2020_ngasub(M = 8.5, T = c(0.035, 6), Rrup = 50, VS30 = 450, event.type = 0,
#' region = "Japan", saturation.region = "Japan_Pac", Z2.5 = "default", hypocentral.depth = 75,
#' Ts = NA, Reg = "Japan_Pac")
#' @references Parker, GA, JP Stewart, DM Boore, GM Atkinson, B Hassani (2020).
#' NGA-Subduction global ground-motion models with regional adjustment factors,
#' PEER Report 2020/03, Pacific Earthquake Engineering Research Center,
#' UC Berkeley (Center Headquarters), 131 pages.
#' @export
#' @importFrom stats approx pchisq
pea_2020_ngasub <- function(M, T = 1000, Rrup, VS30, event.type, region = "Global",
                            saturation.region = "Global", Z2.5 = "default", basin = 0,
                            hypocentral.depth, Ts = NA, Reg = "Global",
                            coefficients1 = pea_2020_coeffs1,
                            coefficients2 = pea_2020_coeffs2, coefficients3 = pea_2020_coeffs3,
                            coefficients4_IF = pea_2020_coeffs4_Interf,
                            coefficients4_Slab = pea_2020_coeffs4_Slab)
{

  pre_defined_periods <- c(-1,0,0.01,0.02,0.025,0.03,0.04,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,
                           0.75,1,1.5,2,2.5,3,4,5,7.5,10)


  if (length(T) == 1 & T[1] == 1000) {

    period_list <- pre_defined_periods

  } else {

    period_list <- T

  }

  for (i_p in 1:length(period_list)) {

    period <- period_list[i_p]

    if( event.type == 0 ){

      if( length( which( abs( pre_defined_periods - period ) < 0.0001 ) ) == 0 ){ # The user defined period requires interpolation

        idx <- which( pre_defined_periods < period )

        period_low = max( pre_defined_periods[ idx ] )

        idx <- which(pre_defined_periods > period )

        period_high = min( pre_defined_periods[ idx ] )

        # comment: low
        Res_sigma_low <- Aleatory_Function(period_low, Rrup, VS30, coefficients3)

        Tmp_low <- GMM_at_VS30_IF_v4(event.type,region,saturation.region,Rrup,M,period_low,VS30,Z2.5,basin,coefficients1)

        Res_sigma_low$Mu <- Tmp_low$mu

        Res_sigma_low$Flin <- Tmp_low$Flin

        Res_sigma_low$Fnl <- Tmp_low$Fnl

        Res_sigma_low$Rock <- Tmp_low$Rock

        Res_sigma_low$Epi <- PSHAB_Epi(event.type,Reg,period_low,coefficients4_IF)

        # comment: high
        Res_sigma_high <- Aleatory_Function(period_high, Rrup, VS30, coefficients3)

        Tmp_high <- GMM_at_VS30_IF_v4(event.type,region,saturation.region,Rrup,M,period_high,VS30,Z2.5,basin,coefficients1)

        Res_sigma_high$Mu <- Tmp_high$mu

        Res_sigma_high$Flin <- Tmp_high$Flin

        Res_sigma_high$Fnl <- Tmp_high$Fnl

        Res_sigma_high$Rock <- Tmp_high$Rock

        Res_sigma_high$Epi <- PSHAB_Epi(event.type,Reg,period_high,coefficients4_IF)

        #comment: interpolation
        periods <- c(  period_low ,  period_high  )

        Res_sigma <- Res_sigma_low

        Res_sigma$period <- period

        Res_sigma$Tau <- approx( periods, c( Res_sigma_low$Tau, Res_sigma_high$Tau ),  period, rule = 2 )$y

        Res_sigma$Phi_total <- approx( periods, c( Res_sigma_low$Phi_total, Res_sigma_high$Phi_total ),  period, rule = 2 )$y

        Res_sigma$Phi_SS <- approx( periods, c( Res_sigma_low$Phi_SS, Res_sigma_high$Phi_SS ), period, rule = 2 )$y

        Res_sigma$PhiS2S <- approx( periods, c( Res_sigma_low$PhiS2S, Res_sigma_high$PhiS2S ), period, rule = 2 )$y

        Res_sigma$Phi_partitioned <- approx( periods, c( Res_sigma_low$Phi_partitioned, Res_sigma_high$Phi_partitioned ), period, rule = 2 )$y

        Res_sigma$Sigma_total <- approx( periods, c( Res_sigma_low$Sigma_total, Res_sigma_high$Sigma_total ), period, rule = 2 )$y

        Res_sigma$Sigma_partitioned <- approx( periods, c( Res_sigma_low$Sigma_partitioned, Res_sigma_high$Sigma_partitioned ), period, rule = 2 )$y

        Res_sigma$Mu <- approx( periods, c( Res_sigma_low$Mu, Res_sigma_high$Mu ), period, rule = 2 )$y

        Res_sigma$Flin <- approx( periods, c( Res_sigma_low$Flin, Res_sigma_high$Flin ), period, rule = 2 )$y

        Res_sigma$Fnl <- approx( periods, c( Res_sigma_low$Fnl, Res_sigma_high$Fnl ), period, rule = 2 )$y

        Res_sigma$Rock <- approx( periods, c( Res_sigma_low$Rock, Res_sigma_high$Rock ), period, rule = 2 )$y

        Res_sigma$Epi <- approx( periods, c( Res_sigma_low$Epi, Res_sigma_high$Epi ), period, rule = 2 )$y

      }else{  # no need interpolation

        period <- pre_defined_periods[ which( abs( pre_defined_periods - period ) < 0.0001 )[ 1 ] ]

        Res_sigma <- Aleatory_Function(period, Rrup, VS30, coefficients3)

        Tmp <- GMM_at_VS30_IF_v4(event.type,region,saturation.region,Rrup,M,period,VS30,Z2.5,basin,coefficients1)

        Res_sigma$Mu <- Tmp$mu

        Res_sigma$Flin <- Tmp$Flin

        Res_sigma$Fnl <- Tmp$Fnl

        Res_sigma$Rock <- Tmp$Rock

        # use period not Ts
        # Res_sigma$Epi <- PSHAB_Epi(event.type,Reg,Ts,coefficients4_IF)
        Res_sigma$Epi <- PSHAB_Epi(event.type,Reg,period,coefficients4_IF)

      }

    }else if( event.type == 1 ){

      if( length( which( abs( pre_defined_periods - period ) < 0.0001 ) ) == 0 ){ # The user defined period requires interpolation

        idx <- which( pre_defined_periods < period )

        period_low = max( pre_defined_periods[ idx ] )

        idx <- which(pre_defined_periods > period )

        period_high = min( pre_defined_periods[ idx ] )

        # comment: low
        Res_sigma_low <- Aleatory_Function(period_low, Rrup, VS30, coefficients3)

        Tmp_low <- GMM_at_VS30_Slab_v4(event.type,region,saturation.region,Rrup,M,hypocentral.depth,period_low,VS30, Z2.5,basin,coefficients2)

        Res_sigma_low$Mu <- Tmp_low$mu

        Res_sigma_low$Flin <- Tmp_low$Flin

        Res_sigma_low$Fnl <- Tmp_low$Fnl

        Res_sigma_low$Rock <- Tmp_low$Rock

        Res_sigma_low$Epi <- PSHAB_Epi(event.type,Reg,period_low,coefficients4_Slab)

        # comment: high
        Res_sigma_high <- Aleatory_Function(period_high, Rrup, VS30, coefficients3)

        Tmp_high <- GMM_at_VS30_Slab_v4(event.type,region,saturation.region,Rrup,M,hypocentral.depth,period_high,VS30, Z2.5,basin,coefficients2)

        Res_sigma_high$Mu <- Tmp_high$mu

        Res_sigma_high$Flin <- Tmp_high$Flin

        Res_sigma_high$Fnl <- Tmp_high$Fnl

        Res_sigma_high$Rock <- Tmp_high$Rock

        Res_sigma_high$Epi <- PSHAB_Epi(event.type,Reg,period_high,coefficients4_Slab)

        #comment: interpolation
        periods <- c( period_low, period_high )

        Res_sigma <- Res_sigma_low

        Res_sigma$period <- period

        Res_sigma$Tau <- approx( periods, c( Res_sigma_low$Tau, Res_sigma_high$Tau ), period, rule = 2 )$y

        Res_sigma$Phi_total <- approx( periods, c( Res_sigma_low$Phi_total, Res_sigma_high$Phi_total ), period, rule = 2 )$y

        Res_sigma$Phi_SS <- approx( periods, c( Res_sigma_low$Phi_SS, Res_sigma_high$Phi_SS ), period, rule = 2 )$y

        Res_sigma$PhiS2S <- approx( periods, c( Res_sigma_low$PhiS2S, Res_sigma_high$PhiS2S ), period, rule = 2 )$y

        Res_sigma$Phi_partitioned <- approx( periods, c( Res_sigma_low$Phi_partitioned, Res_sigma_high$Phi_partitioned ), period, rule = 2 )$y

        Res_sigma$Sigma_total <- approx( periods, c( Res_sigma_low$Sigma_total, Res_sigma_high$Sigma_total ), period, rule = 2 )$y

        Res_sigma$Sigma_partitioned <- approx( periods, c( Res_sigma_low$Sigma_partitioned, Res_sigma_high$Sigma_partitioned ), period, rule = 2 )$y

        Res_sigma$Mu <- approx( periods, c( Res_sigma_low$Mu, Res_sigma_high$Mu ), period, rule = 2 )$y

        Res_sigma$Flin <- approx( periods, c( Res_sigma_low$Flin, Res_sigma_high$Flin ), period, rule = 2 )$y

        Res_sigma$Fnl <- approx( periods, c( Res_sigma_low$Fnl, Res_sigma_high$Fnl ), period, rule = 2 )$y

        Res_sigma$Rock <- approx( periods, c( Res_sigma_low$Rock, Res_sigma_high$Rock ), period, rule = 2 )$y

        Res_sigma$Epi <- approx( periods, c( Res_sigma_low$Epi, Res_sigma_high$Epi ), period, rule = 2 )$y

      }else{  # no need interpolation

        period <- pre_defined_periods[ which( abs( pre_defined_periods - period ) < 0.0001 )[ 1 ] ]

        Res_sigma <- Aleatory_Function(period, Rrup, VS30, coefficients3)

        Tmp <- GMM_at_VS30_Slab_v4(event.type,region,saturation.region,Rrup,M,hypocentral.depth,period,VS30, Z2.5,basin,coefficients2)

        Res_sigma$Mu <- Tmp$mu

        Res_sigma$Flin <- Tmp$Flin

        Res_sigma$Fnl <- Tmp$Fnl

        Res_sigma$Rock <- Tmp$Rock

        # use period not Ts
        # Res_sigma$Epi <- PSHAB_Epi(event.type,Reg,Ts,coefficients4_Slab)
        Res_sigma$Epi <- PSHAB_Epi(event.type,Reg,period,coefficients4_Slab)
      }
    }

    Res_sigma$Mu <- exp(Res_sigma$Mu)

    Res_sigma$Rock <- exp(Res_sigma$Rock)

    ifelse (i_p == 1, Res <- Res_sigma,
            Res <- rbind(Res, Res_sigma))
  }

  return( Res )
}



## helper function
erf <- function (x)
{
  pchisq(2 * x^2, 1) * sign(x)
}


## subroutines
Aleatory_Function <- function(period, Rrup, VS30, coefficients) {
  #Want to generate tau, phi, and total sigma computed from both total and partitioned phi models


  # Coefficients ------------------------------------------------------------
  # coefficients <- read.csv("Table_E3_Aleatory_Coefficients.csv", header = T)
  # coefficients.T <- subset(coefficients, Period..s. == period)
  coefficients.T <- coefficients[coefficients$Period..s. == period, ]

  #Define period-independent coefficients for phi models
  V1 <- 200
  V2 <- 500
  R1 <- 200
  R2 <- 500

  R3 = 200
  R4 = 500
  R5 = 500
  R6 = 800

  # Total Phi ---------------------------------------------------------------

  PhiV <- as.numeric()
  PhiR <- as.numeric()

  for(i in 1:length(VS30)){

    if(Rrup[i] <= R1){
      PhiR[i] <- coefficients.T$Phi21
    }else if(Rrup[i] >= R2){
      PhiR[i] <- coefficients.T$Phi22
    }else{
      PhiR[i] <- ((coefficients.T$Phi22-coefficients.T$Phi21)/(log(R2)-log(R1)))*(log(Rrup[i])-log(R1)) + coefficients.T$Phi21
    }
    if(VS30[i] <= V1){
      PhiV[i] <- coefficients.T$Phi2V*(log(R2/max(R1,min(R2,Rrup[i])))/log(R2/R1))
    }else if(VS30[i] >= V2){
      PhiV[i] <- 0
    }else{
      PhiV[i] <- coefficients.T$Phi2V*((log(V2/min(V2,VS30[i])))/(log(V2/V1)))*(log(R2/max(R1,min(R2,Rrup[i])))/(log(R2/R1)))
    }

  }

  Phi_tot = sqrt(PhiV + PhiR)

  # Partitioned Phi ---------------------------------------------------------

  PhiS2SV <- as.numeric()
  PhiSSR <- as.numeric()
  PhiSSV <- as.numeric()

  for(i in 1:length(VS30)){
    if(VS30[i] <= V1){
      PhiS2SV[i] <- coefficients.T$Phi2S2S.0 + (coefficients.T$a1*log(V1/coefficients.T$VM))*(log(R4/max(R3,min(R4,Rrup[i])))/log(R4/R3))
    }else if(VS30[i] >= V2){
      PhiS2SV[i] <- coefficients.T$Phi2S2S.0 + (coefficients.T$a1*log(V2/coefficients.T$VM))
    }else if(VS30[i] > V1 & VS30[i] < coefficients.T$VM){
      PhiS2SV[i] <- coefficients.T$Phi2S2S.0 + (coefficients.T$a1*log(VS30[i]/coefficients.T$VM))*(log(R4/max(R3,min(R4,Rrup[i])))/log(R4/R3))
    }else if(VS30[i] >=coefficients.T$VM & VS30[i] < V2){
      PhiS2SV[i] <- coefficients.T$Phi2S2S.0 + (coefficients.T$a1*log(VS30[i]/coefficients.T$VM))
    }
  }

  for(i in 1:length(VS30)){
    if(VS30[i] <= V1){
      PhiSSV[i] <- coefficients.T$a2*log(V1/coefficients.T$VM)*(log(R4/max(R3,min(R4,Rrup[i])))/log(R4/R3))
    }else if(VS30[i] >= V2){
      PhiSSV[i] <- coefficients.T$a2*log(V2/coefficients.T$VM)
    }else if(VS30[i] > V1 & VS30[i] < coefficients.T$VM){
      PhiSSV[i] <- coefficients.T$a2*log(VS30[i]/coefficients.T$VM)*(log(R4/max(R3,min(R4,Rrup[i])))/log(R4/R3))
    }else if(VS30[i] >=coefficients.T$VM & VS30[i] < V2){
      PhiSSV[i] <- coefficients.T$a2*log(VS30[i]/coefficients.T$VM)
    }

  }

  for(j in 1:length(Rrup)){
    if(Rrup[i] <= R5){
      PhiSSR[i] <- coefficients.T$Phi2SS.1
    }else if(Rrup[i] >= R6){
      PhiSSR[i] <- coefficients.T$Phi2SS.2
    }else{
      PhiSSR[i] <- ((coefficients.T$Phi2SS.2-coefficients.T$Phi2SS.1)/(log(R6)-log(R5)))*(log(Rrup[i])-log(R5)) + coefficients.T$Phi2SS.1
    }

  }


  PhiC = sqrt(PhiS2SV + PhiSSR + PhiSSV)

  # Define output matrix ----------------------------------------------------

  sigma.matrix <- data.frame(period = period,
                             VS30 = VS30,
                             Rrup = Rrup,
                             Tau = coefficients.T$Tau,
                             Phi_total = Phi_tot,
                             Phi_SS = sqrt(PhiSSR + PhiSSV),
                             PhiS2S = sqrt(PhiS2SV),
                             Phi_partitioned = PhiC)

  sigma.matrix$Sigma_total <- sqrt(sigma.matrix$Tau^2 + sigma.matrix$Phi_total^2)
  sigma.matrix$Sigma_partitioned <- sqrt(sigma.matrix$Tau^2 + sigma.matrix$Phi_partitioned^2)

  return(sigma.matrix)
}

GMM_at_760_IF_v4 <- function(event.type,region,saturation.region,Rrup,M,period,coefficients){

  #Grace Parker
  #Modified February 26 to expand comments
  #Modified March 25, 2020 to call consolidated coefficient tables

  # Input Parameters --------------------------------------------------------

  #Event type: 0 == interface, 1 == slab

  #region corresponds to options in the DatabaseRegion column of the flatfile, plus global. Must be a string. If no matches, default will be global model:
  # "global", "Alaska", "Cascadia", "CAM", "Japan", "SA" or "Taiwan"

  #Saturation Region corresponds to regions defined by by C. Ji and R. Archuleta (2018):
  # "global", "Aleutian","Alaska","Cascadia","CAM_S", "CAM_N", "Japan_Pac","Japan_Phi","SA_N","SA_S", "Taiwan_W","Taiwan_E"

  # Rrup is number in kilometers

  # period can be: (-1,0,0.01,0.02,0.025,0.03,0.04,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1,1.5,2,2.5,3,4,5,7.5,10)
  # where -1 == PGV and 0 == PGA

  # Other pertinent information ---------------------------------------------
  # Coefficient files must be in the active working directory

  # This function has no site term. Can only estimate ground motion at the reference condition VS30 = 760m/s

  # The output is the desired median model prediction in LN units
  # Take the exponential to get PGA, PSA in g or the PGV in cm/s

  #This model does not make predictions for New Zealand. We recommend using the global model for that case.


  if(event.type == 1 | event.type == 5){
    stop("This function is only for IF")
  }
  # Import Coefficient Table -----------------------------------------

  # coefficients <- read.csv("Table_E1_Interface_Coefficients_051920.csv", header =  T)

  #Isolate desired period
  # coefficients.T <- subset(coefficients,Period..s. == period)
  coefficients.T <- coefficients[coefficients$Period..s. == period, ]

  #Define Mb
  if(saturation.region == "global"){
    Mb <- 7.9
  }else{
    Interface.saturation.regions <- data.frame(SBZ = c("Aleutian","Alaska","-999","Cascadia","CAM_S", "CAM_N", "Japan_Pac","Japan_Phi", "New_Zealand_N","New_Zealand_S","SA_N","SA_S", "Taiwan_W","Taiwan_E"),
                                               Mb = c(8,8.6,NA,7.7,7.4,7.4,8.5,7.7,NA,NA,8.5,8.6,7.1,7.1))
    # Mb <- as.numeric(subset(Interface.saturation.regions, SBZ == as.character(saturation.region), select = Mb))
    Mb <- as.numeric(Interface.saturation.regions$Mb[Interface.saturation.regions$SBZ == as.character(saturation.region)])
    # if(is.na(Mb)){
    #   Mb <- 7.9
    # }
    if(length(Mb) == 0){
      Mb <- 7.9
    }
  }

  # Constant ----------------------------------------------------------------

  if(region == "global"){
    # c0 <- as.numeric(subset(coefficients.T, select = Global_c0))
    c0 <- as.numeric(coefficients.T$Global_c0)
  }else{
    c0 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(saturation.region,"_c0",sep=""))))))
  }

  # Path Term ---------------------------------------------------------------
  h <- 10^(-0.82 + 0.252*M)
  Rref <- sqrt(1 + h^2)
  R <- sqrt(Rrup^2 + h^2)
  LogR <- log(R)
  R_Rref <- log(R/Rref)

  #Need  to isolate regional anelastic coefficient, a0
  if(region == "global" | region == "Cascadia"){
    # a0 <- as.numeric(subset(coefficients.T, select = c(Global_a0)))
    a0 <- as.numeric(coefficients.T$Global_a0)
  }else{
    a0 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(region,"_a0",sep=""))))))
  }

  if(is.na(a0)){
    # a0 <- as.numeric(subset(coefficients.T, select = c(Global_a0)))
    a0 <- as.numeric(coefficients.T$Global_a0)
  }

  Fp <- as.numeric(coefficients.T$c1*LogR + (coefficients.T$b4*M)*R_Rref + a0*R)


  # Magnitude Scaling -------------------------------------------------------

  #compute M-scaling term
  func1 <- function(M,Mb){ifelse(M <= Mb,M-Mb,0)}
  func2 <- function(M,Mb){ifelse(M > Mb,M-Mb,0)}
  func3 <- function(M,Mb){ifelse(M <= Mb, (M-Mb)^2, 0)}
  Fm <- coefficients.T$c4*func1(M,Mb) + coefficients.T$c6*func2(M,Mb) + coefficients.T$c5*func3(M,Mb)


  # Add it all up! ----------------------------------------------------------

  mu <-   Fp + Fm + c0
  return(c(mu))
}

GMM_at_760_Slab_v4 <- function(event.type,region,saturation.region,Rrup,M,hypocentral.depth,period,coefficients){

  #Grace Parker
  #Modified February 26 to expand comments
  #Modified March 25, 2020 to call coefficients from master table

  # Input Parameters --------------------------------------------------------

  #Event type: 0 == interface, 1 == slab

  #region corresponds to options in the DatabaseRegion column of the flatfile, plus global. Must be a string. If no matches, default will be global model:
  # "global", "Alaska", "Cascadia", "CAM", "Japan", "SA" or "Taiwan"

  #Saturation Region corresponds to regions defined by R. Archuleta and C. Ji:
  # "global", "Aleutian","Alaska","Cascadia","CAM_S", "CAM_N", "Japan_Pac","Japan_Phi","SA_N","SA_S", "Taiwan_W","Taiwan_E"

  # Rrup is number in kilometers

  #Hypocentral depth in km.
  #To use Ztor value to estimate hypocentral depth, see Ch. 4.3.3/Eqs. 4.16 & 4.17  of Parker et al. (2020) PEER report

  # period can be: (-1,0,0.01,0.02,0.025,0.03,0.04,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1,1.5,2,2.5,3,4,5,7.5,10)
  # where -1 == PGV and 0 == PGA

  # Other pertinent information ---------------------------------------------
  #Coefficient files must be in the active working directory

  # This function has no site term. Can only estimate ground motion at the reference condition VS30 = 760m/s

  # The output is the desired median model prediction in LN units
  # Take the exponential to get PGA, PSA in g or the PGV in cm/s

  #This model does not make predictions for New Zealand. We recommend using the global model for that case.

  #Function to compute GMM predictions at 760m/s for slab and interface


  if(event.type == 0){
    stop("This function is only for slab")
  }

  # Import Coefficient Table -----------------------------------------

  # coefficients <- read.csv("Table_E2_Slab_Coefficients_060720.csv", header =  T, skip = 1)

  #Isolate desired period
  # coefficients.T <- subset(coefficients,Period..s. == period)
  coefficients.T <- coefficients[coefficients$Period..s. == period, ]

  #Define mb based on Archuleta and Ji (2019)
  if(saturation.region == "global"){
    Mb <- 7.6
  }else{
    slab.saturation.regions <- data.frame(SBZ = c("Aleutian","Alaska","-999","Cascadia","CAM_S", "CAM_N", "Japan_Pac","Japan_Phi", "New_Zealand_N","New_Zealand_S","SA_N","SA_S", "Taiwan_W","Taiwan_E"),
                                          Mb = c(7.98,7.2,NA,7.2,7.6,7.4,7.65,7.55,NA,NA,7.3,7.25,7.7,7.7))
    # Mb <- as.numeric(subset(slab.saturation.regions, SBZ == as.character(saturation.region), select = Mb))
    Mb <- as.numeric(slab.saturation.regions$Mb[slab.saturation.regions$SBZ == as.character(saturation.region)])
  }
  # if(is.na(Mb)){
  #   Mb <- 7.6
  # }
  if(length(Mb) == 0){
    Mb <- 7.6
  }
  # Constant ----------------------------------------------------------------

  #Isolate constant
  if(region == "global"){
    # c0 <- as.numeric(subset(coefficients.T,select = Global_c0))
    c0 <- as.numeric(coefficients.T$Global_c0)
  }else if(region == "Alaska" | region == "SA"){
    c0 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(saturation.region,"_c0",sep=""))))))
  }else{
    c0 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(region,"_c0",sep=""))))))
  }

  # Path Term ---------------------------------------------------------------

  #near-source saturation
  if(M <= Mb){
    m <- (log10(35)-log10(3.12))/(Mb-4)
    h <- 10^(m*(M-Mb) + log10(35))
  }else{
    h <- 35
  }

  Rref <- sqrt(1 + h^2)
  R <- sqrt(Rrup^2 + h^2)
  LogR <- log(R)
  R_Rref <- log(R/Rref)

  #Need  to isolate regional anelastic coefficient, a0
  if(region == "global"){
    a0 <- coefficients.T$Global_a0
  }else{
    a0 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(region,"_a0",sep=""))))))
  }
  if(is.na(a0)){
    a0 <- coefficients.T$Global_a0
  }

  Fp <- as.numeric(coefficients.T$c1*LogR + (coefficients.T$b4*M)*R_Rref + a0*R)

  # Magnitude Scaling -------------------------------------------------------

  func1 <- function(M,Mb){ifelse(M <= Mb,M-Mb,0)}
  func2 <- function(M,Mb){ifelse(M > Mb,M-Mb,0)}
  func3 <- function(M,Mb){ifelse(M <= Mb, (M-Mb)^2, 0)}

  Fm <- coefficients.T$c4*func1(M,Mb) + coefficients.T$c6*func2(M,Mb) + coefficients.T$c5*func3(M,Mb)


  # Source Depth Scaling ----------------------------------------------------

  if(hypocentral.depth >= coefficients.T$db..km.){
    Fd <-coefficients.T$d
  }else if(hypocentral.depth <= 20){
    Fd <- coefficients.T$m*(20-coefficients.T$db..km.) + coefficients.T$d
  }else{
    Fd <- coefficients.T$m*(hypocentral.depth - coefficients.T$db..km.) + coefficients.T$d
  }


  mu <- Fp + Fm + Fd + c0
  return(mu)
}

GMM_at_VS30_IF_v4 <- function(event.type,region,saturation.region,Rrup,M,period,VS30,Z2.5,basin,coefficients){

  #Grace Parker
  #Modified February 26, 2020, to expand comments
  #Modified March 25, 2020, to to call consolidated coefficient table

  # Input Parameters --------------------------------------------------------

  #Event type: 0 == interface, 1 == slab

  #region corresponds to options in the DatabaseRegion column of the flatfile, plus global. Must be a string. If no matches, default will be global model:
  # "global", "Alaska", "Cascadia", "CAM", "Japan", "SA" or "Taiwan"

  #Saturation Region corresponds to regions defined by C. Ji and R. Archuleta (2018):
  # "global", "Aleutian","Alaska","Cascadia","CAM_S", "CAM_N", "Japan_Pac","Japan_Phi","SA_N","SA_S", "Taiwan_W","Taiwan_E"

  # Rrup is number in kilometers

  #Hypocentral depth in km.
  #To use Ztor value to estimate hypocentral depth, see Ch. 4.3.3/Eqs. 4.16 & 4.17  of Parker et al. (2020) PEER report

  # period can be: (-1,0,0.01,0.02,0.025,0.03,0.04,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1,1.5,2,2.5,3,4,5,7.5,10)
  # where -1 == PGV and 0 == PGA

  #VS30 in units m/s

  #Z2.5 in units m. Only used if DatabaseRegion == "Japan" or "Cascadia".
  #Can also specify "default" to get no basin term (e.g. delta_Z2.5 = 0)

  #basin is only used if DatabaseRegion == "Cascadia". Value can be 0, 1, or 2,
  #where 0 == having an estimate of Z2.5 outside mapped basin, 1 == Seattle basin,
  #and 0 == other mapped basin (Tacoma, Everett, Georgia, etc.)

  # Other pertinent information ---------------------------------------------

  # Coefficient files must be in the active working directory or have the path specified

  # "GMM_at_VS30_IF_v4.R" calls function "GMM_at_760_IF_v4.R" to compute PGAr in the nonlinear site term.
  #This function must be in the R environment else an error will occur.

  # The output is the desired median model prediction in LN units
  # Take the exponential to get PGA, PSA in g or the PGV in cm/s

  #This model does not make predictions for New Zealand. We recommend using the global model for that case.


  if(event.type == 1 | event.type == 5){
    stop("This function is only for IF")
  }

  # Import Coefficient Table -----------------------------------------

  # coefficients <- read.csv("Table_E1_Interface_Coefficients_051920.csv", header =  T)

  #Isolate desired period
  # coefficients.T <- subset(coefficients,Period..s. == period)
  coefficients.T <- coefficients[coefficients$Period..s. == period, ]

  #Define Mb
  if(saturation.region == "global"){
    Mb <- 7.9
  }else{
    Interface.saturation.regions <- data.frame(SBZ = c("Aleutian","Alaska","-999","Cascadia","CAM_S", "CAM_N", "Japan_Pac","Japan_Phi", "New_Zealand_N","New_Zealand_S","SA_N","SA_S", "Taiwan_W","Taiwan_E"),
                                               Mb = c(8,8.6,NA,7.7,7.4,7.4,8.5,7.7,NA,NA,8.5,8.6,7.1,7.1))
    # Mb <- as.numeric(subset(Interface.saturation.regions, SBZ == as.character(saturation.region), select = Mb))
    Mb <- as.numeric(Interface.saturation.regions$Mb[Interface.saturation.regions$SBZ == as.character(saturation.region)])
    # if(is.na(Mb)){
    #   Mb <- 7.9
    # }
    if(length(Mb) == 0){
      Mb <- 7.9
    }
  }

  # Constant ----------------------------------------------------------------

  if(region == "global"){
    # c0 <- as.numeric(subset(coefficients.T, select = Global_c0))
    c0 <- as.numeric(coefficients.T$Global_c0)
  }else{
    c0 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(saturation.region,"_c0",sep=""))))))
  }

  # Path Term ---------------------------------------------------------------
  h <- 10^(-0.82 + 0.252*M)
  Rref <- sqrt(1 + h^2)
  R <- sqrt(Rrup^2 + h^2)
  LogR <- log(R)
  R_Rref <- log(R/Rref)

  #Need  to isolate regional anelastic coefficient, a0
  if(region == "global" | region == "Cascadia"){
    # a0 <- as.numeric(subset(coefficients.T, select = c(Global_a0)))
    a0 <- as.numeric(coefficients.T$Global_a0)
  }else{
    a0 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(region,"_a0",sep=""))))))
  }

  if(is.na(a0)){
    # a0 <- as.numeric(subset(coefficients.T, select = c(Global_a0)))
    a0 <- as.numeric(coefficients.T$Global_a0)
  }

  Fp <- as.numeric(coefficients.T$c1*LogR + (coefficients.T$b4*M)*R_Rref + a0*R)


  # Magnitude Scaling -------------------------------------------------------

  #compute M-scaling term
  func1 <- function(M,Mb){ifelse(M <= Mb,M-Mb,0)}
  func2 <- function(M,Mb){ifelse(M > Mb,M-Mb,0)}
  func3 <- function(M,Mb){ifelse(M <= Mb, (M-Mb)^2, 0)}
  Fm <- coefficients.T$c4*func1(M,Mb) + coefficients.T$c6*func2(M,Mb) + coefficients.T$c5*func3(M,Mb)

  # Linear Site Amplification ----------------------------------------------

  #Site Coefficients
  V1 <- coefficients.T$V1..m.s.
  V2 <- coefficients.T$V2..m.s.
  Vref <- coefficients.T$Vref..m.s.

  if(region == "global"| region == "CAM"){
    s2 <- coefficients.T$Global_s2
    s1 <- s2
  }else if(region =="Taiwan" | region == "Japan"){
    s2 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(region,"_s2",sep=""))))))
    s1 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(region,"_s1",sep=""))))))
  }else{
    s2 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(region,"_s2",sep=""))))))
    s1 <- s2
  }

  #Compute linear site term
  if(VS30 <= V1){
    Flin <- s1*log(VS30/V1) + s2*log(V1/Vref)
  }else if(VS30 <= V2){
    Flin <- s2*log(VS30/Vref)
  }else{
    Flin <- s2*log(V2/Vref)
  }

  # Nonlinear Site Term -----------------------------------------------------

  PGAr <- exp(GMM_at_760_IF_v4(event.type,region,saturation.region,Rrup,M,0,coefficients))
  f3 <- 0.05
  Vb <- 200
  Vref.Fnl <- 760

  if(period >= 3){
    Fnl = 0;
  }else{
    f2 <- coefficients.T$f4*(exp(coefficients.T$f5*(min(VS30,Vref.Fnl)-Vb))-exp(coefficients.T$f5*(Vref.Fnl-Vb)))
    Fnl = 0 + f2*log((PGAr+f3)/f3)
  }

  # Basin Term --------------------------------------------------------------

  #import coefficients

  if(Z2.5 == "default" | Z2.5 <= 0 | (region != "Japan" & region != "Cascadia")){
    Fb = 0
  }else{
    if(region == "Cascadia"){
      theta0 <- 3.94
      theta1 <- -0.42
      vmu <- 200
      vsig <- 0.2
      e1 <- coefficients.T$C_e1

      if(basin == 0){
        e3 <- coefficients.T$C_e3 + coefficients.T$del_None
        e2 <-  coefficients.T$C_e2 +  coefficients.T$del_None
      }else if(basin == 1){
        e3 <- coefficients.T$C_e3 + coefficients.T$del_Seattle
        e2 <-  coefficients.T$C_e2 +  coefficients.T$del_Seattle
      }else{
        e3 <- coefficients.T$C_e3
        e2 <-  coefficients.T$C_e2
      }

    }else if(region == "Japan"){
      theta0 <- 3.05
      theta1 <- -0.8
      vmu <- 500
      vsig <-0.33
      e3 <- coefficients.T$J_e3
      e2 <- coefficients.T$J_e2
      e1 <- coefficients.T$J_e1
    }


    Z2.5.pred <- 10^(theta0 + theta1*(1 + erf((log10(VS30) - log10(vmu))/(vsig*sqrt(2)))))
    delZ2.5 <- log(Z2.5) - log(Z2.5.pred)

    if(delZ2.5 <= (e1/e3)){
      Fb <- e1
    }else if(delZ2.5 >= (e2/e3)){
      Fb <- e2
    }else{
      Fb <- e3*delZ2.5
    }
  }



  # Add it all up! ----------------------------------------------------------

  mu <- Fp + Fnl + Fb + Fm + c0 + Flin
  res <- list()
  res$mu <- mu
  res$Flin <- Flin
  res$Fnl <- Fnl
  res$Rock <- mu - ( Flin + Fnl )
  return(res)
}

GMM_at_VS30_Slab_v4 <- function(event.type,region,saturation.region,Rrup,M,hypocentral.depth,period,VS30, Z2.5,basin,coefficients){

  #Grace Parker
  #Modified February 26 to expand comments
  #Modified March 25, 2020 to call consolidated coefficient table

  # Input Parameters --------------------------------------------------------

  #Event type: 0 == interface, 1 == slab

  #region corresponds to options in the DatabaseRegion column of the flatfile, plus global. Must be a string. If no matches, default will be global model:
  # "global", "Alaska", "Cascadia", "CAM", "Japan", "SA" or "Taiwan"

  #Saturation Region corresponds to regions defined by C. Ji and R. Archuleta (2018):
  # "global", "Aleutian","Alaska","Cascadia","CAM_S", "CAM_N", "Japan_Pac","Japan_Phi","SA_N","SA_S", "Taiwan_W","Taiwan_E"

  # Rrup is number in kilometers

  #Hypocentral depth in km.
  #To use Ztor value to estimate hypocentral depth, see Ch. 4.3.3/Eqs. 4.16 & 4.17  of Parker et al. (2020) PEER report

  # period can be: (-1,0,0.01,0.02,0.025,0.03,0.04,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1,1.5,2,2.5,3,4,5,7.5,10)
  # where -1 == PGV and 0 == PGA

  #VS30 in units m/s

  #Z2.5 in units m. Only used if DatabaseRegion == "Japan" or "Cascadia".
  #Can also specify "default" to get no basin term (e.g. delta_Z2.5 = 0)

  #basin is only used if DatabaseRegion == "Cascadia". Value can be 0, 1, or 2,
  #where 0 == having an estimate of Z2.5 outside mapped basin, 1 == Seattle basin,
  #and 0 == other mapped basin (Tacoma, Everett, Georgia, etc.)

  # Other pertinent information ---------------------------------------------

  #Coefficient files must be in the active working directory

  # "GMM_at_VS30_Slab_v4.R" calls function "GMM_at_760_Slab_v4.R" to compute PGAr in the nonlinear site term.
  #This function must be in the R environment else an error will occur.

  # The output is the desired median model prediction in LN units
  # Take the exponential to get PGA, PSA in g or  PGV in cm/s

  #Function to compute GMM predictions at various VS30s for slab

  if(event.type == 0){
    stop("This function is only for slab")
  }

  # Import Coefficient Table -----------------------------------------

  # coefficients <- read.csv("Table_E2_Slab_Coefficients_060720.csv", header =  T, skip = 1)

  #Isolate desired period
  # coefficients.T <- subset(coefficients, Period..s. == period)
  coefficients.T <- coefficients[coefficients$Period..s. == period, ]

  #Define mb based on Archuleta and Ji (2019)
  if(saturation.region == "global"){
    Mb <- 7.6
  }else{
    slab.saturation.regions <- data.frame(SBZ = c("Aleutian","Alaska","-999","Cascadia","CAM_S", "CAM_N", "Japan_Pac","Japan_Phi", "New_Zealand_N","New_Zealand_S","SA_N","SA_S", "Taiwan_W","Taiwan_E"),
                                          Mb = c(7.98,7.2,NA,7.2,7.6,7.4,7.65,7.55,NA,NA,7.3,7.25,7.7,7.7))
    # Mb <- as.numeric(subset(slab.saturation.regions, SBZ == as.character(saturation.region), select = Mb))
    Mb <- as.numeric(slab.saturation.regions$Mb[slab.saturation.regions$SBZ == as.character(saturation.region)])
  }
  # if(is.na(Mb)){
  #   Mb <- 7.6
  # }
  if(length(Mb) == 0){
    Mb <- 7.6
  }

  # Constant ----------------------------------------------------------------

  #Isolate constant
  if(region == "global"){
    # c0 <- as.numeric(subset(coefficients.T,select = Global_c0))
    c0 <- as.numeric(coefficients.T$Global_c0)
  }else if(region == "Alaska" | region == "SA"){
    c0 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(saturation.region,"_c0",sep=""))))))
  }else{
    c0 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(region,"_c0",sep=""))))))
  }

  # Path Term ---------------------------------------------------------------

  #near-source saturation
  if(M <= Mb){
    m <- (log10(35)-log10(3.12))/(Mb-4)
    h <- 10^(m*(M-Mb) + log10(35))
  }else{
    h <- 35
  }

  Rref <- sqrt(1 + h^2)
  R <- sqrt(Rrup^2 + h^2)
  LogR <- log(R)
  R_Rref <- log(R/Rref)

  #Need  to isolate regional anelastic coefficient, a0
  if(region == "global"){
    a0 <- coefficients.T$Global_a0
  }else{
    a0 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(region,"_a0",sep=""))))))
  }
  if(is.na(a0)){
    a0 <- coefficients.T$Global_a0
  }

  Fp <- as.numeric(coefficients.T$c1*LogR + (coefficients.T$b4*M)*R_Rref + a0*R)

  # Magnitude Scaling -------------------------------------------------------

  func1 <- function(M,Mb){ifelse(M <= Mb,M-Mb,0)}
  func2 <- function(M,Mb){ifelse(M > Mb,M-Mb,0)}
  func3 <- function(M,Mb){ifelse(M <= Mb, (M-Mb)^2, 0)}

  Fm <- coefficients.T$c4*func1(M,Mb) + coefficients.T$c6*func2(M,Mb) + coefficients.T$c5*func3(M,Mb)


  # Source Depth Scaling ----------------------------------------------------

  if(hypocentral.depth >= coefficients.T$db..km.){
    Fd <-coefficients.T$d
  }else if(hypocentral.depth <= 20){
    Fd <- coefficients.T$m*(20-coefficients.T$db..km.) + coefficients.T$d
  }else{
    Fd <- coefficients.T$m*(hypocentral.depth - coefficients.T$db..km.) + coefficients.T$d
  }

  # Linear Site Amplification ----------------------------------------------

  #Site Coefficients
  V1 <- coefficients.T$V1..m.s.
  V2 <- coefficients.T$V2..m.s.
  Vref <- coefficients.T$Vref..m.s.

  if(region == "global"| region == "CAM"){
    s2 <- coefficients.T$Global_s2
    s1 <- s2
  }else if(region =="Taiwan" | region == "Japan"){
    s2 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(region,"_s2",sep=""))))))
    s1 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(region,"_s1",sep=""))))))
  }else{
    s2 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(region,"_s2",sep=""))))))
    s1 <- s2
  }

  #Compute linear site term
  if(VS30 <= V1){
    Flin <- s1*log(VS30/V1) + s2*log(V1/Vref)
  }else if(VS30 <= V2){
    Flin <- s2*log(VS30/Vref)
  }else{
    Flin <- s2*log(V2/Vref)
  }

  # Nonlinear Site Term -----------------------------------------------------

  PGAr <- exp(GMM_at_760_Slab_v4(event.type,region,saturation.region,Rrup,M,hypocentral.depth,0,coefficients))
  f3 <- 0.05
  Vb <- 200
  Vref.Fnl <- 760

  if(period >= 3){
    Fnl = 0;
  }else{
    f2 <- coefficients.T$f4*(exp(coefficients.T$f5*(min(VS30,Vref.Fnl)-Vb))-exp(coefficients.T$f5*(Vref.Fnl-Vb)))
    Fnl = 0 + f2*log((PGAr+f3)/f3)
  }

  # Basin Term --------------------------------------------------------------

  #import coefficients

  if(Z2.5 == "default" | Z2.5 <= 0 | (region != "Japan" & region != "Cascadia")){
    Fb = 0
  }else{
    if(region == "Cascadia"){
      theta0 <- 3.94
      theta1 <- -0.42
      vmu <- 200
      vsig <- 0.2
      e1 <- coefficients.T$C_e1

      if(basin == 0){
        e3 <- coefficients.T$C_e3 + coefficients.T$del_None
        e2 <-  coefficients.T$C_e2 +  coefficients.T$del_None
      }else if(basin == 1){
        e3 <- coefficients.T$C_e3 + coefficients.T$del_Seattle
        e2 <-  coefficients.T$C_e2 +  coefficients.T$del_Seattle
      }else{
        e3 <- coefficients.T$C_e3
        e2 <-  coefficients.T$C_e2
      }

    }else if(region == "Japan"){
      theta0 <- 3.05
      theta1 <- -0.8
      vmu <- 500
      vsig <-0.33
      e3 <- coefficients.T$J_e3
      e2 <- coefficients.T$J_e2
      e1 <- coefficients.T$J_e1
    }

    Z2.5.pred <- 10^(theta0 + theta1*(1 + erf((log10(VS30) - log10(vmu))/(vsig*sqrt(2)))))
    delZ2.5 <- log(Z2.5) - log(Z2.5.pred)

    if(delZ2.5 <= (e1/e3)){
      Fb <- e1
    }else if(delZ2.5 >= (e2/e3)){
      Fb <- e2
    }else{
      Fb <- e3*delZ2.5
    }
  }

  # Add terms ---------------------------------------------------------------

  mu <-  c0 + Fm +Fp + Fd + Fnl + Fb+ Flin
  res <- list()
  res$mu <- mu
  res$Flin <- Flin
  res$Fnl <- Fnl
  res$Rock <- mu - ( Flin + Fnl )
  return(res)
}

PSHAB_Epi <- function(event.type,Reg,Ts,coefficients){

  #Grace Parker
  #July 23, 2020
  #Want to code up a simple epistemic uncertainty model that is a function of
  #Period, event type, and Reg

  #For interface, event.type = 0
  #For slab, event.type = 1
  # Ts ranges from 0.01-10s
  #Reg options = Global, Alaska, Aleutian, Cascadia, CAM_N, CAM_S, Japan_Pac, Japan_Phi, SA_N, SA_S, Taiwan


  # if(event.type == 0){
  # coefficients <- read.xlsx("Table_E4_Epistemic_Model_Coefficients.xlsx", sheet = "Interface")
  # }else if(event.type == 1){
  #   coefficients <- read.xlsx("Table_E4_Epistemic_Model_Coefficients.xlsx", sheet = "Slab")
  # }

  coefficients.R <- subset(coefficients, coefficients$Region == Reg)

  if(Ts < coefficients.R$T1 | Ts == 0){
    sigma_epi <- coefficients.R$SigEp1
  }else if(Ts > coefficients.R$T2){
    sigma_epi <- coefficients.R$SigEp2
  }else{
    sigma_epi <- coefficients.R$SigEp1 - (coefficients.R$SigEp1-coefficients.R$SigEp2)*(log(Ts/coefficients.R$T1)/log(coefficients.R$T2/coefficients.R$T1))
  }


  return(sigma_epi)
}






