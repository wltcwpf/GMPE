% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bchydro_2019_subduction.R
\name{bchydro_2019}
\alias{bchydro_2019}
\title{The GMPE of BCHydro 2019 for subduction earthquakes (incoporated regional updates for Japan)}
\usage{
bchydro_2019(M, T = 1000, Rrup, Rhypo, F_faba, Vs30, F_event, Zh)
}
\arguments{
\item{M}{Moment magnitude, a numeric value}

\item{T}{Period (sec) (range should be between 0 s and 10 s). PSA at 0.01 s will be return for
any T <= 0.01 s.
Use 1000 (by default) for output the array of Sa with original BCHyrdo per-defined periods}

\item{Rrup}{Closest distance (km) to the ruptured plane. Use for Interface.}

\item{Rhypo}{Hypocentral distance(km). Use for Intraslab.}

\item{F_faba}{Flag for BackArc sites: 0 for Forearc or unknown sites; 1 for Backarc sites}

\item{Vs30}{The time averaged shear wave velocity on the top 30 m.}

\item{F_event}{Flag for interface event or intraslab events:
0 for Interface - use rupture distance; 1 for Intraslab - use hypocentral distance.}

\item{Zh}{The depth of hypocenter (km)}
}
\value{
A list of five elements is returned: Sa - median spectral acceleration prediction (in g);
sigma - totla logarithmic standard deviation (log);
phi - within event residuals logarithmic standard deviation (log);
tau - between event residuals logarithmic standard deviation (log);
period - the corresponding oscillator periods
}
\description{
This function calculates the ground motion median values and standard deviations
}
\examples{
bchydro_2019(M = 7.8, T = 1000, Rrup = 85, Rhypo = 100, F_faba = 0, Vs30 = 300, F_event = 0,
Zh = 25)

bchydro_2019(M = 7.8, T = c(0.010, 0.020, 0.030, 0.050, 0.075, 0.100, 0.150, 0.200, 0.250, 0.300,
0.400, 0.500, 0.600, 0.750, 1.000, 1.500, 2.000, 2.500, 3.000, 4.000, 5.000), Rrup = 85,
Rhypo = 100, F_faba = 1, Vs30 = 300, F_event = 1, Zh = 65)
}
\references{
Abrahamson, N., Gregor, N., and Addo, K. (2019). BC Hydro Ground Motion Prediction Equations
}
