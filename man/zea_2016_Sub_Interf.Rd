% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zea_2016_gmpe.R
\name{zea_2016_Sub_Interf}
\alias{zea_2016_Sub_Interf}
\title{The GMPE of Zea 2016 for interface event}
\usage{
zea_2016_Sub_Interf(mag, T = 1000, fDepth, dist, DistV, SiteClass)
}
\arguments{
\item{mag}{Moment magnitude, a numeric value}

\item{T}{Period (sec) (range should be between 0.005 s and 5 s).
Use 1000 (by default) for output the array of Sa with original Zhao et al per-defined periods}

\item{fDepth}{Depth(km) to the top of ruptured plane (km)}

\item{dist}{The closest distance (km) to rupture plane when is available,
otherwise hypocentral distance.}

\item{DistV}{The distance of volcanic path (km)}

\item{SiteClass}{The site class: 1 for SC I, 2 for SC II, 3 for SC III 4 for SC IV sites. The relationship
between site class and Vs30 is given by Zhao et al., 2016 (https://doi.org/10.1785/0120150056).
SC I: Vs30 > 600 m/s; SC II: 300 < Vs30 =< 600 m/s; SC III: 200 < Vs30 =< 300 m/s; SC IV: Vs30 <= 200 m/s}
}
\value{
A list of six elements is returned: Sa - median spectral acceleration prediction (in g);
RockSa - median spectral acceleration prediction on Rock reference site (in g);
SiteAmp - Site amplication response spectral (subtraction in log(g)); sigma - totla standard deviation;
phiSS - single station standard deviation; period - the corresponding oscillator periods
}
\description{
This function calculates the ground motion median values and standard deviations
}
\examples{
zea_2016_Sub_Interf(mag = 7.5, T = 1000, fDepth = 0.0, dist = 40.0,
DistV = 0.0, SiteClass = 4)

zea_2016_Sub_Interf(mag = 7.5, T = c(0.025, 0.5), fDepth = 0.0, dist = 40.0,
DistV = 0.0, SiteClass = 4)
}
\references{
Zhao, et al. (2016). Ground‐Motion Prediction Equations for Subduction Interface Earthquakes
in Japan Using Site Class and Simple Geometric Attenuation Functions
Bulletin of the Seismological Society of America. 106(4): 1518-1534
}
