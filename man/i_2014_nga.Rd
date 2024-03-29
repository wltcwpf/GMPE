% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/i_2014_nga.R
\name{i_2014_nga}
\alias{i_2014_nga}
\title{The GMPE of CB 2014}
\usage{
i_2014_nga(M, Vs30, T = 1000, Rrup, Fault_Type)
}
\arguments{
\item{M}{Moment magnitude, a numeric value}

\item{Vs30}{Shear wave velocity averaged over top 30 m (in m/s).}

\item{T}{Period (sec); Use 1000 (by default) for output the array of Sa with original NGA West2 periods.
Note no PGV computation}

\item{Rrup}{Closest distance (km) to the ruptured plane}

\item{Fault_Type}{The indicator of fault type: 0 for strike slip & other non-reverse faulting;
1 for reverse}
}
\value{
A list of five elements is returned: med - median spectral acceleration prediction (in g);
sigma - logarithmic standard deviation of spectral acceleration prediction; phi - logarithmic
standard deviation of within event residuals; tau - logarithmic standard deviation of between
event residuals; period - the corresponding oscillator periods
}
\description{
This function calculates the ground motion median values and standard deviations
}
\examples{
i_2014_nga(M = 5.5, Vs30 = 550, T = 1000, Rrup = 90, Fault_Type = 1)
}
\references{
Idriss, I. M. (2014). An NGA-West2 Empirical Model for Estimating
the Horizontal Spectral Values Generated by Shallow Crustal Earthquakes.
Earthquake Spectra, 30(3), 1155-1177.
}
