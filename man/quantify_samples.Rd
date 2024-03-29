% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Function_Documentation.R
\name{quantify_samples}
\alias{quantify_samples}
\title{quantify_samples}
\usage{
quantify_samples(
  sample,
  set.env = F,
  max_calibration,
  use.area,
  quantified_unit,
  IS_deviation,
  use.MINDIST,
  gen.plots,
  Reference_deviation,
  minimum_r2,
  maximum_error_of_quantitation,
  minimum_nr_of_datapoints_calibration
)
}
\arguments{
\item{sample}{the name of the sample / file batch}

\item{set.env}{boolean if the environment of the R session should be setup from scratch}

\item{max_calibration}{to what extent (ratio, i.e. 1.2) are values above the maximum concentration still considered valid}

\item{use.area}{boolean, if TRUE uses peak area, otherwise peak height}

\item{quantified_unit}{the unit of the samples and final unit of quantified values}

\item{IS_deviation}{the value of accepted deviation for the IS in fractions (i.e. 0.2)}

\item{use.MINDIST}{boolean if MINDIST should be applied as peak selection criterion}

\item{gen.plots}{boolean if plots should be generated}

\item{Reference_deviation}{the value of accepted deviation for the IS in fractions (i.e. 0.2)}

\item{minimum_r2}{the lowest accepted coefficient of determination for valid calibrations}

\item{maximum_error_of_quantitation}{the maximal accepted deviation in \% of known and calculated concentrations}

\item{minimum_nr_of_datapoints_calibration}{the lowest number of data points a valid calibration should consist of}
}
\description{
Performs quantitative analysis (recovery/quantitation) of all samples
}
