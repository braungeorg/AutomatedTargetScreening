% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Function_Documentation.R
\name{calculate_system_dependent_shift_and_RT_range_woParallel}
\alias{calculate_system_dependent_shift_and_RT_range_woParallel}
\title{calculate_system_dependent_shift_and_RT_range_woParallel}
\usage{
calculate_system_dependent_shift_and_RT_range_woParallel(
  sample,
  set.env = F,
  minimum_search_window
)
}
\arguments{
\item{sample}{the name of the sample / file batch}

\item{set.env}{boolean if the environment of the R session should be setup from scratch}

\item{minimum_search_window}{the smallest size of the search window in minutes}
}
\description{
Calculates the system-dependent shift and defines the retention time range but without parallelization
}
