% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func.R
\name{plot_lpdid}
\alias{plot_lpdid}
\title{Plot LP-DiD Event Study}
\usage{
plot_lpdid(
  reg,
  conf = 0.95,
  segments = TRUE,
  add = FALSE,
  xlab = NULL,
  ylab = NULL,
  ylim = NULL,
  main = "",
  x.shift = 0,
  pch = 19,
  cex = 1,
  col = "black",
  opacity = 1
)
}
\arguments{
\item{reg}{An object generated via the lpdid function.}

\item{conf}{The confidence level (1-alpha) desired for confidence intervals. Default is 0.95 (95% confidence internval)}

\item{segments}{A boolean (TRUE or FALSE) value for whether confidence intervals should be generated.  Default is TRUE.}

\item{add}{A boolean (TRUE or FALSE) value for whether the user wants to add to an existing figure or generate a new one.  The default is FALSE, meaning a new plot.}

\item{xlab}{The text belonging on the x axis.}

\item{ylab}{The text belonging on the y axis.}

\item{ylim}{A vector of length two that determines the minimum and maximum of the y-axis. Default scales to fit the figure being plotted.}

\item{main}{The text belonging as the figure's title.}

\item{x.shift}{A numeric value that will shift the event study estimates along the x-axis.}

\item{pch}{A numeric value corresponding to the point's shape.}

\item{cex}{A numeric value corresponding to the size of the point.}

\item{col}{The color the points (and confidence intervals) should be.}

\item{opacity}{A numeric value between 0 and 1 that corresponds to the opacity of the color.}
}
\value{
An event study plot.
}
\description{
Plotting LP-DiD event study parameter estimates and confidence intervals.
}
