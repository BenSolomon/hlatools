% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_drb345.R
\name{make_drb345_ratios}
\alias{make_drb345_ratios}
\title{Create DRB345 read ratios}
\usage{
make_drb345_ratios(drb)
}
\arguments{
\item{drb}{A data frame with columns: \cr
"sample" - unique sample name \cr
"DRB1", "DRB3", "DRB4", "DRB5" - raw read counts for each loci}
}
\value{
A data frame with columns: \cr
"sample" - unique sample name \cr
"DRB3", "DRB4", "DRB5" - ratio of read counts for each loci to those of DRB1
}
\description{
Create DRB345 read ratios
}
\examples{
data(drb345_data)
hlatools:::make_drb345_ratios(drb345_data)

}
