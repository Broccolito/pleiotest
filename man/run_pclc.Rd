% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_mpat.R
\name{run_pclc}
\alias{run_pclc}
\title{Run PCLC Pleiotropy Test}
\usage{
run_pclc(pleio)
}
\arguments{
\item{pleio}{An object of class `pleio`, containing summary statistics for genetic variants.}
}
\value{
A `data.frame` with `PCLC`-based p-values for each variant.
The row names correspond to variant identifiers.
}
\description{
Applies the `PCLC` method from the `MPAT` package to assess pleiotropy.
}
\details{
The `PCLC` method (Principal Component-based Likelihood Combination) combines p-values across multiple correlated phenotypes using principal component analysis.
It accounts for the correlation structure among phenotypes and provides an aggregated p-value for each variant.
}
\examples{
result <- run_pclc(pleio)
head(result)

}
