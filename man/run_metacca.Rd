% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_metacca.R
\name{run_metacca}
\alias{run_metacca}
\title{Run metaCCA Pleiotropy Test}
\usage{
run_metacca(pleio)
}
\arguments{
\item{pleio}{An object of class `pleio`, containing summary statistics and cohort details.}
}
\value{
A `data.frame` with the computed p-values using the `metaCCA` method.
}
\description{
This function applies the `metaCCA` method to estimate pleiotropy significance across multiple phenotypes.
}
\examples{
result <- run_metacca(pleio)
head(result)

}
