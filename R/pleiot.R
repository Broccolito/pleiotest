#' Pleiotropy Testing Class
#'
#' Defines a class for storing pleiotropy testing results and associated metadata.
#'
#' @slot pleio_object An object of class `pleio`, which contains simulation details.
#' @slot pleiotest_stats A `data.frame` containing pleiotropy test results for variants.
#' @slot execution_time A `numeric` value indicating the time taken for execution.
#'
#' @return An object of class `pleiot`.
setClass(
  Class = "pleiot",
  slots = list(
    pleio_object = "pleio",
    pleiotest_stats = "data.frame",
    execution_time = "numeric"
  )
)

#' Show Method for pleiot Class
#'
#' Prints a summary of the pleiot object.
#'
#' @param object An object of class `pleiot`.
#'
#' @return Prints a summary of the pleiot object, including phenotypes, SNPs, and applied methods.
#'
#' @examples
#' show(pleiotest)
setMethod(
  f = "show",
  signature = "pleiot",
  definition = function(object){
    cat("-------------------\n")
    cat("[pleiot object]\n")
    cat("Applying pleiotropy testing methods on the simulation of: \n")
    cat(paste0("- ", object@pleio_object@n_phenotype, " phenotypes\n"))
    cat(paste0("- ", (object@pleio_object@n_variant_pleiotropic +
                        sum(object@pleio_object@n_variant_nonpleiotropic) +
                        object@pleio_object@n_variant_null), " SNPs\n"))
    cat(paste0(" - ", object@pleio_object@n_variant_pleiotropic, " pleiotropic SNPs\n"))
    cat(paste0(" - ", sum(object@pleio_object@n_variant_nonpleiotropic), " non-pleiotropic SNPs\n"))
    cat(paste0(" - ", object@pleio_object@n_variant_null, " null SNPs\n"))
    cat(paste0("- Total of ", sum(object@pleio_object@defacto_sample_size_matrix), " participants\n"))
    cat(paste0(ncol(object@pleiotest_stats), " methods applied on a total of ", nrow(object@pleiotest_stats), " variants\n"))
    cat("-------------------\n")
  }
)
