#' Run AMATZ Pleiotropy Test
#'
#' Applies the `AMATZ` method from the `MTARclone` package to evaluate pleiotropy significance.
#' The method uses z-scores computed from effect size estimates and standard errors.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' The function extracts effect size (`beta`) and standard error (`SE`) matrices from the `pleio` object,
#' computes z-scores, and estimates their correlation structure.
#' The `AMATZ` method from `MTARclone` is then applied to test for pleiotropy.
#'
#' @return A `data.frame` with AMATZ-derived p-values for each variant.
#' The row names correspond to variant identifiers.
#'
#' @examples
#' result <- run_amatz(pleio)
#' head(result)
#'
#' @export
run_amatz = function(pleio){

  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  amatz_pvalue = function(mat, sig){
    MTARclone::amatz(mat, sig)$p.value[1]
  }

  pleio_p = apply(zscore_matrix, MARGIN = 1, amatz_pvalue, sig = zscore_sigma)

  pleio_p = data.frame(amatz_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)

}

#' Run CMATS Pleiotropy Test
#'
#' Applies the `CMATS` method from the `MTARclone` package to test for pleiotropy using a covariance structure.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' The function extracts the effect size (`beta`) and standard error (`SE`) matrices from the `pleio` object,
#' computes z-scores, and calculates their correlation structure.
#' The `CMATS` method is then used to estimate pleiotropy p-values.
#'
#' @return A `data.frame` containing CMATS-derived p-values for each variant.
#' The row names correspond to variant identifiers.
#'
#' @examples
#' result <- run_cmats(pleio)
#' head(result)
#'
#' @export
run_cmats = function(pleio){

  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  cmats_pvalue = function(mat, sig){
    MTARclone::cmats(mat, sig)$p.value[1]
  }

  pleio_p = apply(zscore_matrix, MARGIN = 1, cmats_pvalue, sig = zscore_sigma)

  pleio_p = data.frame(cmats_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)

}

#' Run EMATS Pleiotropy Test
#'
#' Applies the `EMATS` method from the `MTARclone` package to estimate pleiotropy significance.
#'
#' @param pleio An object of class `pleio`, containing summary statistics for genetic variants.
#'
#' @details
#' The function retrieves effect size (`beta`) and standard error (`SE`) matrices from the `pleio` object,
#' calculates z-scores, and estimates their correlation structure.
#' The `EMATS` method is then applied to test for pleiotropy.
#'
#' @return A `data.frame` with EMATS-derived p-values for each variant.
#' The row names correspond to variant identifiers.
#'
#' @examples
#' result <- run_emats(pleio)
#' head(result)
#'
#' @export
run_emats = function(pleio){

  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  emats_pvalue = function(mat, sig){
    MTARclone::emats(mat, sig)$p.value[1]
  }

  pleio_p = apply(zscore_matrix, MARGIN = 1, emats_pvalue, sig = zscore_sigma)

  pleio_p = data.frame(emats_p = pleio_p)
  rownames(pleio_p) = rownames(zscore_matrix)

  return(pleio_p)

}
