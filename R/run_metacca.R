#' Run metaCCA Pleiotropy Test
#'
#' This function applies the `metaCCA` method to estimate pleiotropy significance across multiple phenotypes.
#'
#' @param pleio An object of class `pleio`, containing summary statistics and cohort details.
#'
#' @return A `data.frame` with the computed p-values using the `metaCCA` method.
#'
#' @examples
#' result <- run_metacca(pleio)
#' head(result)
#'
#' @export
run_metacca = function(pleio){

  n_phenotype = pleio@n_phenotype
  n_sample = sum(pleio@cohort_makeup_matrix)
  variant_names = rownames(pleio@summary_stats_matrix[[1]])

  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix
  zscore_sigma = cor(zscore_matrix)

  beta_matrix = as.data.frame(beta_matrix)
  names(beta_matrix) = gsub(pattern = "phenotype", replacement = "trait",
                            paste0(names(beta_matrix), "_b"))

  se_matrix = as.data.frame(se_matrix)
  names(se_matrix) = gsub(pattern = "phenotype", replacement = "trait",
                          paste0(names(se_matrix), "_se"))

  index = c()
  for(i in 1:n_phenotype){
    index = c(index, c(i, i+n_phenotype))
  }

  S_XY_study = cbind.data.frame(beta_matrix, se_matrix)[,index]
  S_XY_study[["allele_0"]] = "A"
  S_XY_study[["allele_1"]] = "T"

  S_XY_study = dplyr::select(S_XY_study, allele_0, allele_1, dplyr::everything()) |>
    dplyr::mutate(allele_0 = as.factor(allele_0)) |>
    dplyr::mutate(allele_1 = as.factor(allele_1))

  S_YY_study = metaCCA::estimateSyy(S_XY_study)

  metaCCA_res = metaCCA::metaCcaGp(nr_studies = 1,
                                   S_XY = list(S_XY_study),
                                   std_info = c(0),
                                   S_YY = list(S_YY_study),
                                   N = c(n_sample))

  pleio_p = data.frame(metacca_p = 10^(-metaCCA_res[,2]))
  rownames(pleio_p) = variant_names

  return(pleio_p)

}





