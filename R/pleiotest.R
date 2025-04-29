#' Run Comprehensive Pleiotropy Analysis
#'
#' This function runs multiple pleiotropy tests on a given `pleio` object, aggregating results
#' from various statistical methods. It allows for flexible selection of tests to be executed.
#'
#' @param pleio_object An object of class `pleio`, containing summary statistics for genetic variants.
#' @param run_amatz_test Logical; whether to run the AMATZ test. Default is `TRUE`.
#' @param run_cmats_test Logical; whether to run the CMATS test. Default is `TRUE`.
#' @param run_emats_test Logical; whether to run the EMATS test. Default is `TRUE`.
#' @param run_dsum_test Logical; whether to run the DSUM test. Default is `TRUE`.
#' @param run_metacca_test Logical; whether to run the metaCCA test. Default is `TRUE`.
#' @param run_minp_test Logical; whether to run the MinP test. Default is `TRUE`.
#' @param run_mixada_test Logical; whether to run the MixAda test. Default is `TRUE`.
#' @param run_mixfisher_davies_test Logical; whether to run the MixFisher test (Davies approximation). Default is `TRUE`.
#' @param run_mixfisher_liu_test Logical; whether to run the MixFisher test (Liu approximation). Default is `TRUE`.
#' @param run_mixfisher_liumod_test Logical; whether to run the MixFisher test (Liu-modified approximation). Default is `TRUE`.
#' @param run_mixsd_davies_test Logical; whether to run the MixSD test (Davies approximation). Default is `TRUE`.
#' @param run_mixsd_liu_test Logical; whether to run the MixSD test (Liu approximation). Default is `TRUE`.
#' @param run_mixsd_liumod_test Logical; whether to run the MixSD test (Liu-modified approximation). Default is `TRUE`.
#' @param run_mixtippett_davies_test Logical; whether to run the MixTippett test (Davies approximation). Default is `TRUE`.
#' @param run_mixtippett_liu_test Logical; whether to run the MixTippett test (Liu approximation). Default is `TRUE`.
#' @param run_mixtippett_liumod_test Logical; whether to run the MixTippett test (Liu-modified approximation). Default is `TRUE`.
#' @param run_mixvar_davies_test Logical; whether to run the MixVar test (Davies approximation). Default is `TRUE`.
#' @param run_mixvar_liu_test Logical; whether to run the MixVar test (Liu approximation). Default is `TRUE`.
#' @param run_mixvar_liumod_test Logical; whether to run the MixVar test (Liu-modified approximation). Default is `TRUE`.
#' @param run_pcaq_test Logical; whether to run the PCAQ test. Default is `FALSE`.
#' @param run_pcfisher_test Logical; whether to run the PCFisher test. Default is `TRUE`.
#' @param run_pclc_test Logical; whether to run the PCLC test. Default is `TRUE`.
#' @param run_pcminp_test Logical; whether to run the PCMinP test. Default is `TRUE`.
#' @param run_pco_test Logical; whether to run the PCO test. Default is `FALSE`.
#' @param run_sum_test Logical; whether to run the SUM test. Default is `TRUE`.
#' @param run_tates_test Logical; whether to run the TATES test. Default is `TRUE`.
#' @param run_vc_test Logical; whether to run the Variance Component (VC) test. Default is `TRUE`.
#' @param run_wald_test Logical; whether to run the Wald test. Default is `TRUE`.
#' @param run_wi_test Logical; whether to run the WI (Weighted Integration) test. Default is `TRUE`.
#'
#' @details
#' The `pleiotest` function systematically applies multiple statistical tests to evaluate pleiotropy
#' using a `pleio` object containing summary statistics from a GWAS-like study. The function supports
#' multiple pleiotropy detection methods, including `MixFisher`, `MixVar`, `MixTippett`, and `PCMinP`,
#' among others. Execution time for each test is recorded.
#'
#' The function returns an object of class `pleiot`, containing:
#' - A summary table of results from all selected tests
#' - Execution times for each test
#'
#' @return An object of class `pleiot` containing the aggregated test results and execution times.
#'
#' @examples
#' # Simulating a pleiotropic dataset
#' library(pleiosim)
#' sim_data <- pleiosim::simulate_pleiotropy(n_snps = 1000, n_traits = 5, n_samples = 500)
#'
#' # Running pleiotropy tests
#' results <- pleiotest(sim_data, run_mixfisher_davies_test = TRUE, run_minp_test = TRUE)
#'
#' # View results
#' head(results@pleiotest_stats)
#' results@execution_time
#'
#' @export
pleiotest = function(pleio_object,
                     run_amatz_test = TRUE,
                     run_cmats_test = TRUE,
                     run_emats_test = TRUE,
                     run_dsum_test = TRUE,
                     run_metacca_test = TRUE,
                     run_minp_test = TRUE,
                     run_mixada_test = TRUE,
                     run_mixfisher_davies_test = TRUE,
                     run_mixfisher_liu_test = TRUE,
                     run_mixfisher_liumod_test = TRUE,
                     run_mixsd_davies_test = TRUE,
                     run_mixsd_liu_test = TRUE,
                     run_mixsd_liumod_test = TRUE,
                     run_mixtippett_davies_test = TRUE,
                     run_mixtippett_liu_test = TRUE,
                     run_mixtippett_liumod_test = TRUE,
                     run_mixvar_davies_test = TRUE,
                     run_mixvar_liu_test = TRUE,
                     run_mixvar_liumod_test = TRUE,
                     run_pcaq_test = FALSE,
                     run_pcfisher_test = TRUE,
                     run_pclc_test = TRUE,
                     run_pcminp_test = TRUE,
                     run_pco_test = FALSE,
                     run_sum_test = TRUE,
                     run_tates_test = TRUE,
                     run_vc_test = TRUE,
                     run_wald_test = TRUE,
                     run_wi_test = TRUE){

  pleioanalyze_result = tibble::tibble(snp_name = rownames(pleio_object@summary_stats_matrix[[1]]))
  execution_time = numeric()

  if(run_amatz_test){
    message("Running amatz test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_amatz(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["amatz"] = time_taken
  }

  if(run_cmats_test){
    message("Running cmats test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_cmats(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["cmats"] = time_taken
  }

  if(run_emats_test){
    message("Running emats test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_emats(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["emats"] = time_taken
  }

  if(run_dsum_test){
    message("Running dsum test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_dsum(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["dsum"] = time_taken
  }

  if(run_metacca_test){
    message("Running metacca test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_metacca(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["metacca"] = time_taken
  }

  if(run_minp_test){
    message("Running minp test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_minp(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["minp"] = time_taken
  }

  if(run_mixada_test){
    message("Running mixada test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixada(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixada"] = time_taken
  }

  if(run_mixfisher_davies_test){
    message("Running mixfisher_davies test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixfisher_davies(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixfisher_davies"] = time_taken
  }

  if(run_mixfisher_liu_test){
    message("Running mixfisher_liu test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixfisher_liu(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixfisher_liu"] = time_taken
  }

  if(run_mixfisher_liumod_test){
    message("Running mixfisher_liumod test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixfisher_liumod(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixfisher_liumod"] = time_taken
  }

  if(run_mixsd_davies_test){
    message("Running mixsd_davies test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixsd_davies(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixsd_davies"] = time_taken
  }

  if(run_mixsd_liu_test){
    message("Running mixsd_liu test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixsd_liu(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixsd_liu"] = time_taken
  }

  if(run_mixsd_liumod_test){
    message("Running mixsd_liumod test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixsd_liumod(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixsd_liumod"] = time_taken
  }

  if(run_mixtippett_davies_test){
    message("Running mixtippett_davies test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixtippett_davies(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixtippett_davies"] = time_taken
  }

  if(run_mixtippett_liu_test){
    message("Running mixtippett_liu test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixtippett_liu(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixtippett_liu"] = time_taken
  }

  if(run_mixtippett_liumod_test){
    message("Running mixtippett_liumod test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixtippett_liumod(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixtippett_liumod"] = time_taken
  }

  if(run_mixvar_davies_test){
    message("Running mixvar_davies test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixvar_davies(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixvar_davies"] = time_taken
  }

  if(run_mixvar_liu_test){
    message("Running mixvar_liu test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixvar_liu(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixvar_liu"] = time_taken
  }

  if(run_mixvar_liumod_test){
    message("Running mixvar_liumod test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_mixvar_liumod(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["mixvar_liumod"] = time_taken
  }

  if(run_pcaq_test){
    message("Running pcaq test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_pcaq(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["pcaq"] = time_taken
  }

  if(run_pcfisher_test){
    message("Running pcfisher test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_pcfisher(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["pcfisher"] = time_taken
  }

  if(run_pclc_test){
    message("Running pclc test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_pclc(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["pclc"] = time_taken
  }

  if(run_pcminp_test){
    message("Running pcminp test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_pcminp(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["pcminp"] = time_taken
  }

  if(run_pco_test){
    message("Running pco test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_pco(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["pco"] = time_taken
  }

  if(run_sum_test){
    message("Running sum test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_sum(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["sum"] = time_taken
  }

  if(run_tates_test){
    message("Running tates test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_tates(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["tates"] = time_taken
  }

  if(run_vc_test){
    message("Running vc test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_vc(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["vc"] = time_taken
  }

  if(run_wald_test){
    message("Running wald test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_wald(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["wald"] = time_taken
  }

  if(run_wi_test){
    message("Running wi test...\n")
    start_time = Sys.time()
    pleioanalyze_result = cbind.data.frame(pleioanalyze_result,
                                           run_wi(pleio_object))
    end_time = Sys.time()
    time_taken = as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_time["wi"] = time_taken
  }

  rownames(pleioanalyze_result) = NULL

  pleiot = new("pleiot",
               pleio_object = pleio_object,
               pleiotest_stats = pleioanalyze_result,
               execution_time = execution_time
  )

  return(pleiot)
}
