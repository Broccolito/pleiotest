#' Run TATES Pleiotropy Test
#'
#' Implements the TATES method for pleiotropy testing using correlation structure across multiple phenotypes.
#'
#' @param pleio An object of class `pleio`, containing summary statistics and correlation matrices.
#'
#' @return A `data.frame` containing TATES-based pleiotropy p-values for each variant.
#'
#' @examples
#' result <- run_tates(pleio)
#' head(result)
#'
#' @export
run_tates = function(pleio){
  variant_names = rownames(pleio@summary_stats_matrix[[1]])
  nvar = pleio@n_phenotype
  nsnp = nrow(pleio@eaf_matrix)

  beta_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,1]}) |>
    do.call(what = cbind)
  se_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,2]}) |>
    do.call(what = cbind)
  pvalue_matrix = purrr::map(pleio@summary_stats_matrix,function(x){x[,3]}) |>
    do.call(what = cbind)
  zscore_matrix = beta_matrix / se_matrix

  cormat = cor(zscore_matrix)
  pval = pvalue_matrix
  pval2 = pval

  r = as.matrix(cormat)

  betaa = c(-0.0007609278,-0.0022965148,0.6226249243,0.0148755138,
            0.1095155903,-0.0218930325,0.2178970393)

  synp = matrix(0,nsnp,1)

  for (isnp in 1:nsnp) {
    ps=pval2[isnp,]
    tmp=sort(ps,index.return=T)
    pj=tmp$x
    iorder=tmp$ix
    r2=matrix(0,nvar,nvar)
    r2=r[iorder,iorder]
    ro=diag(nvar)
    for (i1 in 1:nvar)  {
      for (i2 in 1:i1) {
        if (i1>i2) {
          er=r2[i1,i2]
          ro[i1,i2]=ro[i2,i1]= betaa[7]*er^6+betaa[6]*er^5+betaa[5]*er^4+betaa[4]*er^3+betaa[3]*er^2+betaa[2]*er+betaa[1]
        }}}
    alllam=eigen(ro[1:nvar,1:nvar])$values
    mepj=nvar
    for (i1 in 1:nvar) {
      mepj=mepj-(alllam[i1]>1)*(alllam[i1]-1) }
    mej=matrix(c(seq(1,nvar,1)),nvar,1,byrow=T)
    for (j in 1:nvar) {
      sellam=eigen(ro[1:j,1:j])$values
      id=j
      for (i1 in 1:id) {
        mej[j,1]=mej[j,1]-(sellam[i1]>1)*(sellam[i1]-1)
      }
    }
    pg=matrix(0,nvar,1)
    for (i in 1:nvar) {
      pg[i,1]=(mepj/mej[i,1])*pj[i]
    }
    pg=pg[iorder]
    synp[isnp,]=min(pg)
  }

  pleio_p = data.frame(tates_p = synp)
  rownames(pleio_p) = variant_names
  return(pleio_p)

}



