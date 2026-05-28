#' Ancestry-Informed STAAR-SPA procedure using omnibus test
#'
#' The \code{AI-STAAR_Binary_SPA} function takes in genotype, the object from fitting the null
#' model, and functional annotation data to analyze the association between a
#' imbalanced case-control phenotype and a variant-set by using an ancestry-informed STAAR-SPA procedure.
#' For each variant-set, the ancestry-informed STAAR-B p-value is a p-value from an omnibus test
#' that aggregated Burden(1,25) and Burden(1,1) together with p-values of each test weighted by each annotation
#' using Cauchy method, across an user-defined number of base tests. The p-values from each base test
#' are weighted by ancestry-specific ensemble weights estimated independently from the data.
#' @param genotype an n*p genotype matrix (dosage matrix) of the target sequence,
#' where n is the sample size and p is the number of genetic variants.
#' @param obj_nullmodel an object from fitting the null model, which is the
#' output from either \code{\link{fit_null_glm}} function for unrelated samples or
#' \code{\link{fit_null_glmmkin}} function for related samples. Note that \code{\link{fit_null_glmmkin}}
#' is a wrapper of the \code{\link{glmmkin}} function from the \code{\link{GMMAT}} package.
#' @param annotation_phred a data frame or matrix of functional annotation data
#' of dimension p*q (or a vector of a single annotation score with length p).
#' Continuous scores should be given in PHRED score scale, where the PHRED score
#' of j-th variant is defined to be -10*log10(rank(-score_j)/total) across the genome. (Binary)
#' categorical scores should be taking values 0 or 1, where 1 is functional and 0 is
#' non-functional. If not provided, STAAR will perform the
#' Burden(1,25) and Burden(1,1) tests (default = NULL).
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param rv_num_cutoff the cutoff of minimum number of variants of analyzing
#' a given variant-set (default = 2).
#' @param rv_num_cutoff_max the cutoff of maximum number of variants of analyzing
#' a given variant-set (default = 1e+09).
#' @param tol a positive number specifying tolerance, the difference threshold for parameter
#' estimates in saddlepoint apporximation algorithm below which iterations should be stopped (default = ".Machine$double.eps^0.25").
#' @param max_iter a positive integers pecifying the maximum number of iterations for applying the saddlepoint approximation algorithm (default = "1000").
#' @param SPA_p_filter logical: are only the variants with a normal approximation based p-value smaller than a pre-specified threshold use the SPA method to recalculate the p-value, only used for imbalanced case-control setting (default = FALSE).
#' @param p_filter_cutoff threshold for the p-value recalculation using the SPA method, only used for imbalanced case-control setting (default = 0.05).
#' @param find_weight logical: should the ancestry group-specific weights and weighting scenario-specific p-values for each base test be saved as output (default = FALSE).
#' @return A list with the following members:
#' @return \code{num_variant}: the number of variants with minor allele frequency > 0 and less than
#' \code{rare_maf_cutoff} in the given variant-set that are used for performing the
#' variant-set using STAAR.
#' @return \code{cMAC}: the cumulative minor allele count of variants with
#' minor allele frequency > 0 and less than \code{rare_maf_cutoff} in the given variant-set.
#' @return \code{RV_label}: the boolean vector indicating whether each variant in the given
#' variant-set has minor allele frequency > 0 and less than \code{rare_maf_cutoff}.
#' @return \code{results_STAAR_B}: the STAAR-B p-value that aggregated Burden(1,25) and Burden(1,1) together
#' with p-values of each test weighted by each annotation using Cauchy method.
#' @return \code{results_STAAR_B_1_25}: a vector of STAAR-B(1,25) p-values,
#' including Burden(1,25) p-value weighted by MAF, the Burden(1,25)
#' p-values weighted by each annotation, and a STAAR-B(1,25)
#' p-value by aggregating these p-values using Cauchy method.
#' @return \code{results_STAAR_B_1_1}: a vector of STAAR-B(1,1) p-values,
#' including Burden(1,1) p-value weighted by MAF, the Burden(1,1)
#' p-values weighted by each annotation, and a STAAR-B(1,1)
#' p-value by aggregating these p-values using Cauchy method.
#' @return \code{weight_all_1}: a matrix of ancestry-specific weights across
#' \code{B} base tests for scenario 1 (if \code{find_weight} = TRUE).
#' @return \code{weight_all_2}: a matrix of ancestry-specific weights across
#' \code{B} base tests for scenario 2 (if \code{find_weight} = TRUE).
#' @return \code{results_weight}: a list of p-values weighted by
#' MAF, annotations, with aggregated components across \code{B} base tests. The p-values
#' from each base test are combined across weighting scenarios 1 and 2 using Cauchy method
#' (if \code{find_weight} = TRUE).
#' @return \code{results_weight1}: a list of p-values weighted by
#' MAF, annotations, with aggregated components across \code{B} base tests. The p-values
#' from each base test correspond to weighting scenario 1 (if \code{find_weight} = TRUE).
#' @return \code{results_weight2}: a list of p-values weighted by
#' MAF, annotations, with aggregated components across \code{B} base tests. The p-values
#' from each base test correspond to weighting scenario 2 (if \code{find_weight} = TRUE).
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics}, \emph{52}(9), 969-983.
#' (\href{https://doi.org/10.1038/s41588-020-0676-4}{pub})
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting
#' noncoding rare-variant associations of large-scale whole-genome sequencing
#' studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @references Liu, Y., et al. (2019). Acat: A fast and powerful p value combination
#' method for rare-variant analysis in sequencing studies.
#' \emph{The American Journal of Human Genetics}, \emph{104}(3), 410-421.
#' (\href{https://doi.org/10.1016/j.ajhg.2019.01.002}{pub})
#' @references Li, Z., Li, X., et al. (2020). Dynamic scan procedure for
#' detecting rare-variant association regions in whole-genome sequencing studies.
#' \emph{The American Journal of Human Genetics}, \emph{104}(5), 802-814.
#' (\href{https://doi.org/10.1016/j.ajhg.2019.03.002}{pub})
#' @export

AI_STAAR_Binary_SPA <- function(genotype,obj_nullmodel,annotation_phred=NULL,
                                rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                rv_num_cutoff_max=1e9,
                                tol=.Machine$double.eps^0.25,max_iter=1000,
                                SPA_p_filter=FALSE,p_filter_cutoff=0.05,
                                find_weight=FALSE){
  a <- Sys.time()
  if(!inherits(genotype, "matrix") && !inherits(genotype, "Matrix")){
    stop("genotype is not a matrix!")
  }

  if(dim(genotype)[2] == 1){
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }

  annotation_phred <- as.data.frame(annotation_phred)
  if(dim(annotation_phred)[1] != 0 & dim(genotype)[2] != dim(annotation_phred)[1]){
    stop(paste0("Dimensions don't match for genotype and annotation!"))
  }

  if(inherits(genotype, "sparseMatrix")){
    genotype <- as.matrix(genotype)
  }
  genotype <- matrix_flip(genotype)
  genotype_ref <- genotype
  MAF <- genotype$MAF
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  Geno_rare <- genotype$Geno[,RV_label]

  rm(genotype)
  gc()
  annotation_phred <- annotation_phred[RV_label,,drop=FALSE]

  if(sum(RV_label) >= rv_num_cutoff_max){
    stop(paste0("Number of rare variant in the set is more than ",rv_num_cutoff_max,"!"))
  }

  B <- dim(obj_nullmodel$pop_weights_1_1)[2]
  weight_all_1 <- weight_all_2 <- vector("list", B)
  pvalues_1_tot <- pvalues_2_tot <- vector("list", B)

  n_pop <- length(unique(obj_nullmodel$pop.groups))
  pop <- obj_nullmodel$pop.groups
  indices <- list()
  a_p <- matrix(0, nrow = 1, ncol = n_pop)
  for(i in 1:n_pop){
    eth <-  unique(pop)[i]
    indices[[i]] <- which(pop %in% eth)
    a_p[,i] <- mean(apply(as.matrix(Geno_rare[indices[[i]],]), 2, function(x){
      min(mean(x)/2, 1-mean(x)/2)}))
  }
  a_p <- ifelse(a_p > 0, dbeta(a_p,1,25), a_p)

  if(sum(RV_label) >= rv_num_cutoff){
    G <- as(Geno_rare,"dgCMatrix")
    MAF <- MAF[RV_label]

    annotation_rank <- 1 - 10^(-annotation_phred/10)

    ## beta(1,25)
    w_1 <- dbeta(MAF,1,25)
    ## beta(1,1)
    w_2 <- dbeta(MAF,1,1)
    if(dim(annotation_phred)[2] == 0){
      ## Burden, SKAT, ACAT-V
      w_B <- w_S <- as.matrix(cbind(w_1,w_2))
      w_A <- as.matrix(cbind(w_1^2/dbeta(MAF,0.5,0.5)^2,w_2^2/dbeta(MAF,0.5,0.5)^2))
    }else{
      ## Burden
      w_B_1 <- annotation_rank*w_1
      w_B_1 <- cbind(w_1,w_B_1)
      w_B_2 <- annotation_rank*w_2
      w_B_2 <- cbind(w_2,w_B_2)
      w_B <- cbind(w_B_1,w_B_2)
      w_B <- as.matrix(w_B)

      ## SKAT
      w_S_1 <- sqrt(annotation_rank)*w_1
      w_S_1 <- cbind(w_1,w_S_1)
      w_S_2 <- sqrt(annotation_rank)*w_2
      w_S_2 <- cbind(w_2,w_S_2)
      w_S <- cbind(w_S_1,w_S_2)
      w_S <- as.matrix(w_S)

      ## ACAT-V
      w_A_1 <- annotation_rank*w_1^2/dbeta(MAF,0.5,0.5)^2
      w_A_1 <- cbind(w_1^2/dbeta(MAF,0.5,0.5)^2,w_A_1)
      w_A_2 <- annotation_rank*w_2^2/dbeta(MAF,0.5,0.5)^2
      w_A_2 <- cbind(w_2^2/dbeta(MAF,0.5,0.5)^2,w_A_2)
      w_A <- cbind(w_A_1,w_A_2)
      w_A <- as.matrix(w_A)
    }

    w_b_1 <- w_b_2 <- matrix(0, nrow = 1, ncol = n_pop)
    for(b in 1:B){
      print(b)
      #genotype <- genotype_ref
      #MAF <- genotype$MAF

      w_b_1 <- obj_nullmodel$pop_weights_1_1[,b]
      w_b_2 <- t(a_p%*%diag(obj_nullmodel$pop_weights_1_25[,b]))

      if(find_weight == T){
        weight_all_1[[b]] <- w_b_1
        weight_all_2[[b]] <- w_b_2
      }

      Geno_rare_1 <- Geno_rare_2 <- Geno_rare
      w_vec_1 <- w_vec_2 <- rep(NA,nrow(Geno_rare))

      for(i in 1:length(w_b_1)){
        eth <- unique(pop)[i]
        eth_wt_1 <- w_b_1[i]
        eth_wt_2 <- w_b_2[i]

        w_vec_1[indices[[i]]] <- eth_wt_1
        w_vec_2[indices[[i]]] <- eth_wt_2

        Geno_rare_1[indices[[i]],] <- eth_wt_1*Geno_rare[indices[[i]],]
        Geno_rare_2[indices[[i]],] <- eth_wt_2*Geno_rare[indices[[i]],]
      }

      #estimate the G_tilde matrix
      if(!is.null(obj_nullmodel$XSigma_i)){
        XW <- as.matrix(obj_nullmodel$XSigma_i)
        XXWX_inv <- as.matrix(obj_nullmodel$XXSigma_iX_inv)
      }else{
        XW <- obj_nullmodel$XW
        XXWX_inv <- obj_nullmodel$XXWX_inv
      }

      G1_tilde <-  as.matrix((Diagonal(x = w_vec_1) %*% Geno_rare_1) - (Diagonal(x = w_vec_1) %*% XXWX_inv %*% (XW %*% Geno_rare_1)))
      G2_tilde <-  as.matrix((Diagonal(x = w_vec_2) %*% Geno_rare_2) - (Diagonal(x = w_vec_2) %*% XXWX_inv %*% (XW %*% Geno_rare_2)))

      ## use SPA for all p-value calculation
      if(!SPA_p_filter){
        if(obj_nullmodel$relatedness){
          if(!obj_nullmodel$sparse_kins){

            residuals.phenotype <- obj_nullmodel$scaled.residuals
            muhat <- obj_nullmodel$fitted.values

            pvalues1 <- STAAR_B_Binary_SPA_wt(G=G1_tilde,residuals=residuals.phenotype,
                                              muhat=muhat,weights_B=w_B,
                                              tol=tol,max_iter=max_iter)
            pvalues2 <- STAAR_B_Binary_SPA_wt(G=G2_tilde,residuals=residuals.phenotype,
                                              muhat=muhat,weights_B=w_B,
                                              tol=tol,max_iter=max_iter)
          }else{
            residuals.phenotype <- obj_nullmodel$scaled.residuals
            muhat <- obj_nullmodel$fitted.values

            pvalues1 <- STAAR_B_Binary_SPA_wt(G=G1_tilde,residuals=residuals.phenotype,
                                              muhat=muhat,weights_B=w_B,
                                              tol=tol,max_iter=max_iter)
            pvalues2 <- STAAR_B_Binary_SPA_wt(G=G2_tilde,residuals=residuals.phenotype,
                                              muhat=muhat,weights_B=w_B,
                                              tol=tol,max_iter=max_iter)
          }
        }else{
          residuals.phenotype <- obj_nullmodel$y - obj_nullmodel$fitted.values
          muhat <- obj_nullmodel$fitted.values

          pvalues1 <- STAAR_B_Binary_SPA_wt(G=G1_tilde,residuals=residuals.phenotype,
                                            muhat=muhat,weights_B=w_B,
                                            tol=tol,max_iter=max_iter)
          pvalues2 <- STAAR_B_Binary_SPA_wt(G=G2_tilde,residuals=residuals.phenotype,
                                            muhat=muhat,weights_B=w_B,
                                            tol=tol,max_iter=max_iter)
        }
      }else{
        if(obj_nullmodel$relatedness){
          if(!obj_nullmodel$sparse_kins){

            residuals.phenotype <- obj_nullmodel$scaled.residuals
            muhat <- obj_nullmodel$fitted.values

            P <- obj_nullmodel$P

            pvalues1 <- STAAR_B_SPA_SMMAT_wt(G=G1_tilde,residuals=residuals.phenotype,
                                             muhat=muhat,weights_B=w_B,
                                             tol=tol,max_iter=max_iter,
                                             p_filter_cutoff=p_filter_cutoff,P=P)
            pvalues2 <- STAAR_B_SPA_SMMAT_wt(G=G2_tilde,residuals=residuals.phenotype,
                                             muhat=muhat,weights_B=w_B,
                                             tol=tol,max_iter=max_iter,
                                             p_filter_cutoff=p_filter_cutoff,P=P)
          }else{
            residuals.phenotype <- obj_nullmodel$scaled.residuals
            muhat <- obj_nullmodel$fitted.values

            Sigma_i <- obj_nullmodel$Sigma_i
            Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
            cov <- obj_nullmodel$cov

            pvalues1 <- STAAR_B_SPA_SMMAT_sparse_wt(G=G1_tilde,residuals=residuals.phenotype,
                                                    muhat=muhat,weights_B=w_B,
                                                    tol=tol,max_iter=max_iter,
                                                    p_filter_cutoff=p_filter_cutoff,
                                                    Sigma_i=Sigma_i,Sigma_iX=Sigma_iX,cov=cov)
            pvalues2 <- STAAR_B_SPA_SMMAT_sparse_wt(G=G2_tilde,residuals=residuals.phenotype,
                                                    muhat=muhat,weights_B=w_B,
                                                    tol=tol,max_iter=max_iter,
                                                    p_filter_cutoff=p_filter_cutoff,
                                                    Sigma_i=Sigma_i,Sigma_iX=Sigma_iX,cov=cov)
          }
        }else{
          residuals.phenotype <- obj_nullmodel$y - obj_nullmodel$fitted.values
          muhat <- obj_nullmodel$fitted.values

          X <- obj_nullmodel$X
          working <- obj_nullmodel$weights
          sigma <- sqrt(summary(obj_nullmodel)$dispersion)
          if(obj_nullmodel$family[1] == "binomial"){
            fam <- 1
          }else if(obj_nullmodel$family[1] == "gaussian"){
            fam <- 0
          }

          pvalues1 <- STAAR_B_SPA_wt(G=G1_tilde,residuals=residuals.phenotype,
                                     muhat=muhat,weights_B=w_B,
                                     tol=tol,max_iter=max_iter,
                                     p_filter_cutoff=p_filter_cutoff,
                                     X=X,working=working,sigma=sigma,fam=fam)
          pvalues2 <- STAAR_B_SPA_wt(G=G2_tilde,residuals=residuals.phenotype,
                                     muhat=muhat,weights_B=w_B,
                                     tol=tol,max_iter=max_iter,
                                     p_filter_cutoff=p_filter_cutoff,
                                     X=X,working=working,sigma=sigma,fam=fam)
        }
      }
      pvalues_1_tot[[b]] <- pvalues1
      pvalues_2_tot[[b]] <- pvalues2
    }

    pvalues_1_tot <- do.call(cbind, pvalues_1_tot)
    pvalues_2_tot <- do.call(cbind, pvalues_2_tot)

    ind_one_1 <- which(is.na(pvalues_1_tot), arr.ind = TRUE)
    ind_one_2 <- which(is.na(pvalues_2_tot), arr.ind = TRUE)
    pvalues_1_tot[ind_one_1] <- pvalues_2_tot[ind_one_2] <- 1

    pvalues_tot <- cbind(pvalues_1_tot,pvalues_2_tot)

    if(sum(pvalues_1_tot == 1) + sum(pvalues_2_tot == 1) == 0){
      pvalues_aggregate <- apply(pvalues_tot,1,function(x){CCT(x)})
    }else{
      pvalues_aggregate <- indices_NA <- c()
      for(j in 1:ncol(pvalues_tot)){
        dat <- as.numeric(formatC(pvalues_tot[,j], format="e", digits=50))
        if(sum(dat == 1) > 0){
          indices_NA <- c(indices_NA, j)
        }
      }

      for(i in 1:nrow(pvalues_tot)){
        dat <- as.numeric(formatC(pvalues_tot[i,-indices_NA], format="e", digits=50))
        pvalues_aggregate <- c(pvalues_aggregate, CCT(dat[dat != 1]))
      }
    }

    num_variant <- sum(RV_label) #dim(G)[2]
    cMAC <- sum(G)
    num_annotation <- dim(annotation_phred)[2]+1

    if(find_weight == TRUE){

      weight_all_1 <- do.call(cbind, weight_all_1)
      weight_all_2 <- do.call(cbind, weight_all_2)

      pvalues_aggregate_weight <- NULL
      results_weight <- results_weight1 <- results_weight2 <- NULL
      for(i in 1:B){
        pvalues_aggregate_weight <- cbind(pvalues_aggregate_weight,apply(as.matrix(pvalues_tot[,c(i,i+B)]),1,
                                                                         function(x){CCT(x)}))
      }
      for(i in 1:B){

        ## Combined p-values across 2 scenarios ##

        ## STAAR-B
        if(sum(is.na(pvalues_aggregate_weight[,i]))>0){
          ## all NAs
          if(sum(is.na(pvalues_aggregate_weight[,i]))==length(pvalues_aggregate_weight[,i])){
            results_STAAR_B <- 1
          }else{
            ## not all NAs
            pvalues_aggregate_sub <- pvalues_aggregate_weight[,i][!is.na(pvalues_aggregate_weight[,i])]
            if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
              ## not all ones
              results_STAAR_B <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
            }else{
              results_STAAR_B <- 1
            }
          }
        }else{
          if(sum(pvalues_aggregate_weight[,i][pvalues_aggregate_weight[,i]<1])>0){
            results_STAAR_B <- CCT(pvalues_aggregate_weight[,i][pvalues_aggregate_weight[,i]<1])
          }else{
            results_STAAR_B <- 1
          }
        }

        ## STAAR-B-1-25
        pvalues_aggregate_sub <- pvalues_aggregate_weight[,i][1:num_annotation]
        if(sum(is.na(pvalues_aggregate_sub))>0){
          if(sum(is.na(pvalues_aggregate_sub))==length(pvalues_aggregate_sub)){
            pvalues_aggregate_STAAR_B_1_25 <- 1
          }else{
            ## not all NAs
            pvalues_aggregate_sub <- pvalues_aggregate_sub[!is.na(pvalues_aggregate_sub)]
            if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
              ## not all ones
              pvalues_aggregate_STAAR_B_1_25 <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
            }else{
              pvalues_aggregate_STAAR_B_1_25 <- 1
            }
          }
        }else{
          if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
            pvalues_aggregate_STAAR_B_1_25 <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
          }else{
            pvalues_aggregate_STAAR_B_1_25 <- 1
          }
        }
        results_STAAR_B_1_25 <- c(pvalues_aggregate_weight[,i][1:num_annotation],pvalues_aggregate_STAAR_B_1_25)
        results_STAAR_B_1_25 <- data.frame(t(results_STAAR_B_1_25))

        ## STAAR-B-1-1
        pvalues_aggregate_sub <- pvalues_aggregate_weight[,i][(num_annotation+1):(2*num_annotation)]
        if(sum(is.na(pvalues_aggregate_sub))>0){
          if(sum(is.na(pvalues_aggregate_sub))==length(pvalues_aggregate_sub)){
            pvalues_aggregate_STAAR_B_1_1 <- 1
          }else{
            ## not all NAs
            pvalues_aggregate_sub <- pvalues_aggregate_sub[!is.na(pvalues_aggregate_sub)]
            if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
              ## not all ones
              pvalues_aggregate_STAAR_B_1_1 <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
            }else{
              pvalues_aggregate_STAAR_B_1_1 <- 1
            }
          }
        }else{
          if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
            pvalues_aggregate_STAAR_B_1_1 <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
          }else{
            pvalues_aggregate_STAAR_B_1_1 <- 1
          }
        }
        results_STAAR_B_1_1 <- c(pvalues_aggregate_weight[,i][(num_annotation+1):(2*num_annotation)],pvalues_aggregate_STAAR_B_1_1)
        results_STAAR_B_1_1 <- data.frame(t(results_STAAR_B_1_1))

        if(dim(annotation_phred)[2] == 0){
          colnames(results_STAAR_B_1_25) <- c("Burden(1,25)","STAAR-B(1,25)")
          colnames(results_STAAR_B_1_1) <- c("Burden(1,1)","STAAR-B(1,1)")
        }else{
          colnames(results_STAAR_B_1_25) <- c("Burden(1,25)",
                                              paste0("Burden(1,25)-",colnames(annotation_phred)),
                                              "STAAR-B(1,25)")
          colnames(results_STAAR_B_1_1) <- c("Burden(1,1)",
                                             paste0("Burden(1,1)-",colnames(annotation_phred)),
                                             "STAAR-B(1,1)")
        }
        results_weight <- cbind(results_weight, c(num_variant = num_variant,
                                                  cMAC = cMAC,
                                                  results_STAAR_B = results_STAAR_B,
                                                  results_STAAR_B_1_25 = results_STAAR_B_1_25,
                                                  results_STAAR_B_1_1 = results_STAAR_B_1_1))

        ## Scenario 1 p-values ##

        ## STAAR-B
        if(sum(is.na(pvalues_1_tot[,i]))>0){
          ## all NAs
          if(sum(is.na(pvalues_1_tot[,i]))==length(pvalues_1_tot[,i])){
            results_STAAR_B <- 1
          }else{
            ## not all NAs
            pvalues_aggregate_sub <- pvalues_1_tot[,i][!is.na(pvalues_1_tot[,i])]
            if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
              ## not all ones
              results_STAAR_B <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
            }else{
              results_STAAR_B <- 1
            }
          }
        }else{
          if(sum(pvalues_1_tot[,i][pvalues_1_tot[,i]<1])>0){
            results_STAAR_B <- CCT(pvalues_1_tot[,i][pvalues_1_tot[,i]<1])
          }else{
            results_STAAR_B <- 1
          }
        }

        ## STAAR-B-1-25
        pvalues_aggregate_sub <- pvalues_1_tot[,i][1:num_annotation]
        if(sum(is.na(pvalues_aggregate_sub))>0){
          if(sum(is.na(pvalues_aggregate_sub))==length(pvalues_aggregate_sub)){
            pvalues_aggregate_STAAR_B_1_25 <- 1
          }else{
            ## not all NAs
            pvalues_aggregate_sub <- pvalues_aggregate_sub[!is.na(pvalues_aggregate_sub)]
            if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
              ## not all ones
              pvalues_aggregate_STAAR_B_1_25 <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
            }else{
              pvalues_aggregate_STAAR_B_1_25 <- 1
            }
          }
        }else{
          if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
            pvalues_aggregate_STAAR_B_1_25 <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
          }else{
            pvalues_aggregate_STAAR_B_1_25 <- 1
          }
        }
        results_STAAR_B_1_25 <- c(pvalues_1_tot[,i][1:num_annotation],pvalues_aggregate_STAAR_B_1_25)
        results_STAAR_B_1_25 <- data.frame(t(results_STAAR_B_1_25))

        ## STAAR-B-1-1
        pvalues_aggregate_sub <- pvalues_1_tot[,i][(num_annotation+1):(2*num_annotation)]
        if(sum(is.na(pvalues_aggregate_sub))>0){
          if(sum(is.na(pvalues_aggregate_sub))==length(pvalues_aggregate_sub)){
            pvalues_aggregate_STAAR_B_1_1 <- 1
          }else{
            ## not all NAs
            pvalues_aggregate_sub <- pvalues_aggregate_sub[!is.na(pvalues_aggregate_sub)]
            if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
              ## not all ones
              pvalues_aggregate_STAAR_B_1_1 <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
            }else{
              pvalues_aggregate_STAAR_B_1_1 <- 1
            }
          }
        }else{
          if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
            pvalues_aggregate_STAAR_B_1_1 <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
          }else{
            pvalues_aggregate_STAAR_B_1_1 <- 1
          }
        }
        results_STAAR_B_1_1 <- c(pvalues_1_tot[,i][(num_annotation+1):(2*num_annotation)],pvalues_aggregate_STAAR_B_1_1)
        results_STAAR_B_1_1 <- data.frame(t(results_STAAR_B_1_1))

        if(dim(annotation_phred)[2] == 0){
          colnames(results_STAAR_B_1_25) <- c("Burden(1,25)","STAAR-B(1,25)")
          colnames(results_STAAR_B_1_1) <- c("Burden(1,1)","STAAR-B(1,1)")
        }else{
          colnames(results_STAAR_B_1_25) <- c("Burden(1,25)",
                                              paste0("Burden(1,25)-",colnames(annotation_phred)),
                                              "STAAR-B(1,25)")
          colnames(results_STAAR_B_1_1) <- c("Burden(1,1)",
                                             paste0("Burden(1,1)-",colnames(annotation_phred)),
                                             "STAAR-B(1,1)")
        }
        results_weight1 <- cbind(results_weight1, c(num_variant = num_variant,
                                                    cMAC = cMAC,
                                                    results_STAAR_B = results_STAAR_B,
                                                    results_STAAR_B_1_25 = results_STAAR_B_1_25,
                                                    results_STAAR_B_1_1 = results_STAAR_B_1_1))

        ## Scenario 2 p-values ##

        ## STAAR-B
        if(sum(is.na(pvalues_2_tot[,i]))>0){
          ## all NAs
          if(sum(is.na(pvalues_2_tot[,i]))==length(pvalues_2_tot[,i])){
            results_STAAR_B <- 1
          }else{
            ## not all NAs
            pvalues_aggregate_sub <- pvalues_2_tot[,i][!is.na(pvalues_2_tot[,i])]
            if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
              ## not all ones
              results_STAAR_B <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
            }else{
              results_STAAR_B <- 1
            }
          }
        }else{
          if(sum(pvalues_2_tot[,i][pvalues_2_tot[,i]<1])>0){
            results_STAAR_B <- CCT(pvalues_2_tot[,i][pvalues_2_tot[,i]<1])
          }else{
            results_STAAR_B <- 1
          }
        }

        ## STAAR-B-1-25
        pvalues_aggregate_sub <- pvalues_2_tot[,i][1:num_annotation]
        if(sum(is.na(pvalues_aggregate_sub))>0){
          if(sum(is.na(pvalues_aggregate_sub))==length(pvalues_aggregate_sub)){
            pvalues_aggregate_STAAR_B_1_25 <- 1
          }else{
            ## not all NAs
            pvalues_aggregate_sub <- pvalues_aggregate_sub[!is.na(pvalues_aggregate_sub)]
            if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
              ## not all ones
              pvalues_aggregate_STAAR_B_1_25 <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
            }else{
              pvalues_aggregate_STAAR_B_1_25 <- 1
            }
          }
        }else{
          if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
            pvalues_aggregate_STAAR_B_1_25 <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
          }else{
            pvalues_aggregate_STAAR_B_1_25 <- 1
          }
        }
        results_STAAR_B_1_25 <- c(pvalues_2_tot[,i][1:num_annotation],pvalues_aggregate_STAAR_B_1_25)
        results_STAAR_B_1_25 <- data.frame(t(results_STAAR_B_1_25))

        ## STAAR-B-1-1
        pvalues_aggregate_sub <- pvalues_2_tot[,i][(num_annotation+1):(2*num_annotation)]
        if(sum(is.na(pvalues_aggregate_sub))>0){
          if(sum(is.na(pvalues_aggregate_sub))==length(pvalues_aggregate_sub)){
            pvalues_aggregate_STAAR_B_1_1 <- 1
          }else{
            ## not all NAs
            pvalues_aggregate_sub <- pvalues_aggregate_sub[!is.na(pvalues_aggregate_sub)]
            if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
              ## not all ones
              pvalues_aggregate_STAAR_B_1_1 <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
            }else{
              pvalues_aggregate_STAAR_B_1_1 <- 1
            }
          }
        }else{
          if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
            pvalues_aggregate_STAAR_B_1_1 <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
          }else{
            pvalues_aggregate_STAAR_B_1_1 <- 1
          }
        }
        results_STAAR_B_1_1 <- c(pvalues_2_tot[,i][(num_annotation+1):(2*num_annotation)],pvalues_aggregate_STAAR_B_1_1)
        results_STAAR_B_1_1 <- data.frame(t(results_STAAR_B_1_1))

        if(dim(annotation_phred)[2] == 0){
          colnames(results_STAAR_B_1_25) <- c("Burden(1,25)","STAAR-B(1,25)")
          colnames(results_STAAR_B_1_1) <- c("Burden(1,1)","STAAR-B(1,1)")
        }else{
          colnames(results_STAAR_B_1_25) <- c("Burden(1,25)",
                                              paste0("Burden(1,25)-",colnames(annotation_phred)),
                                              "STAAR-B(1,25)")
          colnames(results_STAAR_B_1_1) <- c("Burden(1,1)",
                                             paste0("Burden(1,1)-",colnames(annotation_phred)),
                                             "STAAR-B(1,1)")
        }
        results_weight2 <- cbind(results_weight2, c(num_variant = num_variant,
                                                    cMAC = cMAC,
                                                    results_STAAR_B = results_STAAR_B,
                                                    results_STAAR_B_1_25 = results_STAAR_B_1_25,
                                                    results_STAAR_B_1_1 = results_STAAR_B_1_1))
      }
    }

    ## STAAR-B
    pvalues_aggregate_tot <- as.vector(as.numeric(formatC(pvalues_tot[,-indices_NA], format="e", digits=50)))
    if(sum(is.na(pvalues_aggregate_tot))>0){
      ## all NAs
      if(sum(is.na(pvalues_aggregate_tot))==length(pvalues_aggregate_tot)){
        results_STAAR_B <- 1
      }else{
        ## not all NAs
        pvalues_aggregate_tot_sub <- na.omit(as.vector(pvalues_tot))
        if(sum(pvalues_aggregate_tot_sub[pvalues_aggregate_tot_sub<1])>0){
          ## not all ones
          results_STAAR_B <- CCT(pvalues_aggregate_tot_sub[pvalues_aggregate_tot_sub<1])
        }else{
          results_STAAR_B <- 1
        }
      }
    }else{
      if(sum(pvalues_aggregate_tot[pvalues_aggregate_tot<1])>0){
        results_STAAR_B <- CCT(pvalues_aggregate_tot[pvalues_aggregate_tot<1])
      }else{
        results_STAAR_B <- 1
      }
    }

    ## STAAR-B-1-25
    pvalues_aggregate_sub <- pvalues_aggregate[1:num_annotation]
    if(sum(is.na(pvalues_aggregate_sub))>0){
      if(sum(is.na(pvalues_aggregate_sub))==length(pvalues_aggregate_sub)){
        pvalues_aggregate_STAAR_B_1_25 <- 1
      }else{
        ## not all NAs
        pvalues_aggregate_sub <- pvalues_aggregate_sub[!is.na(pvalues_aggregate_sub)]
        if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
          ## not all ones
          pvalues_aggregate_STAAR_B_1_25 <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
        }else{
          pvalues_aggregate_STAAR_B_1_25 <- 1
        }
      }
    }else{
      if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
        pvalues_aggregate_STAAR_B_1_25 <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
      }else{
        pvalues_aggregate_STAAR_B_1_25 <- 1
      }
    }
    results_STAAR_B_1_25 <- c(pvalues_aggregate[1:num_annotation],pvalues_aggregate_STAAR_B_1_25)
    results_STAAR_B_1_25 <- data.frame(t(results_STAAR_B_1_25))

    ## STAAR-B-1-1
    pvalues_aggregate_sub <- pvalues_aggregate[(num_annotation+1):(2*num_annotation)]
    if(sum(is.na(pvalues_aggregate_sub))>0){
      if(sum(is.na(pvalues_aggregate_sub))==length(pvalues_aggregate_sub)){
        pvalues_aggregate_STAAR_B_1_1 <- 1
      }else{
        ## not all NAs
        pvalues_aggregate_sub <- pvalues_aggregate_sub[!is.na(pvalues_aggregate_sub)]
        if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
          ## not all ones
          pvalues_aggregate_STAAR_B_1_1 <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
        }else{
          pvalues_aggregate_STAAR_B_1_1 <- 1
        }
      }
    }else{
      if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
        pvalues_aggregate_STAAR_B_1_1 <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
      }else{
        pvalues_aggregate_STAAR_B_1_1 <- 1
      }
    }
    results_STAAR_B_1_1 <- c(pvalues_aggregate[(num_annotation+1):(2*num_annotation)],pvalues_aggregate_STAAR_B_1_1)
    results_STAAR_B_1_1 <- data.frame(t(results_STAAR_B_1_1))

    if(dim(annotation_phred)[2] == 0){
      colnames(results_STAAR_B_1_25) <- c("Burden(1,25)","STAAR-B(1,25)")
      colnames(results_STAAR_B_1_1) <- c("Burden(1,1)","STAAR-B(1,1)")
    }else{
      colnames(results_STAAR_B_1_25) <- c("Burden(1,25)",
                                          paste0("Burden(1,25)-",colnames(annotation_phred)),
                                          "STAAR-B(1,25)")
      colnames(results_STAAR_B_1_1) <- c("Burden(1,1)",
                                         paste0("Burden(1,1)-",colnames(annotation_phred)),
                                         "STAAR-B(1,1)")
    }

    if(find_weight == T){
      return(list(num_variant = num_variant,
                  cMAC = cMAC,
                  RV_label = RV_label,
                  results_STAAR_B = results_STAAR_B,
                  results_STAAR_B_1_25 = results_STAAR_B_1_25,
                  results_STAAR_B_1_1 = results_STAAR_B_1_1,
                  weight_all_1 = weight_all_1,
                  weight_all_2 = weight_all_2,
                  results_weight = results_weight,
                  results_weight1 = results_weight1,
                  results_weight2 = results_weight2))
    }else{
      return(list(num_variant = num_variant,
                  cMAC = cMAC,
                  RV_label = RV_label,
                  results_STAAR_B = results_STAAR_B,
                  results_STAAR_B_1_25 = results_STAAR_B_1_25,
                  results_STAAR_B_1_1 = results_STAAR_B_1_1))
    }
  }else{
    stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
  }

}

