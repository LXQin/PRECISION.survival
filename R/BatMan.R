#' Build survival prediction model based on training data using BatMan or ComBat along with normalization and calculate concordance index based on test data.
#'
#' @param data_train training dataset containing miRNA expression data and survival outcome. The first three columns need to be subject id, survival time, and censoring status (1=event, 0=censor), respectively. Expression data need to start at the fourth column.
#' @param data_test test dataset containing miRNA expression data and survival outcome. The first three columns need to be subject id, survival time, and censoring status, respectively. Expression data need to start at the fourth column.
#' @param nonzero_beta_position positions in the true \eqn{\beta} that have nonzero values. Only provide this argument if the true model is known.
#' @param b0value vector of nonzero values corresponding to \code{nonzero_beta_position}. Only provide this argument if the true model is known.
#' @param batch_id_train vector of batch id for all subjects in the training dataset indicating their batch membership. The length of the vector should be the size of the training dataset.
#' @param batch_id_test vector of batch id for all subjects in the test dataset indicating their batch membership. The length of the vector should be the size of the test dataset.
#' @param ps_mean_cutoff cutoff of mean expression level during pre-screening (default=8). Features with mean expression in the training dataset below this cutoff will be excluded during pre-screening.
#' @param ps_rho_cutoff cutoff of correlation coefficient during pre-screening (default=0.9).
#' @param BatMan_train 0 if no BatMan is used on training data; 1 otherwise.
#' @param BatMan_test 0 if no BatMan is used to calculate C-index from test data; 1 otherwise.
#' @param cmbt_train 0 if no ComBat is used on training data; 1 otherwise.
#' @param cmbt_test 0 if no ComBat is used to calculate C-index from test data; 1 otherwise.
#' @param norm_type normalization methods (0: no normalization; 1: median normalization; 2: quantile normalization; 3: Variance stabilizing transformation (VSN).
#' @param norm_train 0 if no normalization applied to training data; 1 otherwise. Ignored if \code{norm_type=0}.
#' @param norm_test 0 if no normalization applied to test data; 1 otherwise. Ignored if \code{norm_type=0}.
#' @param univ_cutoff_range range of cutoffs of p values for univariate analysis (default=\code{c(0,0.05)}).
#' @param nuniv_cutoff number of equal-spaced cutoff values in \code{univ_cutoff_range} (default=20).
#' @param lambda_max_glmnet maximum value of grid-search range of tuning parameter lambda in lasso and adaptive lasso (default=0.3).
#' @param nlambda number of equal-sapced lambda values in \code{lambda_max_glmnet} (default=30).
#' @param nfold number of folds for the cross-validation for tuning parameter selection in lasso and adaptive lasso (default=6).
#' @description
#' \code{BatMan} builds a survival prediction model based on training data using BatMan or ComBat along with three optional normalization methods 
#' (median normalization, quantile normalization, Variance stabilizing transformation). Four variable selection methods (oracle, univariate screening, 
#' lasso, adaptive lasso) are available to build a parsimonious model. Oracle method is only available when \code{nonzero_beta_position} and \code{b0value}
#' are provided. This is usually the case when the data is simulated so that the true model is known.
#' The survival prediction model is applied to the test data to calculate a concordance index (C-index) as a measure of prediction accuracy.
#' @return An object of class \code{BatMan.result} containing estimated regression coefficients and C-indices for each variable selection method.
#'   \item{bhat_oracle}{estimated regression coefficients from oracle method (only available for the true model is known)}
#'   \item{bhat_univ}{estimated regression coefficients from univariate method}
#'   \item{bhat_lasso}{estimated regression coefficients from lasso method}
#'   \item{bhat_alasso}{estimated regression coefficients from adaptive lasso method}
#'   \item{c_stats}{concordance indices (C-indices) based on test data from all variable selection methods}
#'   \item{lambdas}{selected tuning parameter lambda for Lasso and adaptive Lasso methods}
#'   \item{subset_size}{the sizes of data used for variable selection after pre-screening}
#'   \item{bhat_oracle_BatMan}{estimated regression coefficients from oracle method with BatMan}
#'   \item{bhat_univ_BatMan}{estimated regression coefficients from univariate method with BatMan}
#'   \item{bhat_lasso_BatMan}{estimated regression coefficients from lasso method with BatMan}
#'   \item{bhat_alasso_BatMan}{estimated regression coefficients from adaptive with BatMan}
#'   \item{c_stats_BatMan}{C-indices from all four variable selection methods based on test data with BatMan}
#'   \item{lambdas_BatMan}{selected tuning parameter lambda for Lasso and adaptive Lasso methods with BatMan}
#'   \item{bhat_oracle_ComBat}{estimated regression coefficients from oracle method with ComBat}
#'   \item{bhat_univ_ComBat}{estimated regression coefficients from univariate method with ComBat}
#'   \item{bhat_lasso_ComBat}{estimated regression coefficients from lasso method with ComBat}
#'   \item{bhat_alasso_ComBat}{estimated regression coefficients from adaptive lasso method with ComBat}
#'   \item{c_stats_ComBat}{C-indices from all four methods based on test data with ComBat}
#'   \item{lambdas_ComBat}{selected tuning parameter lambda for Lasso and adaptive Lasso methods with ComBat}
#'   \item{subset_size_ComBat}{the sizes of data used for variable selection after pre-screening with ComBat}
#' @import "MASS","survival","preprocessCore","glmnet","sva","vsn"
#' @examples
#' ## Build a survival prediction model using a simulated training dataset. BatMan and median normalization are used to correct for batch effects. Prediction accuracy
#' is measured by C-index calculated from applying the survival prediction model to a simulated test dataset. The test dataset is also subjected to BatMan and median
#' normalization.
#' data_sim=data_simulation(be.survival, he, nonzero_beta_position=c(385,866,1010,2218,2660,3026), 
#' b0value=c(1.04,1.40,1.45,1.62,3.13,1.76), batch_level="slides", t_sort_train=0, t_sort_test=0, 
#' t_rev=0, he_train=1, he_test=0, he_train_scale=1, he_test_scale=1, seed=123)
#' data_train=data_sim[[1]]
#' batch_id_train=data_sim[[2]]
#' data_test=data_sim[[3]]
#' batch_id_test=data_sim[[4]]
#' fit1=BatMan(data_train, data_test, nonzero_beta_position=c(385,866,1010,2218,2660,3026), b0value=c(1.04,1.40,1.45,1.62,3.13,1.76), 
#' batch_id_train=batch_id_train, batch_id_test=batch_id_test, BatMan_train=1, BatMan_test=1, cmbt_train=0, cmbt_test=0, norm_type=1, norm_train=1, norm_test=1)
#' @references Ni, A., Qin, L. (2021) "Performance Evaluation of Transcriptomics Data Normalization for Survival Risk Prediction". \emph{Briefings in Bioinformatics}, 22(6), bbab257.
#' @export

BatMan=function(data_train, data_test, nonzero_beta_position=NULL, b0value=NULL, batch_id_train, batch_id_test, ps_mean_cutoff=8, ps_rho_cutoff=0.9,
                BatMan_train, BatMan_test, cmbt_train, cmbt_test, norm_type, norm_train, norm_test, univ_cutoff_range=c(0,0.05), nuniv_cutoff=20, 
                lambda_max_glmnet=0.3, nlambda=30, nfold=6){
  
  names(data_train)[2:3]=c("t","delta")
  x_train_full=as.matrix(data_train[,4:ncol(data_train)])
  x_train_full_raw=x_train_full
  x_train_full_norm=x_train_full
  data_train_full_raw=data_train
  p=ncol(x_train_full)
  n=nrow(x_train_full)
  ncase=sum(data_train$delta)
  geneid=colnames(x_train_full)
  
  names(data_test)[2:3]=c("t","delta")
  x_test=as.matrix(data_test[,4:ncol(data_test)])
  x_test_raw=x_test
  x_test_norm=x_test
  data_test_raw=data_test
  n_test=nrow(x_test)

  ####### normalization ######
  
  if(norm_type==1){
    if(norm_train==1){
      x_train_full_rowmedian=apply(x_train_full_raw,1,median)
      x_train_full_rowmedian_diff=x_train_full_rowmedian-mean(x_train_full_rowmedian)
      x_train_full_norm=as.matrix(x_train_full_raw-x_train_full_rowmedian_diff%*%t(rep(1,ncol(x_train_full_raw))))
    }
    if(norm_test==1){
      x_test_rowmedian=apply(x_test_raw,1,median)
      x_test_rowmedian_diff=x_test_rowmedian-mean(x_train_full_rowmedian)
      x_test_norm=as.matrix(x_test_raw-x_test_rowmedian_diff%*%t(rep(1,ncol(x_test_raw))))
    }
  }
      
  if(norm_type==2){
    if(norm_train==1){
      x_train_full_norm=t(normalize.quantiles(t(x_train_full_raw)))
      x_quantile=sort(x_train_full_norm[1,])           
    }
    if(norm_test==1){
      x_test_norm=t(apply(x_test_raw, 1, function(u){x_quantile[rank(u)]}))
    }
  }
      
  if(norm_type==3){
    if(norm_train==1){
      data.vsn=vs.norm(train=t(x_train_full_raw))
      x_train_full_norm=t(data.vsn[[1]])
    }
    if(norm_test==1){
      x_test_norm=t(vs.norm(test=t(x_test_raw), ref.dis=data.vsn$ref.dis)[[2]])
    }
  }

  colnames(x_train_full_norm)=colnames(x_train_full_raw)
  data_train_full=cbind(data_train_full_raw[,1:3],x_train_full_norm) 
  colnames(x_test_norm)=colnames(x_test_raw)
  data_test=cbind(data_test_raw[,1:3],x_test_norm)

  
  #### pre-screen

  means_train=colMeans(x_train_full_norm)
  x_train_sub=x_train_full_norm[,means_train>ps_mean_cutoff]
  xcorr_train=cor(x_train_sub)
  exclude_train=c()
  for(i in 1:ncol(x_train_sub)){
    for(j in i:ncol(x_train_sub)){
      if(xcorr_train[i,j]>ps_rho_cutoff & xcorr_train[i,j]!=1){
        exclude_train_1=ifelse(means_train[geneid==rownames(xcorr_train)[i]]>=means_train[geneid==colnames(xcorr_train)[j]],
                       colnames(xcorr_train)[j], rownames(xcorr_train)[i])
        exclude_train=c(exclude_train, exclude_train_1)
      }
    }
  }
  if(length(unique(exclude_train))==0){
    x_train=x_train_sub
  }else{
    x_train=x_train_sub[,!(colnames(x_train_sub) %in% unique(exclude_train))]
  }
  p_sub=ncol(x_train)
  subsize=p_sub
  data_train=cbind(data_train_full[,1:3],x_train)


  ################ analysis without BatMan or ComBat ##############

  if(BatMan_train==0 & cmbt_train==0){
    # oracle analysis without BatMan or ComBat
    if(!(is.null(nonzero_beta_position))){
      b0=rep(0,p)
      b0[nonzero_beta_position]=b0value
      x0_train=x_train_full_norm[,b0!=0]
      coxfit_o=tryCatch(coxph(Surv(data_train$t,data_train$delta)~as.matrix(x0_train)), error=function(e) e, warning=function(w) w)
      if(is(coxfit_o, "warning") | is(coxfit_o, "error")){
        coxfit_o=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x0_train), theta=1))
      }
      b_o=rep(0,p)
      b_o[b0!=0]=summary(coxfit_o)$coefficient[,1]
    
      x0_test=x_test_norm[,b0!=0]
      coxfit_o_test=coxph(Surv(data_test$t,data_test$delta)~as.matrix(x0_test), init=summary(coxfit_o)$coefficient[,1], control=coxph.control(iter.max=0))
      c_o_test=summary(coxfit_o_test)$concordance[1]
    }
    
    # univariate analysis without BatMan or ComBat
    ps_u=lkhd_u=rep(0,p_sub)
    for(i in 1:p_sub){
      coxfit_ui=coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,i]))
      ps_u[i]=summary(coxfit_ui)$coefficient[1,5]
      lkhd_u[i]=coxfit_ui$loglik[2]
    }
    cutoff_grid=seq(univ_cutoff_range[1], univ_cutoff_range[2], length.out=nuniv_cutoff)
    aics=rep(NA, nuniv_cutoff)
    k=1
    for(cut in cutoff_grid){
      selected=1*(ps_u<=cut)
      if(sum(selected)!=0){
        coxfit_u=tryCatch(coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,selected==1])), error=function(e) e, warning=function(w) w)
        if(!(is(coxfit_u, "warning") | is(coxfit_u, "error"))){
          if(max(abs(summary(coxfit_u)$coefficient[,1]))<15){
            aics[k]=extractAIC(coxfit_u)[2]
          }
        }
      }
      k=k+1
    }
    if(sum(is.na(aics))==nuniv_cutoff){
      if(sum(selected)==0){
        sel_genes=colnames(x_train)[ps_u==min(ps_u)]
        coxfit_u=coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,sel_genes]))
      }else{
        sel_genes=colnames(x_train)[ps_u<=univ_cutoff_range[2]]
        x_selu=as.matrix(x_train[,sel_genes])
        pselu=ps_u[ps_u<=univ_cutoff_range[2]]
        coxfit_u=tryCatch(coxph(Surv(data_train$t,data_train$delta)~ridge(x_selu, theta=1)), error=function(e) e, warning=function(w) w)
        while(is(coxfit_u, "warning") | is(coxfit_u, "error")){
          sel_genes=sel_genes[-which(pselu==max(pselu))]
          x_selu=as.matrix(x_train[,sel_genes])
          pselu=pselu[-which(pselu==max(pselu))]
          coxfit_u=tryCatch(coxph(Surv(data_train$t,data_train$delta)~ridge(x_selu, theta=1)), error=function(e) e, warning=function(w) w)
        }
      }
    }else{
      cut_sel=cutoff_grid[which(aics==min(aics, na.rm=T))[1]]
      sel_genes=colnames(x_train)[ps_u<=cut_sel]
      coxfit_u=coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,ps_u<=cut_sel]))
    }
    b_u=rep(0,p)
    b_u[geneid %in% sel_genes]=summary(coxfit_u)$coefficient[,1]

    coxfit_u_test=coxph(Surv(data_test$t,data_test$delta)~as.matrix(x_test_norm[,geneid %in% sel_genes]), init=summary(coxfit_u)$coefficient[,1], control=coxph.control(iter.max=0))
    c_u_test=summary(coxfit_u_test)$concordance[1]      
    

  
    lkhd_p_sort=sort(lkhd_u, decreasing=TRUE)
    lkhd_p_thres=lkhd_p_sort[round(min(ncase/4,p_sub/4))]
    sel_p_genes=colnames(x_train)[lkhd_u>=lkhd_p_thres]
    x_p_train=as.matrix(x_train[,lkhd_u>=lkhd_p_thres])
    inifit_p=tryCatch(coxph(Surv(data_train$t, data_train$delta)~as.matrix(x_train[,lkhd_u>=lkhd_p_thres])), error=function(e) e, warning=function(w) w)
    if(is(inifit_p, "warning") | is(inifit_p, "error") | max(abs(coef(inifit_p)))>10){
      inifit_p=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x_train[,lkhd_u>=lkhd_p_thres]), theta=1))
    }
    x_p_test=as.matrix(x_test_norm[,geneid %in% sel_p_genes])
  
    # Lasso-penalized analysis without BatMan or ComBat
    alpha_grid=seq(1,1,0.1)
    cvs_alpha=lambdas_alpha=c()
    for(alpha in alpha_grid){
      cv_l=cv.glmnet(as.matrix(x_train[,lkhd_u>=lkhd_p_thres]), Surv(data_train$t, data_train$delta), nfolds=nfold, family="cox", alpha=alpha, standardize=F)
      lambda_min=cv_l$lambda.min
      cvm_min=cv_l$cvm[cv_l$lambda==lambda_min]
      cvs_alpha=c(cvs_alpha, cvm_min)
      lambdas_alpha=c(lambdas_alpha, lambda_min)
    }
    lambda_sel=lambdas_alpha[which(cvs_alpha==min(cvs_alpha,na.rm=T))][1]
    alpha_sel=alpha_grid[which(cvs_alpha==min(cvs_alpha,na.rm=T))][1]
    glmnet_l=glmnet(as.matrix(x_train[,lkhd_u>=lkhd_p_thres]), Surv(data_train$t, data_train$delta), family="cox", alpha=alpha_sel, lambda=lambda_sel, standardize=F)
    b_l0temp=as.vector(glmnet_l$beta)
    b_l=rep(0,p)
    b_l[geneid %in% sel_p_genes]=b_l0temp
    l_l=lambda_sel
  
    coxfit_l_test=coxph(Surv(data_test$t,data_test$delta)~x_p_test, init=b_l0temp, control=coxph.control(iter.max=0))
    c_l_test=summary(coxfit_l_test)$concordance[1]
  
  
    # Adaptive Lasso-penalized analysis without BatMan or ComBat
    w_a=1/abs(coef(inifit_p))
  
    cvs_alpha=lambdas_alpha=c()
    for(alpha in alpha_grid){
      cv_a=cv.glmnet(as.matrix(x_train[,lkhd_u>=lkhd_p_thres]), Surv(data_train$t, data_train$delta), nfolds=nfold, family="cox", alpha=alpha, standardize=F, penalty.factor=w_a)
      lambda_min=cv_a$lambda.min
      cvm_min=cv_a$cvm[cv_a$lambda==lambda_min]
      cvs_alpha=c(cvs_alpha, cvm_min)
      lambdas_alpha=c(lambdas_alpha, lambda_min)
    }
    lambda_sel=lambdas_alpha[which(cvs_alpha==min(cvs_alpha,na.rm=T))][1]
    alpha_sel=alpha_grid[which(cvs_alpha==min(cvs_alpha,na.rm=T))][1]
    glmnet_a=glmnet(as.matrix(x_train[,lkhd_u>=lkhd_p_thres]), Surv(data_train$t, data_train$delta), family="cox", alpha=alpha_sel, lambda=lambda_sel, standardize=F, penalty.factor=w_a)
    b_a0temp=as.vector(glmnet_a$beta)
    b_a=rep(0,p)
    b_a[geneid %in% sel_p_genes]=b_a0temp
    l_a=lambda_sel
  
    coxfit_a_test=coxph(Surv(data_test$t,data_test$delta)~x_p_test, init=b_a0temp, control=coxph.control(iter.max=0))
    c_a_test=summary(coxfit_a_test)$concordance[1]    
  }



  ################ BatMan analysis ##############

  if(BatMan_train==1){

  # oracle analysis with BatMan
  if(!(is.null(nonzero_beta_position))){
    b0=rep(0,p)
    b0[nonzero_beta_position]=b0value
    x0_train=x_train_full_norm[,b0!=0]
    coxfit_os=tryCatch(coxph(Surv(data_train$t,data_train$delta)~as.matrix(x0_train)+strata(batch_id_train)), error=function(e) e, warning=function(w) w)
    if(is(coxfit_os, "warning") | is(coxfit_os, "error")){
      coxfit_os=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x0_train), theta=1)+strata(batch_id_train))
    }
    b_os=rep(0,p)
    b_os[b0!=0]=summary(coxfit_os)$coefficient[,1]
  
    x0_test=x_test_norm[,b0!=0]
    if(BatMan_test==0){
      coxfit_os_test=coxph(Surv(data_test$t,data_test$delta)~as.matrix(x0_test), init=summary(coxfit_os)$coefficient[,1], control=coxph.control(iter.max=0))
    }
    if(BatMan_test==1){
      coxfit_os_test=coxph(Surv(data_test$t,data_test$delta)~as.matrix(x0_test)+strata(batch_id_test), init=summary(coxfit_os)$coefficient[,1], control=coxph.control(iter.max=0))
    }
    c_os_test=summary(coxfit_os_test)$concordance[1]
  }

    # univariate analysis with BatMan
    ps_us=lkhd_us=rep(0,p_sub)
    for(i in 1:p_sub){
      coxfit_usi=coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,i]) + strata(batch_id_train))
      ps_us[i]=summary(coxfit_usi)$coefficient[1,5]
      lkhd_us[i]=coxfit_usi$loglik[2]
    }
    cutoff_grid=seq(univ_cutoff_range[1], univ_cutoff_range[2], length.out=nuniv_cutoff)
    aics=rep(NA, nuniv_cutoff)
    k=1
    for(cut in cutoff_grid){
      selected=1*(ps_us<=cut)
      if(sum(selected)!=0){
        coxfit_us=tryCatch(coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,selected==1])+strata(batch_id_train)), error=function(e) e, warning=function(w) w)
        if(!(is(coxfit_us, "warning") | is(coxfit_us, "error"))){
          if(max(abs(summary(coxfit_us)$coefficient[,1]))<15){
            aics[k]=extractAIC(coxfit_us)[2]
          }
        }
      }
      k=k+1
    }
    if(sum(is.na(aics))==nuniv_cutoff){
      sel_genes_strat=colnames(x_train)[ps_us==min(ps_us)]
      coxfit_us=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x_train[,ps_us==min(ps_us)]), theta=1)+strata(batch_id_train))
    }else{
      cut_sel=cutoff_grid[which(aics==min(aics, na.rm=T))[1]]
      sel_genes_strat=colnames(x_train)[ps_us<=cut_sel]
      coxfit_us=coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,ps_us<=cut_sel])+strata(batch_id_train))
    }
    b_us=rep(0,p)
    b_us[geneid %in% sel_genes_strat]=summary(coxfit_us)$coefficient[,1]

    if(BatMan_test==0){
      coxfit_us_test=coxph(Surv(data_test$t,data_test$delta)~as.matrix(x_test_norm[,geneid %in% sel_genes_strat]), init=summary(coxfit_us)$coefficient[,1], control=coxph.control(iter.max=0))
    }
    if(BatMan_test==1){
      coxfit_us_test=coxph(Surv(data_test$t,data_test$delta)~as.matrix(x_test_norm[,geneid %in% sel_genes_strat]) + strata(batch_id_test), init=summary(coxfit_us)$coefficient[,1], control=coxph.control(iter.max=0))
    }
    c_us_test=summary(coxfit_us_test)$concordance[1]


    lkhd_ps_sort=sort(lkhd_us, decreasing=TRUE)
    lkhd_ps_thres=lkhd_ps_sort[round(min(ncase/4,p_sub/4))]
    sel_p_genes_strat=colnames(x_train)[lkhd_us>=lkhd_ps_thres]
    x_ps_train=as.matrix(x_train[,lkhd_us>=lkhd_ps_thres])
    inifit_ps=tryCatch(coxph(Surv(data_train$t, data_train$delta)~as.matrix(x_train[,lkhd_us>=lkhd_ps_thres])+strata(batch_id_train)),
                       error=function(e) e, warning=function(w) w)
    if(is(inifit_ps, "warning") | is(inifit_ps, "error")){
      inifit_ps=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x_train[,lkhd_us>=lkhd_ps_thres]), theta=1)+strata(batch_id_train))
    }else{
      if(max(abs(coef(inifit_ps)))>10){
        inifit_ps=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x_train[,lkhd_us>=lkhd_ps_thres]), theta=1)+strata(batch_id_train))
      }
    }

    x_nps_test=as.matrix(x_test_norm[,geneid %in% sel_p_genes_strat])

    # Lasso-penalized analysis with BatMan

    pfit_ls=glmnet_strat_cv(data=cbind(data_train[,1:3], batch_id_train, x_train[,lkhd_us>=lkhd_ps_thres]), lambda_max=lambda_max_glmnet, nlambda=nlambda, nfold=nfold, penalty_wt=rep(1,n), threshold=10^(-5))
    b_ls=rep(0,p)
    b_ls0temp=pfit_ls[[1]]
    b_ls[geneid %in% sel_p_genes_strat]=b_ls0temp
    l_ls=pfit_ls[[2]]

    if(BatMan_test==0){
      coxfit_ls_test=coxph(Surv(data_test$t,data_test$delta)~x_nps_test, init=b_ls0temp, control=coxph.control(iter.max=0))
    }
    if(BatMan_test==1){
      coxfit_ls_test=coxph(Surv(data_test$t,data_test$delta)~x_nps_test + strata(batch_id_test), init=b_ls0temp, control=coxph.control(iter.max=0))
    }
    c_ls_test=summary(coxfit_ls_test)$concordance[1]

    # Adaptive Lasso-penalized analysis with BatMan
    w_a=1/abs(coef(inifit_ps))

    pfit_as=glmnet_strat_cv(data=cbind(data_train[,1:3], batch_id_train, x_train[,lkhd_us>=lkhd_ps_thres]), lambda_max=lambda_max_glmnet, nlambda=nlambda, nfold=nfold, penalty_wt=w_a, threshold=10^(-5))
    b_as=rep(0,p)
    b_as0temp=pfit_as[[1]]
    b_as[geneid %in% sel_p_genes_strat]=b_as0temp
    l_as=pfit_as[[2]]

    if(BatMan_test==0){
      coxfit_as_test=coxph(Surv(data_test$t,data_test$delta)~x_nps_test, init=b_as0temp, control=coxph.control(iter.max=0))
    }
    if(BatMan_test==1){
      coxfit_as_test=coxph(Surv(data_test$t,data_test$delta)~x_nps_test + strata(batch_id_test), init=b_as0temp, control=coxph.control(iter.max=0))
    }
    c_as_test=summary(coxfit_as_test)$concordance[1]
  }


  ################ analysis with ComBat ##############
  
  if(cmbt_train==1){
    x_train_full=t(ComBat(dat=t(x_train_full_norm), batch=batch_id_train, mod=NULL, par.prior=F, mean.only=T))  
    data_train_full=cbind(data_train_full_raw[,1:3],x_train_full)

    if(cmbt_test==1){
      x_test=t(ComBat(dat=t(x_test_norm), batch=batch_id_test, mod=NULL, par.prior=F, mean.only=T))
      data_test=cbind(data_test_raw[,1:3],x_test)
    }else{
      x_test=x_test_norm
      data_test=data_test
    }
    
    ## pre-screen
    
    means_train=colMeans(x_train_full)
    x_train_sub=x_train_full[,means_train>ps_mean_cutoff]
    xcorr_train=cor(x_train_sub)
    exclude_train=c()
    for(i in 1:ncol(x_train_sub)){
      for(j in i:ncol(x_train_sub)){
        if(xcorr_train[i,j]>ps_rho_cutoff & xcorr_train[i,j]!=1){
          exclude_train_1=ifelse(means_train[geneid==rownames(xcorr_train)[i]]>=means_train[geneid==colnames(xcorr_train)[j]],
                                 colnames(xcorr_train)[j], rownames(xcorr_train)[i])
          exclude_train=c(exclude_train, exclude_train_1)
        }
      }
    }
    if(length(unique(exclude_train))==0){
      x_train=x_train_sub
    }else{
      x_train=x_train_sub[,!(colnames(x_train_sub) %in% unique(exclude_train))]
    }
    p_sub=ncol(x_train)
    subsize_c=p_sub
    data_train=cbind(data_train_full[,1:3],x_train)
    
    
    # oracle analysis with ComBat
    if(!(is.null(nonzero_beta_position))){
      b0=rep(0,p)
      b0[nonzero_beta_position]=b0value
      x0_train=x_train_full[,b0!=0]
      coxfit_oc=tryCatch(coxph(Surv(data_train$t,data_train$delta)~as.matrix(x0_train)),error=function(e) e, warning=function(w) w)
      if(is(coxfit_oc, "warning") | is(coxfit_oc, "error")){
        coxfit_oc=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x0_train), theta=1))
      }
      b_oc=rep(0,p)
      b_oc[b0!=0]=summary(coxfit_oc)$coefficient[,1]
  
      x0_test=x_test[,b0!=0]
      coxfit_oc_test=coxph(Surv(data_test$t,data_test$delta)~as.matrix(x0_test), init=summary(coxfit_oc)$coefficient[,1], control=coxph.control(iter.max=0))
      c_oc_test=summary(coxfit_oc_test)$concordance[1]
    }
    
    
    # univariate analysis with ComBat
    ps_uc=lkhd_uc=rep(0,p_sub)
    for(i in 1:p_sub){
      coxfit_uci=coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,i]))
      ps_uc[i]=summary(coxfit_uci)$coefficient[1,5]
      lkhd_uc[i]=coxfit_uci$loglik[2]
    }
    cutoff_grid=seq(univ_cutoff_range[1], univ_cutoff_range[2], length.out=nuniv_cutoff)
    aics=rep(NA, nuniv_cutoff)
    k=1
    for(cut in cutoff_grid){
      selected=1*(ps_uc<=cut)
      if(sum(selected)!=0){
        coxfit_uc=tryCatch(coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,selected==1])), error=function(e) e, warning=function(w) w)
        if(!(is(coxfit_uc, "warning") | is(coxfit_uc, "error"))){
          if(max(abs(summary(coxfit_uc)$coefficient[,1]))<15){
            aics[k]=extractAIC(coxfit_uc)[2]
          }
        }
      }
      k=k+1
    }
    if(sum(is.na(aics))==nuniv_cutoff){
      sel_genes=colnames(x_train)[ps_uc==min(ps_uc)]
      coxfit_uc=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x_train[,ps_uc==min(ps_uc)]), theta=1))
    }else{
      cut_sel=cutoff_grid[which(aics==min(aics, na.rm=T))[1]]
      sel_genes=colnames(x_train)[ps_uc<=cut_sel]
      coxfit_uc=coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,ps_uc<=cut_sel]))
    }
    b_uc=rep(0,p)
    b_uc[geneid %in% sel_genes]=summary(coxfit_uc)$coefficient[,1]
    
    coxfit_uc_test=coxph(Surv(data_test$t,data_test$delta)~as.matrix(x_test[,geneid %in% sel_genes]), init=summary(coxfit_uc)$coefficient[,1], control=coxph.control(iter.max=0))
    c_uc_test=summary(coxfit_uc_test)$concordance[1]
   

    lkhd_pc_sort=sort(lkhd_uc, decreasing=TRUE)
    lkhd_pc_thres=lkhd_pc_sort[round(min(ncase/4,p_sub/4))]
    sel_pc_genes=colnames(x_train)[lkhd_uc>=lkhd_pc_thres]
    x_pc_train=as.matrix(x_train[,lkhd_uc>=lkhd_pc_thres])
    inifit_pc=tryCatch(coxph(Surv(data_train$t, data_train$delta)~as.matrix(x_train[,lkhd_uc>=lkhd_pc_thres])), error=function(e) e, warning=function(w) w)
    if(is(inifit_pc, "warning") | is(inifit_pc, "error")){
      inifit_pc=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x_train[,lkhd_uc>=lkhd_pc_thres]), theta=1))
    }else{
      if(max(abs(coef(inifit_pc)))>10){
        inifit_pc=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x_train[,lkhd_uc>=lkhd_pc_thres]), theta=1))
      }
    }
    x_pc_test=as.matrix(x_test[,geneid %in% sel_pc_genes])

    
    # Lasso-penalized analysis with ComBat
    alpha_grid=seq(1,1,0.1)
    cvs_alpha=lambdas_alpha=c()
    for(alpha in alpha_grid){
      cv_l=cv.glmnet(as.matrix(x_train[,lkhd_uc>=lkhd_pc_thres]), Surv(data_train$t, data_train$delta), nfolds=nfold, family="cox", alpha=alpha, standardize=F)
      lambda_min=cv_l$lambda.min
      cvm_min=cv_l$cvm[cv_l$lambda==lambda_min]
      cvs_alpha=c(cvs_alpha, cvm_min)
      lambdas_alpha=c(lambdas_alpha, lambda_min)
    }
    lambda_sel=lambdas_alpha[which(cvs_alpha==min(cvs_alpha,na.rm=T))][1]
    alpha_sel=alpha_grid[which(cvs_alpha==min(cvs_alpha,na.rm=T))][1]
    glmnet_lc=glmnet(as.matrix(x_train[,lkhd_uc>=lkhd_pc_thres]), Surv(data_train$t, data_train$delta), family="cox", alpha=alpha_sel, lambda=lambda_sel, standardize=F)
    b_lc0temp=as.vector(glmnet_lc$beta)
    b_lc=rep(0,p)
    b_lc[geneid %in% sel_pc_genes]=b_lc0temp
    l_lc=lambda_sel
    coxfit_lc_test=coxph(Surv(data_test$t,data_test$delta)~x_pc_test, init=b_lc0temp, control=coxph.control(iter.max=0))
    c_lc_test=summary(coxfit_lc_test)$concordance[1]    
    
    
    # Adaptive Lasso-penalized analysis with ComBat
    w_ac=1/abs(coef(inifit_pc))
    
    cvs_alpha=lambdas_alpha=c()
    for(alpha in alpha_grid){
      cv_a=cv.glmnet(as.matrix(x_train[,lkhd_uc>=lkhd_pc_thres]), Surv(data_train$t, data_train$delta), nfolds=nfold, family="cox", alpha=alpha, standardize=F, penalty.factor=w_ac)
      lambda_min=cv_a$lambda.min
      cvm_min=cv_a$cvm[cv_a$lambda==lambda_min]
      cvs_alpha=c(cvs_alpha, cvm_min)
      lambdas_alpha=c(lambdas_alpha, lambda_min)
    }
    lambda_sel=lambdas_alpha[which(cvs_alpha==min(cvs_alpha,na.rm=T))][1]
    alpha_sel=alpha_grid[which(cvs_alpha==min(cvs_alpha,na.rm=T))][1]
    glmnet_ac=glmnet(as.matrix(x_train[,lkhd_uc>=lkhd_pc_thres]), Surv(data_train$t, data_train$delta), family="cox", alpha=alpha_sel, lambda=lambda_sel, standardize=F, penalty.factor=w_ac)
    b_ac0temp=as.vector(glmnet_ac$beta)
    b_ac=rep(0,p)
    b_ac[geneid %in% sel_pc_genes]=b_ac0temp
    l_ac=lambda_sel
    coxfit_ac_test=coxph(Surv(data_test$t,data_test$delta)~x_pc_test, init=b_ac0temp, control=coxph.control(iter.max=0))
    c_ac_test=summary(coxfit_ac_test)$concordance[1]
  }

  cat(" done \n")
  
  out=NULL
  out$call=match.call()
  
  if(BatMan_train==0 & cmbt_train==0){
    if(!(is.null(nonzero_beta_position))){
      out$bhat_oracle=b_o
    }
    out$bhat_univ=b_u
    out$bhat_lasso=b_l
    out$bhat_alasso=b_a
    if(!(is.null(nonzero_beta_position))){
      out$c_stats=cbind(c_o_test, c_u_test, c_l_test, c_a_test)
      colnames(out$c_stats)=c("Oracle","Univ","Lasso","Alasso")
    }else{
      out$c_stats=cbind(c_u_test, c_l_test, c_a_test)
      colnames(out$c_stats)=c("Univ","Lasso","Alasso")
    }
    out$lambdas=cbind(l_l, l_a)
    colnames(out$lambdas)=c("lambda_lasso","lambda_alasso")
    out$subset_size=subsize  
  }
    
  if(BatMan_train==1){
    if(!(is.null(nonzero_beta_position))){
      out$bhat_oracle_BatMan=b_os
    }
    out$bhat_univ_BatMan=b_us
    out$bhat_lasso_BatMan=b_ls
    out$bhat_alasso_BatMan=b_as
    if(!(is.null(nonzero_beta_position))){
      out$c_stats_BatMan=cbind(c_os_test, c_us_test, c_ls_test, c_as_test)
      colnames(out$c_stats_BatMan)=c("Oracle","Univ","Lasso","Alasso")
    }else{
      out$c_stats_BatMan=cbind(c_us_test, c_ls_test, c_as_test)
      colnames(out$c_stats_BatMan)=c("Univ","Lasso","Alasso")      
    }
    out$lambdas_BatMan=cbind(l_ls, l_as)
    colnames(out$lambdas_BatMan)=c("lambda_lasso","lambda_alasso")
  }
  
  if(cmbt_train==1){
    if(!(is.null(nonzero_beta_position))){
      out$bhat_oracle_ComBat=b_oc
    }
    out$bhat_univ_ComBat=b_uc
    out$bhat_lasso_ComBat=b_lc
    out$bhat_alasso_ComBat=b_ac
    if(!(is.null(nonzero_beta_position))){
      out$c_stats_ComBat=cbind(c_oc_test, c_uc_test, c_lc_test, c_ac_test)
      colnames(out$c_stats_ComBat)=c("Oracle","Univ","Lasso","Alasso")
    }else{
      out$c_stats_ComBat=cbind(c_uc_test, c_lc_test, c_ac_test)
      colnames(out$c_stats_ComBat)=c("Univ","Lasso","Alasso")      
    }
    out$lambdas_ComBat=cbind(l_lc, l_ac)
    colnames(out$lambdas_ComBat)=c("lambda_lasso","lambda_alasso")
    out$subset_size_ComBat=subsize_c
  }
  
  class(out)="BatMan.result"
  return(out)
}





