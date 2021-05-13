#' Comparison of prediction accuracy of normalization methods in survival analysis.
#'
#' @param be.surv.data dataset containing true miRNA expression data (i.e. biological effects) and the survival outcome.
#' @param he.data dataset containing pure handling effects and batch variables (e.g. technician, batch, slides).
#' @param nonzero_position positions in true \eqn{\beta} that have nonzero values.
#' @param b0value vector of nonzero values corresponding to \code{nonzero_position}.
#' @param ps_rho_cutoff cutoff of correlation coefficient during pre-screen (default=0.9).
#' @param t_sort_train 0 (default) if training samples are randomly assigned to microarray slides; 1 if training samples are sorted by survival time and then assigned to slides.
#' @param t_sort_test 0 (default) if test samples are randomly assigned to microarray slides; 1 if test samples are sorted by survival time and then assigned to slides.
#' @param t_rev 0 if training and test samples are sorted by survival time by the same direction; 1 if sorted by the opposite direction. Ignored if \code{t_sort_train=0} or \code{t_sort_test=0}.
#' @param he_train 0 if no handling effect in training data; 1 otherwise.
#' @param he_test 0 if no handling effect in test data; 1 otherwise.
#' @param he_train_scale multiplicative scalar applied to the original handling effect for training data (default=1).
#' @param he_test_scale multiplicative scalar applied to the original handling effect for test data (default=1).
#' @param norm_type normalization methods (0: no normalization; 1: median normalization; 2: quantile normalization; 3: Variance stabilizing transformation (VSN).
#' @param norm_train 0 if no normalization applied to training data; 1 otherwise. Ignored if \code{norm_type=0}.
#' @param norm_test 0 if no normalization applied to test data; 1 otherwise. Ignored if \code{norm_type=0}.
#' @param univ_cutoff_range range of cutoffs of p values for univariate analysis (default=\code{c(0,0.05)}).
#' @param nuniv_cutoff number of equal-spaced cutoff values in \code{univ_cutoff_range} (default=20).
#' @param lambda_max_glmnet maximum value of grid-search range of tuning parameter lambda in lasso and adaptive lasso (default=0.3).
#' @param nlambda number of equal-sapced lambda values in \code{lambda_max_glmnet} (default=30).
#' @param nfold number of folds for the cross-validation for tuning parameter selection in lasso and adaptive lass (default=6).
#' @param nsim number of simulation replicates (default=400).
#' @param seed seed for random number generation (default=123).
#' @description
#' \code{simulation} performs nonparametric simulation to evaluate prediction accuracy of three normalization methods (median normalization, quantile normalization,
#' Variance stabilizing transformation) using a pair of well-prepared datasets with one containing true miRNA expression and the other containing pure handling effects.
#' @return An object of class \code{sim.result} containing estimated regression coefficients and C-indices from each method.
#'   \item{bhat_oracle}{estimated regression coefficients from Oracle method (each row represents one simulation replicate; each column represents one marker)}
#'   \item{bhat_univ}{estimated regression coefficients from univariate method}
#'   \item{bhat_lasso}{estimated regression coefficients from lasso method}
#'   \item{bhat_alasso}{estimated regression coefficients from adaptive lasso method}
#'   \item{c_stats}{C-indices from all four methods based on test data (each column represents one method)}
#'   \item{lambdas}{selected tuning parameter lambda for Lasso and adaptive Lasso methods}
#'   \item{subset_size}{the sizes of data used for variable selection after pre-screening}
#' @import "MASS","survival","preprocessCore","glmnet","sva","vsn"
#' @examples
#' ## conduct a simulation with six nonzero true beta values, 
#' ## where both training and test data contain handling effects 
#' ## that are not associated with the survival outcome, and both 
#' ## are subject to median normalization.
#' sim1=simulation(be.surv, he, nonzero_position=c(385,866,1010,2218,2660,3026), 
#' b0value=c(1.04,1.40,1.45,1.62,3.13,1.76), t_sort_train=0, t_sort_test=0, 
#' t_rev=0, he_train=1, he_test=1, norm_type=1, norm_train=1, norm_test=1)
#' @references Ni, A., Qin, L. (2021+) "Performance Evaluation of Transcriptomics Data Normalization for Survival Risk Prediction". accepted by \emph{Briefings in Bioinformatics}.
#' @export


survival.simulation=function(be.surv.data, he.data, nonzero_position, b0value, ps_rho_cutoff=0.9, t_sort_train=0, t_sort_test=0, t_rev, he_train, he_test, he_train_scale=1, he_test_scale=1, 
                         norm_type, norm_train, norm_test, univ_cutoff_range=c(0,0.05), nuniv_cutoff=20, lambda_max_glmnet=0.3, nlambda=30, nfold=6, nsim=400, seed=123){

  he=he.data[order(he.data$slide),]

  be.PFS2=be.surv.data[order(be.surv.data$t),]
  x=as.matrix(be.PFS2[,4:ncol(be.PFS2)])
  p=ncol(x)
  n=nrow(x)
  b0=rep(0,p)
  b0[nonzero_position]=b0value
  xb=x%*%b0
  ncase=sum(be.PFS2$delta)
  n_test=n
  geneid=colnames(x)
  
  set.seed(seed)
  
  # letter after "_": o - oracle analysis, u - univariate analysis, l - lasso-penalized reg, a - adaptive lasso-penalized reg
  
  b_o=b_u=b_l=b_a=matrix(0, nrow=nsim, ncol=p)
  c_o_test=c_u_test=c_l_test=c_a_test=l_l=l_a=al_l=al_a=rep(0, nsim)
  subsize=rep(0, nsim)
  
  cat("Iteration ")
  for(iter in 1:nsim){
    ## permute covariates based on hazards to induce association between six genes and the survival outcome
    ## training data
    x2=x
    xb2=xb
    x_train=matrix(NA,nrow=n,ncol=p)
    for(i in 1:(n-1)){
      ps=exp(xb2)/sum(exp(xb2))
      sel=rmultinom(1, 1, ps)
      x_train[i,]=x2[sel==1,]
      xb2=xb2[sel!=1]
      x2=x2[sel!=1,]
    }
    x_train[n,]=x2
    colnames(x_train)=geneid
    data_train=cbind(be.PFS2[,1:3],x_train)
    
    ## test data
    x2=x
    xb2=xb
    x_test=matrix(NA,nrow=n,ncol=p)
    for(i in 1:(n-1)){
      ps=exp(xb2)/sum(exp(xb2))
      sel=rmultinom(1, 1, ps)
      x_test[i,]=x2[sel==1,]
      xb2=xb2[sel!=1]
      x2=x2[sel!=1,]
    }
    x_test[n,]=x2
    colnames(x_test)=geneid
    data_test=cbind(be.PFS2[,1:3],x_test)
    
    
    # assign handling effects to covariates
    
    he1=he[he$technician==1,]
    he2=he[he$technician==2,]
    hedata_train1=he1[he1$slide %in% 1:5,] 
    hedata_test1=he1[he1$slide %in% 6:10,]
    hedata_train2=he2[he2$slide %in% 11:17,] 
    hedata_test2=he2[he2$slide %in% 18:24,]    
    hedata_train=rbind(hedata_train1, hedata_train2)
    hedata_train[,2:3524]=hedata_train[,2:3524]*he_train_scale
    hedata_test=rbind(hedata_test1, hedata_test2)
    hedata_test[,2:3524]=hedata_test[,2:3524]*he_test_scale
    
    if(he_train==1){
      if(t_sort_train==0){
        data_train_full=data_train[sample(1:n),]
      }
      if(t_sort_train==1){
        train_ns=table(hedata_train$technician)
        data_train1=data_train[1:train_ns[1],]
        data_train1=data_train1[sample(1:nrow(data_train1)),]
        data_train2=data_train[(train_ns[1]+1):n,]
        data_train2=data_train2[sample(1:nrow(data_train2)),] 
        data_train_full=rbind(data_train1, data_train2)
      }
      data_train_full[,4:(p+3)]=data_train_full[,4:(p+3)]+hedata_train[,2:(p+1)]
      data_train_full$strat8id=hedata_train$slide
      data_train_full$stratbid=hedata_train$batch
      x_train_full=data_train_full[,4:(p+3)]
      data_train_full=cbind(data_train_full[,c(1:3,ncol(data_train_full)-1,ncol(data_train_full))],x_train_full)
    }
    if(he_train==0){
      data_train_full=data_train[sample(1:n),]
      data_train_full$strat8id=hedata_train$slide
      data_train_full$stratbid=hedata_train$batch
      x_train_full=data_train_full[,4:(p+3)]
      data_train_full=cbind(data_train_full[,c(1:3,ncol(data_train_full)-1,ncol(data_train_full))],x_train_full)
    }

    if(he_test==1){
      if(t_sort_test==0){
        data_test=data_test[sample(1:n_test),]
      }
      if(t_sort_test==1){
        test_ns=table(hedata_test$technician)
        if(t_rev==0){
          data_test1=data_test[1:test_ns[1],]
          data_test1=data_test1[sample(1:nrow(data_test1)),]
          data_test2=data_test[(test_ns[1]+1):n_test,]
          data_test2=data_test2[sample(1:nrow(data_test2)),] 
        }        
        if(t_rev==1){
          data_test1=data_test[(test_ns[1]+1):n_test,]
          data_test1=data_test1[sample(1:nrow(data_test1)),]
          data_test2=data_test[1:test_ns[1],]
          data_test2=data_test2[sample(1:nrow(data_test2)),]           
        }
        data_test=rbind(data_test1, data_test2)
      }
      data_test[,4:(p+3)]=data_test[,4:(p+3)]+hedata_test[,2:(p+1)]
      data_test$strat8id=hedata_test$slide
      data_test$stratbid=hedata_test$batch
      x_test=data_test[,4:(p+3)]
      data_test=cbind(data_test[,c(1:3,ncol(data_test)-1,ncol(data_test))],x_test)
    }
    if(he_test==0){
      data_test=data_test[sample(1:n_test),]
      data_test$strat8id=hedata_test$slide
      data_test$stratbid=hedata_test$batch
      x_test=data_test[,4:(p+3)]
      data_test=cbind(data_test[,c(1:3,ncol(data_test)-1,ncol(data_test))],x_test)
    }
    
    x_train_full_raw=x_train_full
    data_train_full_raw=data_train_full
    x_test_raw=x_test
    data_test_raw=data_test
  

    ####### normalization ######
    
    if(norm_type==1){
      if(norm_train==1){
        x_train_full_rowmedian=apply(x_train_full_raw,1,median)
        x_train_full_rowmedian_diff=x_train_full_rowmedian-mean(x_train_full_rowmedian)
        x_train_full=as.matrix(x_train_full_raw-x_train_full_rowmedian_diff%*%t(rep(1,ncol(x_train_full_raw))))
      }
      if(norm_test==1){
        x_test_rowmedian=apply(x_test_raw,1,median)
        x_test_rowmedian_diff=x_test_rowmedian-mean(x_train_full_rowmedian)
        x_test=as.matrix(x_test_raw-x_test_rowmedian_diff%*%t(rep(1,ncol(x_test_raw))))
      }
    }
        
    if(norm_type==2){
      if(norm_train==1){
        x_train_full=t(normalize.quantiles(t(x_train_full_raw)))
        x_quantile=sort(x_train_full[1,])           
      }
      if(norm_test==1){
        x_test=t(apply(x_test_raw, 1, function(u){x_quantile[rank(u)]}))
      }
    }
        
    if(norm_type==3){
      if(norm_train==1){
        data.vsn=vs.norm(train=t(x_train_full_raw))
        x_train_full=t(data.vsn[[1]])
      }
      if(norm_test==1){
        x_test=t(vs.norm(test=t(x_test_raw), ref.dis=data.vsn$ref.dis)[[2]])
      }
    }

    colnames(x_train_full)=colnames(x_train_full_raw)
    data_train_full=cbind(data_train_full_raw[,1:5],x_train_full) 
    colnames(x_test)=colnames(x_test_raw)
    data_test=cbind(data_test_raw[,1:5],x_test)

    #### pre-screen

    means_train=colMeans(x_train_full)
    x_train_sub=x_train_full[,means_train>8]
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
    subsize[iter]=p_sub
    data_train=cbind(data_train_full[,1:5],x_train)


    ################ analysis ##############

    # oracle analysis
    x0_train=x_train_full[,b0!=0]
    coxfit_o=tryCatch(coxph(Surv(data_train$t,data_train$delta)~as.matrix(x0_train)), error=function(e) e, warning=function(w) w)
    if(is(coxfit_o, "warning") | is(coxfit_o, "error")){
      coxfit_o=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x0_train), theta=1))
    }
    b_o0=rep(0,p)
    b_o0[b0!=0]=summary(coxfit_o)$coefficient[,1]
    b_o[iter,]=b_o0

    x0_test=x_test[,b0!=0]
    coxfit_o_test=coxph(Surv(data_test$t,data_test$delta)~as.matrix(x0_test), init=summary(coxfit_o)$coefficient[,1], control=coxph.control(iter.max=0))
    c_o_test[iter]=summary(coxfit_o_test)$concordance[1]

  
    # univariate analysis with handling effect
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
        coxfit_u=tryCatch(coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,selected==1])),
                          error=function(e) e, warning=function(w) w)
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
        coxfit_u=tryCatch(coxph(Surv(data_train$t,data_train$delta)~ridge(x_selu, theta=1)),
                          error=function(e) e, warning=function(w) w)
        while(is(coxfit_u, "warning") | is(coxfit_u, "error")){
          sel_genes=sel_genes[-which(pselu==max(pselu))]
          x_selu=as.matrix(x_train[,sel_genes])
          pselu=pselu[-which(pselu==max(pselu))]
          coxfit_u=tryCatch(coxph(Surv(data_train$t,data_train$delta)~ridge(x_selu, theta=1)),
                            error=function(e) e, warning=function(w) w)
        }
      }
    }else{
      cut_sel=cutoff_grid[which(aics==min(aics, na.rm=T))[1]]
      sel_genes=colnames(x_train)[ps_u<=cut_sel]
      coxfit_u=coxph(Surv(data_train$t,data_train$delta)~as.matrix(x_train[,ps_u<=cut_sel]))
    }
    b_u0=rep(0,p)
    b_u0[geneid %in% sel_genes]=summary(coxfit_u)$coefficient[,1]
    b_u[iter,]=b_u0
    coxfit_u_test=coxph(Surv(data_test$t,data_test$delta)~as.matrix(x_test[,geneid %in% sel_genes]), init=summary(coxfit_u)$coefficient[,1], control=coxph.control(iter.max=0))
    c_u_test[iter]=summary(coxfit_u_test)$concordance[1]


    lkhd_p_sort=sort(lkhd_u, decreasing=TRUE)
    lkhd_p_thres=lkhd_p_sort[round(min(ncase/4,p_sub/4))]
    sel_p_genes=colnames(x_train)[lkhd_u>=lkhd_p_thres]
    x_p_train=as.matrix(x_train[,lkhd_u>=lkhd_p_thres])
    inifit_p=tryCatch(coxph(Surv(data_train$t, data_train$delta)~as.matrix(x_train[,lkhd_u>=lkhd_p_thres])), error=function(e) e, warning=function(w) w)
    if(is(inifit_p, "warning") | is(inifit_p, "error") | max(abs(coef(inifit_p)))>10){
      inifit_p=coxph(Surv(data_train$t,data_train$delta)~ridge(as.matrix(x_train[,lkhd_u>=lkhd_p_thres]), theta=1))
    }
    x_p_test=as.matrix(x_test[,geneid %in% sel_p_genes])

    # Lasso-penalized analysis with handling effect
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
    b_l0=rep(0,p)
    b_l0[geneid %in% sel_p_genes]=b_l0temp
    b_l[iter,]=b_l0
    l_l[iter]=lambda_sel
    al_l[iter]=alpha_sel

    coxfit_l_test=coxph(Surv(data_test$t,data_test$delta)~x_p_test, init=b_l0temp, control=coxph.control(iter.max=0))
    c_l_test[iter]=summary(coxfit_l_test)$concordance[1]


    # Adaptive Lasso-penalized analysis with handling effect
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
    b_a0=rep(0,p)
    b_a0[geneid %in% sel_p_genes]=b_a0temp
    b_a[iter,]=b_a0
    l_a[iter]=lambda_sel
    al_a[iter]=alpha_sel

    coxfit_a_test=coxph(Surv(data_test$t,data_test$delta)~x_p_test, init=b_a0temp, control=coxph.control(iter.max=0))
    c_a_test[iter]=summary(coxfit_a_test)$concordance[1]
    
    cat(paste(iter," "))
  }
  cat(" done \n")
  
  out=NULL
  out$call=match.call()
  out$bhat_oracle=b_o
  out$bhat_univ=b_u
  out$bhat_lasso=b_l
  out$bhat_alasso=b_a
  c_stats=cbind(c_o_test, c_u_test, c_l_test, c_a_test)
  colnames(c_stats)=c("c_oracle","c_univ","c_lasso","c_alasso")
  out$c_stats=c_stats
  lambdas=cbind(l_l, l_a)
  colnames(lambdas)=c("lambda_lasso","lambda_alasso")
  out$lambdas=lambdas
  out$subset_size=subsize
  class(out)="sim.result"
  return(out)
}






