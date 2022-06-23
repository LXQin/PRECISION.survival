#' Simulate batch-effect-contaminated miRNA expression data and survival outcome based on the paired datasets from MSKCC.
#'
#' @param be.surv.data dataset containing true miRNA expression data (i.e. biological effects) and the survival outcome. Expression data starts at the fourth column.
#' @param he.data dataset containing pure handling effects and batch variables (e.g. technician, batch, slides).
#' @param nonzero_beta_position positions in true \eqn{\beta} that have nonzero values.
#' @param b0value vector of nonzero values corresponding to \code{nonzero_beta_position}.
#' @param batch_level at what level is batch defined. Either "slides" or "batches". The default is "slides".
#' @param t_sort_train 0 (default) if training samples are randomly assigned to microarray slides; 1 if training samples are sorted by survival time and then assigned to slides.
#' @param t_sort_test 0 (default) if test samples are randomly assigned to microarray slides; 1 if test samples are sorted by survival time and then assigned to slides.
#' @param t_rev 0 if training and test samples are sorted by survival time by the same direction; 1 if sorted by the opposite direction. Ignored if \code{t_sort_train=0} or \code{t_sort_test=0}.
#' @param he_train 0 if no handling effect in training data; 1 otherwise.
#' @param he_test 0 if no handling effect in test data; 1 otherwise.
#' @param he_train_scale multiplicative scalar applied to the original handling effect for training data (default=1).
#' @param he_test_scale multiplicative scalar applied to the original handling effect for test data (default=1).
#' @param seed seed for random number generation (default=123).
#' @description
#' \code{data_simulation} performs nonparametric simulation to generate training and test datasets from a pair of well-prepared MSKCC datasets with one containing 
#' true miRNA expression and the other containing pure batch effects. Batch-effect-contaminated miRNA expression data is simulated by virtual hybridization. Survival
#' outcome is simulated by permutation method. The sizes of the simulated training and test datasets are both 96.
#' @return A list that contains two datasets and two vectors of batch ids. The first element is the training dataset. The second element is the batch id vector
#' for the training dataset. The third element is the test dataset. The fourth element is the batch id vector for the test dataset.
#' @import "MASS"
#' @examples
#' ## Simulate a training dataset and a test dataset with six nonzero true beta values, 
#' ## where both training and test data contain batch effects 
#' ## that are not associated with the survival outcome.
#' sim1=data_simulation(be.survival, he, nonzero_beta_position=c(385,866,1010,2218,2660,3026), 
#' b0value=c(1.04,1.40,1.45,1.62,3.13,1.76), t_sort_train=0, t_sort_test=0, 
#' t_rev=0, he_train=1, he_test=1, he_train_scale=1, he_test_scale=1, seed=123)
#' @references Ni, A., Qin, L. (2021) "Performance Evaluation of Transcriptomics Data Normalization for Survival Risk Prediction". \emph{Briefings in Bioinformatics}, 22(6), bbab257.
#' @export

data_simulation=function(be.surv.data, he.data, nonzero_beta_position, b0value, batch_level="slides", t_sort_train=0, t_sort_test=0, t_rev, 
                         he_train, he_test, he_train_scale=1, he_test_scale=1, seed=123){
  
  if(!(batch_level %in% c("slides","batches"))){
    stop("Invaid batch_level. Must be either 'slides' or 'batches'")
  }
  be.PFS2=be.surv.data[order(be.surv.data$t),]
  x=as.matrix(be.PFS2[,4:ncol(be.PFS2)])
  p=ncol(x)
  n=nrow(x)
  n_test=n
  b0=rep(0,p)
  b0[nonzero_beta_position]=b0value
  xb=x%*%b0
  geneid=colnames(x)
  
  set.seed(seed)

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
  hedata_train[,2:(p+1)]=hedata_train[,2:(p+1)]*he_train_scale
  hedata_test=rbind(hedata_test1, hedata_test2)
  hedata_test[,2:(p+1)]=hedata_test[,2:(p+1)]*he_test_scale
  
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
    if(batch_level=="slides"){
      batch_id_train=hedata_train$slide
    }
    if(batch_level=="batches"){
      batch_id_train=hedata_train$batch
    }    
  }
  if(he_train==0){
    data_train_full=data_train[sample(1:n),]
    if(batch_level=="slides"){
      batch_id_train=hedata_train$slide
    }
    if(batch_level=="batches"){
      batch_id_train=hedata_train$batch
    }
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
    if(batch_level=="slides"){
      batch_id_test=hedata_test$slide
    }
    if(batch_level=="batches"){
      batch_id_test=hedata_test$batch
    }
  }
  if(he_test==0){
    data_test=data_test[sample(1:n_test),]
    if(batch_level=="slides"){
      batch_id_test=hedata_test$slide
    }
    if(batch_level=="batches"){
      batch_id_test=hedata_test$batch
    }
  }
  
  out=list(data_train_full, batch_id_train, data_test, batch_id_test)
  return(out)
}





