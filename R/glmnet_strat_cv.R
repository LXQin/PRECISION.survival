#' Perform GLMnet for stratified Cox model

"glmnet_strat_cv"=function(data, lambda_max, nlambda, nfold, penalty_wt, threshold){
  data=data[order(data[,4]),]
  options(warn=-1)
  temp=as.vector(unlist(by(data[,c("t","delta")], data[,4], function(x){x[,1]<min(x[x[,2]==1,1])})))
  options(warn=0)
  data=data[!temp,]
  n=nrow(data)
  x=as.matrix(data[,5:ncol(data)])
  p=ncol(x)
  te=data$t[data$delta==1]
  N=length(te)
  I0full=ifelse(data$t%*%t(rep(1,N))>=rep(1,n)%*%t(te), 1, 0)
  alpha_grid=seq(1,1,0.1)
  # alpha_grid=1
  lambda_grid=lambda_max*0.05^((0:nlambda)/nlambda)
  stratid=data[,4]
  fold_size=length(unique(stratid))/nfold
  ustratid=unique(stratid)
  cvs_alpha=c()
  lambdas_alpha=c()
  for(alpha in alpha_grid){
    bint=rep(0,p)
    cvs=c()
    for(lambda in lambda_grid){
      flag=0
      cvl=rep(NA,nfold)
      dev_conv=rep(0,nfold)
      for(f in 1:nfold){
        data_traincv=data[!(stratid %in% ustratid[((1:fold_size)+fold_size*(f-1))]),]
        n_traincv=nrow(data_traincv)
        te_traincv=data_traincv$t[data_traincv$delta==1]
        N_traincv=length(te_traincv)
        x_traincv=as.matrix(data_traincv[,5:ncol(data_traincv)])
        x_meanscv=colMeans(x_traincv)
        x_traincv=x_traincv-rep(1,n_traincv)%*%t(x_meanscv)
        data_traincv=cbind(data_traincv[,1:4],x_traincv)
        
        diff=1
        step=0
        while(diff>threshold){
          bintold=bint
          eta_traincv=x_traincv%*%bint
          el=ws=rep(NA,n_traincv)
          stratid_train=data_traincv[,4]
          for(s in unique(stratid_train)){
            datas_traincv=data_traincv[stratid_train==s,]
            ns_traincv=nrow(datas_traincv)
            etas_traincv=eta_traincv[stratid_train==s]
            ts_traincv=datas_traincv$t
            tes_traincv=ts_traincv[datas_traincv$delta==1]
            Ns_traincv=length(tes_traincv)
            I0s_traincv=ifelse(ts_traincv%*%t(rep(1,Ns_traincv))>=rep(1,ns_traincv)%*%t(tes_traincv), 1, 0)  ## risk set
            exp_etas_traincv=exp(etas_traincv%*%t(rep(1,Ns_traincv)))
            exp_etas_sum_traincv=colSums(exp_etas_traincv*I0s_traincv)
            exp_etas_sum_n_traincv=rep(1,ns_traincv)%*%t(exp_etas_sum_traincv)
            el[stratid_train==s]=datas_traincv$delta-rowSums(exp_etas_traincv*I0s_traincv/exp_etas_sum_n_traincv)
            exp_etas_sum2_traincv=exp_etas_sum_traincv^2
            exp_etas_sum2_n_traincv=rep(1,ns_traincv)%*%t(exp_etas_sum2_traincv)
            ws[stratid_train==s]=-rowSums((exp_etas_traincv*I0s_traincv*exp_etas_sum_n_traincv-(exp_etas_traincv*I0s_traincv)^2)/exp_etas_sum2_n_traincv)
          }
          z=eta_traincv-el/ws
          change=1
          iter=0
          while(max(change)>threshold){
            change=rep(NA,p)
            for(k in 1:p){
              xk=x_traincv[,k]
              xnk=x_traincv[,-k]
              numer=mean(ws*xk*(z-xnk%*%bint[-k]))
              S=sign(numer)*ifelse((abs(numer)-lambda*alpha*penalty_wt[k])>0, abs(numer)-lambda*alpha*penalty_wt[k], 0)
              bhatk=S/(mean(ws*xk^2)+lambda*(1-alpha))
              change[k]=abs(bint[k]-bhatk)
              bint[k]=bhatk
            }
            iter=iter+1
            # print(max(change))
          # }
            # print("p ")
            if(abs(max(bint))>30 | NaN %in% bint | NA %in% bint | Inf %in% bint | -Inf %in% bint | iter>200){
              flag=1
              # cat(paste("flag=", flag,sep=""))
              break
            }
            # cat(paste("flag=", flag,sep=""))
          }
          if(flag==1){break}
          diff=max(abs(bintold-bint))
          step=step+1
          if(step>200){
            flag=1
            break
          }
        } # end of while loop
        if(flag==1){break}
        xbeta_traincv=x_traincv%*%bint
        expx_traincv=exp(xbeta_traincv)
        lkhd_traincv_term2=term2_null=0
        for(s in unique(stratid_train)){
          datas_traincv=data_traincv[stratid_train==s,]
          ns_traincv=nrow(datas_traincv)
          ts_traincv=datas_traincv$t
          tes_traincv=ts_traincv[datas_traincv$delta==1]
          Ns_traincv=length(tes_traincv)
          I0s_traincv=ifelse(ts_traincv%*%t(rep(1,Ns_traincv))>=rep(1,ns_traincv)%*%t(tes_traincv), 1, 0)  ## risk set
          lkhd_traincv_term2=lkhd_traincv_term2+sum(log(colSums((expx_traincv[stratid_train==s]%*%t(rep(1,Ns_traincv)))*I0s_traincv)))
          term2_null=term2_null+sum(log(colSums(I0s_traincv)))
        }
        lkhd_traincv=sum(xbeta_traincv[data_traincv$delta==1])-lkhd_traincv_term2
        D_traincv=-2*lkhd_traincv
        D_nulltraincv=2*term2_null
        dev_conv[f]=1*((D_nulltraincv-D_traincv)>=0.99*D_nulltraincv)
        
        x_stdcv=x-rep(1,n)%*%t(x_meanscv)
        xbeta=x_stdcv%*%bint
        expx_full=exp(xbeta)
        lkhd_full_term2=0
        for(s in unique(stratid)){
          datas=data[stratid==s,]
          ns=nrow(datas)
          ts=datas$t
          tes=ts[datas$delta==1]
          Ns=length(tes)
          I0s=ifelse(ts%*%t(rep(1,Ns))>=rep(1,ns)%*%t(tes), 1, 0)  ## risk set
          lkhd_full_term2=lkhd_full_term2+sum(log(colSums((expx_full[stratid==s]%*%t(rep(1,Ns)))*I0s)))        
        }
        lkhd_full=sum(xbeta[data$delta==1])-lkhd_full_term2
        cvi=lkhd_full-lkhd_traincv
        cvl[f]=cvi
        # cat(f)
      } # end for nfold loop 
      if(flag==1){
        cvs=c(cvs, NA)
        # cat("lambda=NA ")
        # cat(paste("flag=",flag,sep=""))
        bint=rep(0,p)
      }else{
        cvs=c(cvs, mean(cvl))
        # cat(paste("lambda=",lambda," "))        
      }
      if(sum(dev_conv)==nfold){break}
    } # end lambda loop
    if(sum(is.na(cvs))==length(cvs)){
      cvs_alpha=c(cvs_alpha, NA)
      lambdas_alpha=c(lambdas_alpha, NA)
    }else{
      lambda0=lambda_grid[which(cvs==max(cvs,na.rm=T))][1]  
      cvs_alpha=c(cvs_alpha, max(cvs,na.rm=T))
      lambdas_alpha=c(lambdas_alpha, lambda0)      
    }
    # cat(paste("alpha=",alpha," "))
  } # end alpha loop
  
  lambda=lambdas_alpha[which(cvs_alpha==max(cvs_alpha,na.rm=T))][1]
  alpha=alpha_grid[which(cvs_alpha==max(cvs_alpha,na.rm=T))][1]
  bint=rep(0,p)
  
  diff=1
  x_means=colMeans(x)
  x_std=x-rep(1,n)%*%t(x_means)
  while(diff>threshold){
    bintold=bint
    eta=x_std%*%bint
    el=ws=rep(NA,n)
    for(s in unique(stratid)){
      datas=data[stratid==s,]
      ns=nrow(datas)
      etas=eta[stratid==s]
      ts=datas$t
      tes=ts[datas$delta==1]
      Ns=length(tes)
      I0s=ifelse(ts%*%t(rep(1,Ns))>=rep(1,ns)%*%t(tes), 1, 0)  ## risk set
      exp_etas=exp(etas%*%t(rep(1,Ns)))
      exp_etas_sum=colSums(exp_etas*I0s)
      exp_etas_sum_n=rep(1,ns)%*%t(exp_etas_sum)
      el[stratid==s]=datas$delta-rowSums(exp_etas*I0s/exp_etas_sum_n)
      exp_etas_sum2=exp_etas_sum^2
      exp_etas_sum2_n=rep(1,ns)%*%t(exp_etas_sum2)
      ws[stratid==s]=-rowSums((exp_etas*I0s*exp_etas_sum_n-(exp_etas*I0s)^2)/exp_etas_sum2_n)
    }
    z=eta-el/ws
    change=1
    while(max(change)>threshold){
      change=rep(NA,p)
      for(k in 1:p){
        xk=x_std[,k]
        xnk=x_std[,-k]
        numer=mean(ws*xk*(z-xnk%*%bint[-k]))
        S=sign(numer)*ifelse((abs(numer)-lambda*alpha*penalty_wt[k])>0, abs(numer)-lambda*alpha*penalty_wt[k], 0)
        bhatk=S/(mean(ws*xk^2)+lambda*(1-alpha))
        change[k]=abs(bint[k]-bhatk)
        bint[k]=bhatk
      }
    }
    diff=max(abs(bintold-bint))
  }
  
  return(list(bint, lambda, alpha))
}