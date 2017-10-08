library("SQUAREM")

source('../code/lik.R')
source('../code/set_data.R')
source('../code/ash.R')
source('../code/mix.R')
source('../code/normalmix.R')
source('../code/mix_opt.R')
source('../code/output.R')
source('../code/ashutility.R')

generic.cormotif=function(betahat,sebetahat,mixcompdist = "normal",optmethod = "w_mixEM",df = NULL
       ,nullweight = 10,pointmass = TRUE,mixsd = NULL,gridmult = sqrt(2),outputlevel = 2,g = NULL
       ,fixg = FALSE,mode = 0,alpha = 0,grange = c(-Inf,Inf),control = list(),lik = NULL
       ,weights=NULL,pi_thresh = 1e-10,K=1,BIC=TRUE, mess = TRUE) {
  
  n = dim(betahat)[1]
  R = dim(betahat)[2]
  
  # set likelihood based on defaults if missing
  if(is.null(lik)){ 
    if(is.null(df)){
      lik = lik_normal()
    } else {lik = lik_t(df)}
  }
  
  ## determine mixsd
  betahat_v = as.vector(betahat)
  sebetahat_v = as.vector(sebetahat)
  data_v = set_data(betahat_v, sebetahat_v, lik, alpha)
  if(is.null(mixsd)){
    mixsd = autoselect.mixsd(data_v,gridmult,mode,grange,mixcompdist)
  }
  if(pointmass){
    mixsd = c(0,mixsd)
  }
  L = length(mixsd)
  null.comp = which.min(mixsd)
  w = initpi(L,n,null.comp)
  if(mixcompdist == "normal") g=normalmix(w,rep(mode,L),mixsd)
  
  # compute conv_dens matrix
  matrix_llik = list()
  data = list()
  for(r in 1:R){
    data[[r]] = set_data(betahat[,r], sebetahat[,r], lik, alpha)
    matrix_llik[[r]] =  t(log_comp_dens_conv(g,data[[r]]))
  }
  
  ## fit the model via EM 
  fitresult<-list()
  for(i in 1:length(K)){
    fitresult[[i]] = gcmfit(matrix_llik,K=K[i],w_init=w,pi_thresh=pi_thresh,mess=mess)
  }
  
  bic<-rep(0,length(K))
  aic<-rep(0,length(K))
  loglike<-rep(0,length(K))
  for(i in 1:length(K)){
    bic[i] = fitresult[[i]]$BIC
    aic[i] = fitresult[[i]]$AIC
    loglike[i] = fitresult[[i]]$loglike
  }
  if(BIC==TRUE){
    bestflag=which(bic==min(bic))
  }else{
    bestflag=which(aic==min(aic))
  }
  bestmotif=fitresult[[bestflag]]
  ## lfdr and lsdr
  lfdr = list()
  lfsr = list()
  post_mean = list()
  for(i in 1:length(K)){
    ds = comp_lfdsr(fitresult[[i]], g, data)
    lfdr[[i]] = ds$lfdr
    lfsr[[i]] = ds$lfsr
    post_mean[[i]] = ds$post_mean
  }
  #ds = comp_lfdsr(bestmotif, g, data)
  # output = set_output(outputlevel) #sets up flags for what to output
  # resfns = set_resfns(output)
  # result = list()
  # if(length(resfns)>0){
  #   for(r in 1:R){
  #     data = set_data(betahat[,r], sebetahat[,r], lik, alpha)
  #     what = t(bestmotif$W[[r]])%*%bestmotif$pi
  #     if(mixcompdist == "normal") ghat=normalmix(what,rep(mode,L),mixsd)
  #     result[[r]] = as.data.frame(lapply(resfns,do.call,list(g=prune(ghat,pi_thresh),data=data)))
  #   }
  #  
  # }
  return(list(bestmotif=fitresult[[bestflag]],bic=cbind(K,bic),
                aic=cbind(K,aic),loglike=cbind(K,loglike),mixsd=mixsd,allmotif=fitresult,lfdr=lfdr,lfsr=lfsr,post_mean=post_mean))
}

gcmfit = function(matrix_llik,K=1,max.iter=300,tol = 1e-4,w_init,pi_thresh,mess=TRUE){
  R = length(matrix_llik)
  n =  dim(matrix_llik[[1]])[1]
  L = dim(matrix_llik[[1]])[2]
  
  ## initialize
  p = rep(1,K)/K
  W = list()
  for(r in 1:R){
    W[[r]] = MCMCpack::rdirichlet(K,alpha=rep(1,L))
    #W[[r]] = matrix(1/L,nrow=K,ncol=L)
    #W[[r]] = matrix(rep(w_init,each=K),nrow=K)
    W[[r]][1,] = w_init
  }
  
  for(i.iter in 1:max.iter) {
    if(mess){
      if((i.iter%%5) == 0) {
        print(paste("We have run the first ", i.iter, " iterations for K=", K,sep=""))
        #print(loglike.old)
      }
    }
    clustlike<-matrix(0,n,K)
    templike <- matrix(0,n,L)
    for(k in 1:K) {
      for(r in 1:R) {
        templike = t(log(W[[r]][k,])+t(matrix_llik[[r]]))
        tempmax<-apply(templike,1,max)
        templike<-exp(templike-tempmax)
        tempsum<-apply(templike,1,sum)
        clustlike[,k]<-clustlike[,k]+tempmax+log(tempsum)
      }
      #clustlike[,k]<-clustlike[,k]+log(p[k])
    }
    ## update pi
    matrix_llik2 = clustlike
    tempmax = apply(matrix_llik2,1,max)
    matrix_lik2 = exp(matrix_llik2-tempmax)
    fit.p=mixEM(matrix_lik=matrix_lik2, prior=1,pi_init = p)
    p.new = normalize(fit.p$pihat + 1e-8)
    #print(p.new)
    
    ## update W
    ### compute p_{jk}
    for(k in 1:K){
      clustlike[,k]<-clustlike[,k]+log(p[k])
    }
    tempmax2 = apply(clustlike,1,max)
    clustlike<-exp(clustlike-tempmax2)
    tempsum<-apply(clustlike,1,sum)
    for(k in 1:K) {
      clustlike[,k]<-clustlike[,k]/tempsum
    }
    clustlike[clustlike<pi_thresh] = 0
    for(j in 1:n){
      clustlike[j,] = normalize(clustlike[j,])
    }
    W.new = W
    for(k in 1:K){
      for(r in 1:R){
        if(sum(clustlike[,k])>0){
          ww = normalize((clustlike[,k]))
          #print((clustlike[,k]))
          matrix_llik3 = matrix_llik[[r]]
          lnorm = apply(matrix_llik3,1,max)
          matrix_lik3 = exp(matrix_llik3 - lnorm)
          fit.w = w_mixEM(matrix_lik=matrix_lik3, prior=1,pi_init = W[[r]][k,],weights = ww)
          W.new[[r]][k,] = fit.w$pihat
        }
      }
    }
    ## compute error
    err.p = max(abs(p.new-p))
    err.w = rep(0,R)
    for(r in 1:R){
      err.w[r] = max(abs(W.new[[r]]-W[[r]]))
    }
    err.w1 = max(err.w)
    err = max(err.p,err.w1)
    
    p = p.new
    W = W.new
    
    if(err < tol) {
      break;
    }
  }
  ## compute log-likelihood
  loglike = penloglik(p,matrix_lik2,prior=1)+sum(tempmax)
  BIC = -2*loglike+(K*R*(L-1)+K-1)*log(n)
  AIC = -2*loglike+2*(K*R*(L-1)+K-1)
  return(list(K = K,pi=p,W=W,loglike=loglike,BIC=BIC,AIC=AIC,clustlike=clustlike))
}

comp_lfdsr= function(bestmotif, g, data){
  R = length(data)
  n = length(data[[1]]$x)
  L = length(g$pi)
  K0 = bestmotif$K
  W0 = bestmotif$W
  clustlike0 = bestmotif$clustlike
  lfdr = matrix(0,n,R)
  lfsr = matrix(0,n,R)
  post_mean = matrix(0,n,R)
  Theta = array(0,dim=c(R,n,L))
  for(r in 1:R){
    we = array(0,dim=c(K0,n,L))
    for(k in 1:K0){
      g$pi = W0[[r]][k,]
      we[k,,]= t(comp_postprob(g,data[[r]]))
    }
    for(j in 1:n){
      if(K0>1) Theta[r,j,] = t(we[,j,])%*%clustlike0[j,]
      else Theta[r,j,] = we[,j,]
    }
    lfdr[,r] = Theta[r,,1]
    NegativeProb = rowSums(t(comp_cdf_post(g,0,data[[r]]))*Theta[r,,])-lfdr[,r]
    lfsr[,r] = compute_lfsr(NegativeProb,lfdr[,r])
    post_mean[,r] = rowSums(t(comp_postmean(g,data[[r]]))*Theta[r,,])
  }
  
  list(lfdr=lfdr,lfsr=lfsr,post_mean = post_mean,Theta=Theta)
}


#results = generic.cormotif(betahat,sebetahat,K=1:4)

