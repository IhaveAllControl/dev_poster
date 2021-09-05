
cd_enet <- function(X,y,lambda1,lambda2, tol,niter,init,loc=getwd())
{ 
  n <- as.integer(nrow(X))
  p <- as.integer(ncol(X))
  initer <- as.integer(niter)
  iriter <- as.integer(0)
  
  dX <- as.double(X)
  dy <- as.double(y)
  d_lam1 <- as.double(lambda1)
  d_lam2 = as.double(lambda2)
  dbe <- as.double(rep(0,p))
  dtol <- as.double(tol)
  ddiff <- as.double(0)
  dres <- as.double(rep(0,n))
  dinit <- as.double(init)
  
  dyn.load(paste0(loc,"/cd_enet",.Platform$dynlib.ext))
  out <- .C("cd_enet_ftn",n=n,p=p,y=dy,X=dX,beta=dbe,lam1=d_lam1,lam2=d_lam2, niter=initer,tol=dtol,riter=iriter,diff=ddiff,resid=dres,init=dinit)
  dyn.unload(paste0(loc,"/cd_enet",.Platform$dynlib.ext))
   
  return(list(est=out$beta,resid=out$resid,max_diff=out$diff,riter=out$riter))
  
}

cal_enet_debias <- function(X,y,est,lambda1,lambda2)
{
  n = nrow(X)
  p = ncol(X)
  
  indx = which(est!=0)
  if(length(indx)==0)
     return(list(bias=rep(0,p),index=1:p,length=p))
  len = length(indx)
  if(len==1)
  {
     XI = matrix(X[,indx],n,len,byrow=F)
  
   } else {
       XI = X[,indx]
   }
  lam2 = lambda2
  
  XItY = (t(XI)%*%y)
  XItXI = t(XI)%*%XI
  est_c1 = est*(1+lam2)
  if(len < n)
  {
      est_c2 = rep(0,p)
      est_c2_temp = solve(XItXI)%*%XItY
      est_c2[indx] = est_c2_temp
  } else {
      est_c2 = NULL
  }
  
  if(len < n | lam2 >0)
  {
      est_c3 = rep(0,p)
      est_c3_temp = solve(XItXI+lam2*diag(len))%*%XItY
      est_c3[indx] = est_c3_temp
      est_c4 = est_c3*(1+lam2)
  } else {
      est_c3 = NULL
      est_c4 = NULL
  }
  
  
  if(lambda2==0)
    df = len
  if(lambda2>0)
    df = sum(diag(XI%*%solve(XItXI+lambda2*diag(len))%*%t(XI)))
  

  return(list(est_c0=est, est_c1=est_c1, est_c2=est_c2, est_c3=est_c3, est_c4=est_c4, df=df))
}


info_crit  = function(X, y, est, df, lambda1,lambda2)
{
  n = nrow(X)
  p = ncol(X)
  
  if(is.null(est))
  {
    return(c(Inf,Inf,Inf))
  }
    
  s = sum(est!=0)
  RSS = sum((y - X%*%est)^2)
  nRSS = RSS/n
  
  crit_base = n*log(nRSS)
  aic = crit_base + df*2
  bic = crit_base + df * log(n)
  ebic = bic + lchoose(p,s)
  
  res = c(aic,bic,ebic)
  names(res) = c("AIC","BIC","EBIC")

  return(res)
}
