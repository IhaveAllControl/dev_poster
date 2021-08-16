library(foreach)
library(doParallel)
cl <- makeCluster(11)
registerDoParallel(cl)
setwd("D:/Down/0115_feedback")
source("D:/Down/0115_feedback/cd_Lasso_C.r")
simul_data <- function(n, p,sd=1){
  X = matrix(rnorm(n*p, 0, 1), n, p)
  #B = c(rep(1, 0.05*p), rep(-1, 0.05*p), rep(2,0.05*p),rep(0, 0.85*p))
  c1_be <- rep(0,p)
  c1_be[1:20] = 2
  c1_be[21:25] = -1
  c1_be[26:50] = 3
  
  B = c1_be
  e = rnorm(n, 0, sd)
  y = X %*% B + e
  
  test_X = matrix(rnorm(n*p, 0, 1), n, p)
  test_e = rnorm(n, 0, sd)
  test_y = test_X %*% B + test_e
  return(list(X = scale(X,T,F), y=scale(y,T,F), testx = scale(test_X,T,F), testy = scale(test_y,T,F)))
  
}
cvcompare_no_bias_v2 = function(n, p, X_data, y_data, X.test, y.test, k, using_data, st=5,len=2,lv=NULL) {
  tol=1e-6;
  cd_iter = 1e6
  niter = 100
  init=rep(0,p)
  
  if(is.null(lv)) lv <- max(abs(t(dat$X)%*%dat$y))
  num.fol <- floor(dim(X_data)[1]/k)
  lam.vec <- seq(from=st, to=lv, by=len)
  nl = length(lam.vec)
  res.matrix <- vector()
  
  for (fold in 1:k) {
    
    if (fold!=k) {
      i.except.train <- (1 + (fold-1) * num.fol):(fold * num.fol)
    }
    else{
      i.except.train <- (1 + (fold-1) * num.fol):dim(X_data)[1]
    }
    
    if(using_data == 'train')
    {
      
      print(paste0("train","_",fold,"/",k))
      fr = foreach(lam = lam.vec, .combine = rbind )%dopar% {
        cd_lasso = function(X,y,lambda,tol,niter,init,loc="D:/Down/0115_feedback")
        { 
          n <- as.integer(nrow(X))
          p <- as.integer(ncol(X))
          initer <- as.integer(niter)
          iriter <- as.integer(0)
          
          dX <- as.double(X)
          dy <- as.double(y)
          dlam <- as.double(lambda)
          dbe <- as.double(rep(0,p))
          dtol <- as.double(tol)
          ddiff <- as.double(0)
          dres <- as.double(rep(0,n))
          dinit <- as.double(init)
          
          dyn.load(paste0(loc,"/cd_lasso",.Platform$dynlib.ext))
          out <- .C("cd_lasso",n=n,p=p,y=dy,X=dX,beta=dbe,lam=dlam,niter=initer,tol=dtol,riter=iriter,diff=ddiff,resid=dres,init=dinit)
          dyn.unload(paste0(loc,"/cd_lasso",.Platform$dynlib.ext))
          
          return(list(est=out$beta,resid=out$resid,max_diff=out$diff,riter=out$riter))
          
        }
        l =  cd_lasso(X_data[-i.except.train,],y_data[-i.except.train],lam,1e-6,1e6,rep(0,p))
        mean((y_data[i.except.train] - X_data[i.except.train,] %*% l$est)^2)
      }
      res.matrix = cbind(res.matrix, fr[,1])
    }
    
  }
  
  rownames(res.matrix) <- lam.vec
  res.matrix <- apply(res.matrix, 1, mean, na.rm=T)
  return(res.matrix)
}


make_iter = function(n, p, k, using_data,
                     st,len,lv, iteration) {
  cand_alpha = seq(0.5,2,by=0.01)
  result = 0
  result2 = vector()
  set.seed(2021)
  for (i in 1:iteration) {
    dat = simul_data(n,p,sd=0.5)
    X = dat$X
    y = dat$y
    tX = dat$testx
    ty = dat$testy
    out = cvcompare_no_bias_v2(n, p, X_data=X, y_data=y, X.test=tX, y.test=ty, 
                               k, using_data='train',st=st,len=len,lv=lv)
    result = cbind(result, out)
    l2 =  cd_lasso(X, y, as.numeric(names(which.min(out))),1e-6,1e6,rep(0,p))
    result2[i] = mean((ty- tX %*% l2$est)^2)
    
    print(paste0(i,"--iteration")) 
  }
  return(list(result[,-1], result2))
}

nm2 = make_iter(200,500,2,"train", st = 0, len =0.4, lv =500, iteration = 50)

nm5 = make_iter(200,500,5,"train", st = 0, len =0.4, lv =500, iteration = 50)
nm10 = make_iter(200,500,10,"train", st = 0, len =0.4, lv =500, iteration = 50)
nm20 = make_iter(200,500,20,"train", st = 0, len =0.4, lv =500, iteration = 50)

nm50 = make_iter(200,500,50,"train", st = 0, len =0.4, lv =500, iteration = 50)
save(list=ls(), file = "C:/Users/pc/OneDrive/¹ÙÅÁ È­¸é/cv/May/lasso200500(until50).RData")




a1 = apply(nm2[[1]], 2, function(x) seq(0,500,0.4)[which.min(x)])
    

    