r_tot=function(data){
  r_temp=NULL
  data=as.matrix(data)
  for (i in 1:(length(data[1,])-1)){
    for (j in (i+1):(length(data[1,]))){
      r_temp=c(r_temp,cor(data[,c(i,j)])[1,2])
    }
  }
  y=mean(r_temp)
  return(y)
}
############################boootstrap CIs#########################################################
bootstrap_CIs=function(data,fxn=r_tot,fxn_dim=1,B=1000,alfa=0.05){
  ##euresh toy gnwstou 'r'
  
  if (fxn_dim==1){ #monometablhth analush
    
    r=fxn(data)
    n=length(data)
    set.seed(1)
    y=replicate(B,data[sample(1:n,n,replace = TRUE)])
    
    #euresh tou bootstrap CI
    bootstrap.cint=NULL
    r.star=NULL
    T.star=NULL
    for (b in 1:B){
      r.star[b]=fxn(y[,b])
    }
    bootstrap.cint[1]=unname(2*r-rev(quantile(r.star,c(alfa/2,1-alfa/2),type=1)))[1]
    bootstrap.cint[2]=unname(2*r-rev(quantile(r.star,c(alfa/2,1-alfa/2),type=1)))[2]
    
    #euresh bootstrap-t
    #euresh sd sto kathena bootstrap
    
    bootstrap_t.cint=NULL
    for (b in 1:B){
      
      r.star2=NULL
      for (k in 1:100) {
        l=sample(1:n,n,replace=TRUE)
        r.star2[k]=fxn(y[,b][l])
      }
      T.star[b]=(r.star[b]-r)/sd(r.star2)
    }
    bootstrap_t.cint[1]=r-quantile(T.star,c(1-alfa/2,alfa/2),type=1,names=F)[1]*sd(r.star)
    bootstrap_t.cint[2]=r-quantile(T.star,c(1-alfa/2,alfa/2),type=1,names=F)[2]*sd(r.star)
    
    #euresh quantile bootstrap
    quant_bootstrap.cint=NULL
    quant_bootstrap.cint[1]=quantile(r.star,c(alfa/2,1-alfa/2),type=1,names=F)[1]
    quant_bootstrap.cint[2]=quantile(r.star,c(alfa/2,1-alfa/2),type=1,names=F)[2]
    
    #euresh bca bootstrap
    r.jack = NULL
    bca_bootstrap.cint=NULL
    for ( k in 1:n ) r.jack[k]=fxn(y[-k,b])
    a = sum((r.jack-mean(r.jack))^3)/(6*sum((r.jack-mean(r.jack))^2)^(3/2))
    z0 = qnorm(mean(r.star<r))
    z = qnorm(0.975)
    alfa1 = pnorm(z0+(z0-z)/(1-alfa*(z0-z)))
    alfa2 = pnorm(z0+(z0+z)/(1-a*(z0+z)))
    bca_bootstrap.cint[1]=quantile(r.star,c(alfa1,alfa2),type=1,names=F)[1]
    bca_bootstrap.cint[2]=quantile(r.star,c(alfa1,alfa2),type=1,names=F)[2]
    pinakas=matrix(c(bootstrap.cint[1],bootstrap_t.cint[1],quant_bootstrap.cint[1],bca_bootstrap.cint[1],bootstrap.cint[2],bootstrap_t.cint[2],quant_bootstrap.cint[2],bca_bootstrap.cint[2]),ncol = 2)
    row.names(pinakas)=c('Bootstrap','Bootstrap-t','Quantile bootstrap','Bca bootstrap')
    colnames(pinakas)=c('Lower','Upper')
    list(est=mean(r.star),se=sd(r.star),CIs=pinakas)
  }
  else {   #fxn_dim!=1 dhladh polumetablhth analush
    data=as.matrix(data)
    n=length(data[,1])
    r=fxn(data)
    set.seed(1)
    y=replicate(B,data[sample(1:n,n,replace = TRUE),])
    
    #euresh tou bootstrap CI
    bootstrap.cint=NULL
    r.star=NULL
    T.star=NULL
    for (b in 1:B){
      r.star[b]=fxn(y[,,b])
    }
    bootstrap.cint[1]=unname(2*r-rev(quantile(r.star,c(alfa/2,1-alfa/2),type=1)))[1]
    bootstrap.cint[2]=unname(2*r-rev(quantile(r.star,c(alfa/2,1-alfa/2),type=1)))[2]
    
    #euresh bootstrap-t
    #euresh sd sto kathena bootstrap
    
    bootstrap_t.cint=NULL
    for (b in 1:B){
      
      r.star2=NULL
      for (k in 1:100) {
        l=sample(1:n,n,replace=TRUE)
        r.star2[k]=fxn(y[,,b][l,])
      }
      T.star[b]=(r.star[b]-r)/sd(r.star2)
    }
    bootstrap_t.cint[1]=r-quantile(T.star,c(1-alfa/2,alfa/2),type=1,names=F)[1]*sd(r.star)
    bootstrap_t.cint[2]=r-quantile(T.star,c(1-alfa/2,alfa/2),type=1,names=F)[2]*sd(r.star)
    
    #euresh quantile bootstrap
    quant_bootstrap.cint=NULL
    quant_bootstrap.cint[1]=quantile(r.star,c(alfa/2,1-alfa/2),type=1,names=F)[1]
    quant_bootstrap.cint[2]=quantile(r.star,c(alfa/2,1-alfa/2),type=1,names=F)[2]
    
    #euresh bca bootstrap
    r.jack = NULL
    bca_bootstrap.cint=NULL
    for ( k in 1:n ) r.jack[k]=fxn(y[-k,,b])
    a = sum((r.jack-mean(r.jack))^3)/(6*sum((r.jack-mean(r.jack))^2)^(3/2))
    z0 = qnorm(mean(r.star<r))
    z = qnorm(0.975)
    alfa1 = pnorm(z0+(z0-z)/(1-alfa*(z0-z)))
    alfa2 = pnorm(z0+(z0+z)/(1-a*(z0+z)))
    bca_bootstrap.cint[1]=quantile(r.star,c(alfa1,alfa2),type=1,names=F)[1]
    bca_bootstrap.cint[2]=quantile(r.star,c(alfa1,alfa2),type=1,names=F)[2]
    pinakas=matrix(c(bootstrap.cint[1],bootstrap_t.cint[1],quant_bootstrap.cint[1],bca_bootstrap.cint[1],bootstrap.cint[2],bootstrap_t.cint[2],quant_bootstrap.cint[2],bca_bootstrap.cint[2]),ncol = 2)
    row.names(pinakas)=c('Bootstrap','Bootstrap-t','Quantile bootstrap','Bca bootstrap')
    colnames(pinakas)=c('Lower','Upper')
    list(est=mean(r.star),se=sd(r.star),CIs=pinakas)
 }
}


##################################### jack-knife############################################
jknife=function(x,theta_f=sd,alfa=0.02,fxn_dim=1) {
  if (fxn_dim==1) {
    theta=theta_f(x)
    n=length(x)
    partials=NULL
    
    for (i in 1:n) {
      partials[i]=theta_f(x[-i])
    }
    pseudos=(n*theta)-(n-1)*partials
    jack.bias=mean((n-1)*(partials-theta))
    jack.est=mean(pseudos)
    jack.se=sd(pseudos)/sqrt(n)
    CI=qt(alfa/2,n-1,lower.tail = FALSE)*jack.se
    jack.ci=c(jack.est-CI,jack.est+CI)
    list(est=jack.est,se=jack.se,ci=jack.ci,pseudos=pseudos,partials=partials,bias=jack.bias)
    
  }
  else {
    theta=theta_f(x)
    n=nrow(x)
    partials=NULL
    for (i in 1:n) {
      partials[i]=theta_f(x[-i,])
    }
    pseudos=(n*theta)-(n-1)*partials
    jack.bias=mean((n-1)*(partials-theta))
    jack.est=mean(pseudos)
    jack.se=sd(pseudos)/sqrt(n)
    CI=qt(alfa/2,n-1,lower.tail = FALSE)*jack.se
    jack.ci=c(jack.est-CI,jack.est+CI)
    list(est=jack.est,se=jack.se,ci=jack.ci,pseudos=pseudos,partials=partials,bias=jack.bias)
  }
}

mu_test=function(data,B=2000,alfa=0.05,mu0=0){
  data=as.matrix(data)
  x=data[,1];y=data[,2];
  n=length(x)
  m=length(y)
  xy.bar=mean(c(x,y))
  x.tilde = x-mean(x)+xy.bar
  y.tilde = y-mean(y)+xy.bar
  set.seed(1)
  x.star=replicate(B,x.tilde[sample(1:n,n,replace = TRUE)])
  y.star=replicate(B,y.tilde[sample(1:n,n,replace = TRUE)])
  
  #euresh tou bootstrap CI
  bootstrap_cint=NULL
  
  t.star = NULL
  for (b in 1:B){
    t.star[b] = (mean(x.star[,b])-mean(y.star[,b]))/sqrt(var(x.star[,b])/n+var(y.star[,b])/m)
  }
  
  bootstrap_cint[1]=t-mu0-quantile(t.star,c(1-alfa/2,alfa/2),type=1,names=F)[1]*sd(t.star)
  bootstrap_cint[2]=t-mu0-quantile(t.star,c(1-alfa/2,alfa/2),type=1,names=F)[2]*sd(t.star)
  list(est=mean(t.star),p_value=2*min(mean(t.star<=t+mu0),mean(t.star>=t+mu0)),CIs=bootstrap_cint)
}
